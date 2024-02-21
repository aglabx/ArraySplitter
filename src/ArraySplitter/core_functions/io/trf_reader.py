#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Classes:

- TRFFileIO(AbstractBlockFileIO)

"""

from .block_file import AbstractBlockFileIO
from ..models.trf_model import TRModel


class TRFFileIO(AbstractBlockFileIO):
    """Working with raw output from TRF."""

    def __init__(self):
        """Overridden. Hardcoded start token."""
        token = "Sequence:"
        super(TRFFileIO, self).__init__(token)

    def iter_parse(self, trf_file, filter=True):
        """Iterate over raw trf data and yield TRFObjs.

        Args:
            trf_file (str): The path to the TRF file.
            filter (bool, optional): Whether to filter the TRF objects. Defaults to True.

        Yields:
            TRModel: The parsed TRF object.

        """
        for ii, (head, body, start, next) in enumerate(self.read_online(trf_file)):
            head = head.replace("\t", " ")
            obj_set = []
            n = body.count("\n")
            for i, line in enumerate(self._gen_data_line(body)):
                trf_obj = TRModel()
                trf_obj.set_raw_trf(head, None, line)
                obj_set.append(trf_obj)
            yield from obj_set

    def _gen_data_line(self, data):
        """Generate lines of data from the raw TRF data.

        Args:
            data (str): The raw TRF data.

        Yields:
            str: A line of data.

        """
        for line in data.split("\n"):
            line = line.strip()
            if line.startswith("Sequence"):
                continue
            if line.startswith("Parameters"):
                continue
            if not line:
                continue
            yield line


def sc_iter_arrays_trf_file(file_name):
    """Iter over trf file."""
    reader = TRFFileIO()
    for trf_obj in reader.iter_parse(file_name):
        yield trf_obj.trf_array

def sc_iter_trf_file(file_name):
    """Iter over trf file."""
    reader = TRFFileIO()
    for trf_obj in reader.iter_parse(file_name):
        yield trf_obj.trf_id, trf_obj.trf_array
