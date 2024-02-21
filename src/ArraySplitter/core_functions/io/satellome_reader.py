#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from .tab_file import sc_iter_tab_file
from ..models.trf_model import TRModel


def sc_iter_arrays_satellome_file(trf_file):
    """Iter over satellome file."""
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        yield trf_obj.trf_array

def sc_iter_satellome_file(trf_file):
    """Iter over satellome file."""
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        yield trf_obj.trf_id, trf_obj.trf_array
