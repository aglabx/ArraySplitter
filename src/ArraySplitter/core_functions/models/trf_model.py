#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import re

from PyExp import AbstractModel

from ..tools.parsers import (parse_chromosome_name,
                                                        parse_fasta_head,
                                                        trf_parse_head,
                                                        trf_parse_line)


def get_gc(sequence):
    """Count GC content."""
    length = len(sequence)
    if not length:
        return 0
    count_c = sequence.count("c") + sequence.count("C")
    count_g = sequence.count("g") + sequence.count("G")
    gc = float(count_c + count_g) / float(length)
    return float(gc)


def clear_sequence(sequence):
    """Clear sequence (full alphabet):

    - lower case
    - \s -> ""
    - [^actgn] -> ""
    """
    sequence = sequence.strip().upper()
    sequence = re.sub("\s+", "", sequence)
    return re.sub("[^actgnuwsmkrybdhvACTGNUWSMKRYBDHV\-]", "", sequence)


class TRModel(AbstractModel):
    """Class for tandem repeat wrapping

    Public properties:

    Indexes:

    - id,
    - giid,
    - trf_id,
    - project,
    - trf_gi

    Coordinates:

    - trf_l_ind,
    - trf_r_indel

    TR:

    - trf_period,
    - trf_n_copy,
    - trf_pmatch,
    - trf_pvar,
    - trf_entropy,
    - trf_array_length,
    - trf_joined,
    - trf_chr

    Sequence:

    - trf_consensus,
    - trf_array

    GC%:

    - trf_array_gc,
    - trf_consensus_gc

    Other:

    - trf_head,
    - trf_params

    Annotation:

    - trf_family,
    - trf_subfamily,
    - trf_subsubfamily,
    - trf_hor,
    - trf_n_chrun,
    - trf_n_refgenome,
    - trf_ref_annotation,
    - trf_bands_refgenome,
    - trf_repbase,
    - trf_strand

    Dumpable attributes:

    - "project",
    - "id" (int),
    - "trf_id" (int),
    - "trf_type",
    - "trf_family",
    - "trf_family_prob",
    - "trf_l_ind" (int),
    - "trf_r_ind" (int),
    - "trf_period" (int),
    - "trf_n_copy" (float),
    - "trf_pmatch" (float),
    - "trf_pvar" (float),
    - "trf_consensus",
    - "trf_array",
    - "trf_array_gc" (float),
    - "trf_consensus_gc" (float),
    - "trf_gi",
    - "trf_head",
    - "trf_param",
    - "trf_array_length" (int),
    - "trf_chr",
    - "trf_joined" (int),
    - "trf_superfamily",
    - "trf_superfamily_ref",
    - "trf_superfamily_self",
    - "trf_subfamily",
    - "trf_subsubfamily",
    - "trf_family_network",
    - "trf_family_self",
    - "trf_family_ref",
    - "trf_hor" (int),
    - "trf_n_chrun" (int),
    - "trf_ref_annotation",
    - "trf_bands_refgenome",
    - "trf_repbase",
    - "trf_strand",

    """

    dumpable_attributes = [
        "project",
        "id",
        "trf_id",
        "trf_type",
        "trf_family",
        "trf_family_prob",
        "trf_l_ind",
        "trf_r_ind",
        "trf_period",
        "trf_n_copy",
        "trf_pmatch",
        "trf_pvar",
        "trf_entropy",
        "trf_consensus",
        "trf_array",
        "trf_array_gc",
        "trf_consensus_gc",
        "trf_gi",
        "trf_head",
        "trf_param",
        "trf_array_length",
        "trf_chr",
        "trf_joined",
        "trf_superfamily",
        "trf_superfamily_ref",
        "trf_superfamily_self",
        "trf_subfamily",
        "trf_subsubfamily",
        "trf_family_network",
        "trf_family_self",
        "trf_family_ref",
        "trf_hor",
        "trf_n_chrun",
        "trf_ref_annotation",
        "trf_bands_refgenome",
        "trf_repbase",
        "trf_strand",
    ]

    int_attributes = [
        "trf_l_ind",
        "trf_r_ind",
        "trf_period",
        "trf_array_length",
        "trf_joined",
        "trf_hor",
        "trf_n_chrun",
    ]

    float_attributes = [
        "trf_n_copy",
        "trf_family_prob",
        "trf_entropy",
        "trf_pmatch",
        "trf_pvar",
        "trf_array_gc",
        "trf_consensus_gc",
    ]

    def set_project_data(self, project):
        """Add project data to self.project."""
        self.project = project

    def set_raw_trf(self, head, body, line):
        """Init object with data from parsed trf ouput.

        Parameters:

        - trf_obj: TRFObj instance
        - head: parsed trf head
        - body: parsed trf body
        - line: parsed trf line
        """
        self.trf_param = 0
        self.trf_head = trf_parse_head(head).strip()
        self.trf_gi = parse_fasta_head(self.trf_head)[0]
        self.trf_chr = parse_chromosome_name(self.trf_head)

        (
            self.trf_l_ind,
            self.trf_r_ind,
            self.trf_period,
            self.trf_n_copy,
            self.trf_l_cons,
            self.trf_pmatch,
            self.trf_indels,
            self.trf_score,
            self.trf_n_a,
            self.trf_n_c,
            self.trf_n_g,
            self.trf_n_t,
            self.trf_entropy,
            self.trf_consensus,
            self.trf_array,
        ) = trf_parse_line(line)

        self.trf_pmatch = float(self.trf_pmatch)
        self.trf_pvar = int(100 - float(self.trf_pmatch))

        try:
            self.trf_l_ind = int(self.trf_l_ind)
        except:
            print("Error:", self)

        self.trf_r_ind = int(self.trf_r_ind)
        self.trf_period = int(self.trf_period)
        self.trf_n_copy = float(self.trf_n_copy)

        self.trf_consensus = clear_sequence(self.trf_consensus)
        self.trf_array = clear_sequence(self.trf_array)

        self.trf_array_gc = get_gc(self.trf_array)
        self.trf_consensus_gc = get_gc(self.trf_consensus)
        self.trf_chr = parse_chromosome_name(self.trf_head)
        self.trf_array_length = len(self.trf_array)

    def get_header_string(self):
        """Get header string for tsv file."""
        data = "\t".join(self.dumpable_attributes)
        return f"#{data}\n"

    def get_fasta_repr(self, add_project=False):
        """Get array fasta representation, head - trf_id."""
        if add_project:
            return ">%s_%s\n%s\n" % (self.trf_id, self.project, self.trf_array)
        else:
            return ">%s\n%s\n" % (self.trf_id, self.trf_array)

    @property
    def fasta(self):
        return self.get_fasta_repr()

    def get_gff3_string(
        self,
        chromosome=True,
        trs_type="complex_tandem_repeat",
        probability=1000,
        tool="PySatDNA",
        prefix=None,
        properties=None,
        force_header=False,
    ):
        """Return TR in gff format."""
        if chromosome and self.trf_chr and self.trf_chr != "?":
            seqid = self.trf_chr
        elif self.trf_gi and self.trf_gi != "Unknown":
            seqid = self.trf_gi
        else:
            seqid = self.trf_head
        features = []
        if not properties:
            properties = {}
        for name, attr in properties.items():
            features.append("%s=%s" % (name, getattr(self, attr)))
        features = ";".join(features)
        if prefix:
            seqid = prefix + seqid
        if self.trf_l_ind < self.trf_r_ind:
            strand = "+"
        else:
            strand = "-"
            self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

        if force_header:
            seqid = self.trf_head
        d = (
            seqid,
            tool,
            trs_type,
            self.trf_l_ind,
            self.trf_r_ind,
            probability,
            strand,
            ".",
            features,
        )
        return "%s\n" % "\t".join(map(str, d))

    def get_bed_string(self):
        """Return TR in bed format."""
        if self.trf_l_ind < self.trf_r_ind:
            strand = "+"
        else:
            strand = "-"
            self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

        seqid = self.trf_head
        d = (
            seqid,
            self.trf_l_ind,
            self.trf_r_ind,
        )
        return "%s\n" % "\t".join(map(str, d))
