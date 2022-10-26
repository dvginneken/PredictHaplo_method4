
# HaploHIV pipeline for retrieving haplotypes from paired-end reads from HIV
# Copyright (C) 2019  Gijs Schroeder, Terry Huisman, Kobus Bosman,
#                     Monique Nijhuis, Rob de Boer, Aridaman Pandit
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#
# E-mail us at these addresses:
#
# gijsschroder@gmail.com    <Gijs Schroeder>
# terry@mhuisman.cumail.nl  <Terry Huisman>
# A.Pandit@umcutrecht.nl    <Aridaman Pandit>


"""This module exports a function that finds the best supported
scaffold to map reads to using MIRA as part of the HaploHIV pipeline.

"""


import logging
import operator
import sys

from externals import samutils


_log = logging.getLogger(__name__)

# column 0 in idxstats TSV file contains name, column 2 number of mapped reads
_COL_NAME = 0
_COL_NREADS = 2


def choose_scaffold(*, reference_alignment, ec_reads, min_ratio):
    """Return name of scaffold with most support or None if little support."""
    try:
        nr_ec_reads = number_ec_reads(ec_reads)
    except AssertionError:
        _log.error("Cannot read FASTQ records in %s", ec_reads)
        sys.exit(1)

    name, n_reads = best_supported(reference_alignment)

    if n_reads / nr_ec_reads < min_ratio:
        _log.critical("Best supported read %s has %d reads out of"
                      " a total of %d, which is fewer than the"
                      " threshold proportion %f. This threshold"
                      " proportion is set via the option 'min_mapping_ratio'.",
                      name, n_reads, nr_ec_reads, min_ratio)
        result = None
    else:
        result = name

    return result


def best_supported(reference_alignment):
    """To what reference do most reads map and how many do so?"""
    idxstats_table = samutils.idxstats(reference_alignment)
    best_record = max(idxstats_table, key=operator.itemgetter(_COL_NREADS))

    name, n_reads = [best_record[index] for index in [_COL_NAME, _COL_NREADS]]
    return name, n_reads


def number_ec_reads(ec_reads):
    """Count number of error corrected reads."""
    with open(ec_reads, "r") as f_in:
        n_records = 0

        for index, line in enumerate(f_in):
            line_nr_mod4 = index % 4

            if line_nr_mod4 == 0:
                assert line[0] == "@"
                n_records += 1

            elif line_nr_mod4 == 2:
                assert line[0] == "+"

    return n_records
