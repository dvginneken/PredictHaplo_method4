
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


"""This module exports functions that help interact with
PredictHaplo-Paired and is part of the HaploHIV pipeline. The exported
functions are write_conf(...), which writes a configuration file for
PredictHaplo-Paired and also checks whether the files that it refers
to exist; greatest_recon_window(...), which selects from the passed
list of paths to PredictHaplo-Paired output files the path to the file
that encompasses the greatest reconstruction window; and
fas_records(...), which produces an iterator over the records in these
files.

"""


import dataclasses
import enum
import logging
import os
import re

from implementation import util


_SELF_DIRNAME = util.absdirname(__file__)
_PH_CONF_TEMPLATE_PATH = os.path.join(_SELF_DIRNAME,
                                      "templates",
                                      "predicthaplo.conf.template")

_WINDOW_RE = re.compile(r"global_(?P<start>\d+)_(?P<end>\d+)\.fas$")

_log = logging.getLogger(__name__)


class _ParserState(enum.Enum):
    """States that a PredictHaploFasParser can have."""
    START_FILE = 0
    NEW_RECORD = 1
    COMMENT = 2
    CONFIDENCE = 3
    SEQUENCE = 4


class _PredictHaploFasParser:
    """Parses FASTA files output by PredictHaplo."""
    _END_OF_COMMENT = ";EndOfComments"
    _START_OF_CONFIDENCE = ";Confidence scores"
    _OVERLAP = ";Overlap quality scores "
    _FREQ_RE = re.compile(r";Freq:(?P<frequency>\d+(?:\.\d+)?)")

    def __init__(self, file_in):
        """Initialize also current_record and state, which are both stateful."""
        self._record = _PredictHaploFasRecord()
        self._file = file_in
        self._state = _ParserState.START_FILE
        self._maybe_finished = False

    def _consume_comment(self):
        """Advance file and update records until end of comment is seen."""
        # caller's responsibility to make current line be start of record
        assert self._state == _ParserState.NEW_RECORD
        self._state = _ParserState.COMMENT

        legal_states = _ParserState.COMMENT, _ParserState.CONFIDENCE
        for line in self._file:
            assert self._state in legal_states
            assert line[0] == ";"

            if line.startswith(self._END_OF_COMMENT):
                self._state = _ParserState.SEQUENCE
                break

            elif line.startswith(self._START_OF_CONFIDENCE):
                # we don't need confidence string and null byte is only easy
                # terminator for confidence string, so ignore for now
                self._update_record(confidence_score=True)
                self._state = _ParserState.CONFIDENCE

            elif line.startswith(self._OVERLAP):
                # no clue how to represent overlap quality score
                qual = line[len(self._OVERLAP):]
                self._update_record(overlap_quality_score=qual)

            else:
                match = self._FREQ_RE.match(line)
                if match is not None:
                    self._state = _ParserState.COMMENT
                    frequency = float(match.group("frequency"))
                    self._record = self._record.update(frequency=frequency)
#HACK
#                else:
#                    # we are ignoring confidence strings
#                    assert self._state == _ParserState.CONFIDENCE

        # check whether the loop exited prematurely
        if not self._state == _ParserState.SEQUENCE:
            # since before we could not check for end of file,
            # this has to do instead
            if self._maybe_finished and self._state == _ParserState.COMMENT:
                return False

            # only way to get here is if file is malformed
            assert False
        return True

    def _consume_sequence(self):
        """Advance file until end of sequence, updating record."""
        assert self._state == _ParserState.SEQUENCE
        sequence = ""

        for line in self._file:
            if line[0] == ">":
                break
            sequence += line.rstrip()

        # EOF also marks end of sequence
        self._state = _ParserState.NEW_RECORD
        self._update_record(sequence=sequence)

    def _update_record(self, **kwargs):
        """Update the record, which is a frozen dataclass."""
        self._record = self._record.update(**kwargs)

    def is_record_left(self):
        """Advance lines of file until start of record is found."""
        if self._state == _ParserState.START_FILE:
            for line in self._file:
                if line[0] == ">":
                    self._state = _ParserState.NEW_RECORD
                    break

            if self._state == _ParserState.START_FILE:
                # no new record found
                return False

        elif self._state == _ParserState.NEW_RECORD:
            # once next() is called on file, there is no way to know where
            # the cursor is in the file
            self._maybe_finished = True

        else:
            # parser is in illegal state
            assert False

        return True

    def next_record(self):
        """Parse and advance file until next record is complete."""
        legal_states = _ParserState.START_FILE, _ParserState.NEW_RECORD
        assert self._state in legal_states

        if not self.is_record_left():
            return None

        if not self._consume_comment():
            return None

        self._consume_sequence()
        assert self._record.is_complete()
        return self._record


@dataclasses.dataclass(frozen=True)
class _PredictHaploFasRecord:
    frequency: float = dataclasses.field(default=None)
    overlap_quality_score: str = dataclasses.field(default=None)
    confidence_score: bool = dataclasses.field(default=None)
    sequence: str = dataclasses.field(default=None)

    def update(self, **kwargs):
        """Return new record with fields updated as per keywords."""
        self_dict = dataclasses.asdict(self)
        self_ctor = type(self)

        return self_ctor(**{**self_dict, **kwargs})

    def is_complete(self):
        """Consider a record complete when all of its fields are non-None."""
        self_dict = dataclasses.asdict(self)
        self_values = self_dict.values()

        return all(x is not None for x in self_values)


def _check_file(string):
    """Log and raise if the argument is not a file."""
    if not os.path.isfile(string):
        _log.warning("Does not exist or not a file: %s", string)

        raise FileNotFoundError(string)


def _reconstruction_window(global_fas_path):
    """Compute reconstruction window length from filename."""
    match = _WINDOW_RE.search(global_fas_path)

    if match is None:
        _log.warning("Filename does not match PredictHaplo output: %s",
                     global_fas_path)
        raise ValueError(f"Cannot match filename: {global_fas_path}")

    start, end = [int(match.group(x)) for x in ["start", "end"]]
    return end - start


def fas_records(global_fas_file):
    """Return an iterator over records in a sample_global_start_end.fas file."""
    parser = _PredictHaploFasParser(global_fas_file)
    while True:
        maybe_record = parser.next_record()

        if maybe_record is None:
            break
        yield maybe_record


def greatest_recon_window(global_fas_paths):
    """Select sample_global_start_end.fas file with greatest window."""
    return max(global_fas_paths, key=_reconstruction_window)


def write_conf(*, f_out, sample_name, consensus_path, aligned_reads_path):
    """Write a PredictHaplo configuration to a file object."""
    replacement_dict = locals()
    del replacement_dict["f_out"]

    for path in consensus_path, aligned_reads_path:
        _check_file(path)

    with open(_PH_CONF_TEMPLATE_PATH, "r") as f_in:
        util.write_formatted_template(f_out, f_in, replacement_dict)
