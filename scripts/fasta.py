
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


"""This module is part of the HaploHIV pipeline and exports the FastaRecord
class, which represents records in Fasta files, and the read_records(...) and
write_records(...) functions, that respectively read and write collections of
Fasta records to files.

"""


import dataclasses

from implementation import util
import iupac


def _match_case(case_leader, case_follower):
    """Output a string with the contents of follower, but the case of leader."""
    return "".join(follow.upper() if lead.isupper() else follow.lower()
                   for lead, follow in zip(case_leader, case_follower))


@dataclasses.dataclass
class FastaRecord:
    """Represents a Fasta record, with sequence identifier and sequence."""
    __slots__ = "identifier", "sequence"

    identifier: str
    sequence: str

    def __repr__(self):
        """Output identifier: assume that it is unique."""
        return self.identifier

    def __getitem__(self, index):
        """Ensure easy slicing."""
        return FastaRecord(identifier=self.identifier,
                           sequence=self.sequence[index])

    def __len__(self):
        return len(self.sequence)

    def is_complete(self):
        """Are all fields set to a non-empty value?"""
        field_values = dataclasses.asdict(self).values()
        return all([bool(x) for x in field_values])

    def revcomp(self):
        """Reverse complement of this record."""
        rev = reversed(self.sequence)
        revcomped = "".join(iupac.COMPLEMENT[x.lower()] for x in rev)
        matched = _match_case(revcomped, rev)
        rv_identifier = self.identifier + "_reversecomplement"

        return FastaRecord(identifier=rv_identifier, sequence=matched)

    def show(self, chunk_size=80):
        """FASTA text representation of this record."""
        if not chunk_size or chunk_size <= 0:
            return f">{self.identifier}\n{self.sequence}"

        chunked_sequence = "\n".join(util.chunks(self.sequence, chunk_size))
        return f">{self.identifier}\n{chunked_sequence}\n"

    def with_deletion(self, region):
        """Return this record's sequence with some region removed."""
        start_region, end_region = region
        starting_slice = slice(0, start_region)
        ending_slice = slice(end_region, len(self))

        starting_seq = self.sequence[starting_slice]
        ending_seq = self.sequence[ending_slice]
        sliced_seq = "".join([starting_seq, ending_seq])

        return sliced_seq

    def with_insertion(self, position, insertion):
        """Return this record's sequence with an insertion before position."""
        seq = self.sequence
        inserted_seq = "".join([seq[:position], insertion, seq[position:]])

        return inserted_seq


def read_records(f_in):
    """Yield FastA records from file (should iterate over lines)."""
    current_record = FastaRecord(identifier="", sequence="")
    no_record_yet = True

    for line in f_in:
        stripped = line.rstrip()

        if stripped[0] == ">":
            if current_record.is_complete():
                yield current_record
                no_record_yet = False
                current_record = FastaRecord(identifier="", sequence="")

            elif not no_record_yet:
                # this path is unreachable in valid FASTA files.
                assert False

            current_record.identifier = stripped[1:]
        else:
            current_record.sequence += stripped

    if current_record.is_complete():
        yield current_record


def write_records(f_out, records):
    """Write fasta records to file object."""
    # do not ensure_iterator here, as it will generate
    # one record per base.
    for record in records:
        representation = record.show()
        f_out.write(representation)
