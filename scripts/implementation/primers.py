
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


"""This module is part of the HaploHIV pipeline and exports two functions that
are used for cutting primers from reference or consensus sequences:

 - primer_loci(...)
   For each passed sequence, yield the best primer location given a collection
   of primers and what end we are cutting from -- are these primers head primers
   or tail primers? We assume that tail primers are reverse complemented, such
   their sequence aligns to its proper location.

 - cut_primers(...)
   Given (sequence, primer location) pairs and the end that we are cutting from,
   yield sequences from which the primers are cut.

Both the sequences that are cut and the primers themselves may contain IUPAC
ambiguous nucleotide codes. These codes are interpreted differently depending on
whether they are encountered in a sequence or in a primer.

In a sequence a IUPAC ambiguous nucleotide code indicates uncertainty as to what
nucleotide should be represented in that particular location. In contrast, in a
primer a IUPAC ambiguous nucleotide code indicates certainty that for each
ambiguously identified nucleotide, there is a primer in the mixture that has
that nucleotide in that position. For example, the following fasta file
containing primers:

>collection_of_primers
atgckatgc
>another_primer
cgtacgta

would indicate that all of the atgcgatgc, atgctatgc, and cgtacgta primers are
present in the mixture.

This difference in meaning between IUPAC ambiguous nucleotides in sequences and
primers has consequences for aligning. Since ambiguity in primers indicate
certainty as to the presence of certain nucleotide, all of the following align-
ments between primer and sequence score equally:

sequence atgcgatgc
primer   atgckatgc

sequence atgctatgc
primer   atgckatgc

sequence atgcgatgc
primer   atgcgatgc

sequence atgctatgc
primer   atgctatgc

On the other hand, the alignment

sequence atgcgatgc
primer   atgcgatgc

has a better score than the alignment

sequence atgckatgc
primer   atgcgatgc

The similarity formula for each IUPAC code pair is that we take the number of
all nucleotide possibilities in the sequence that are covered by any possibility
in the primer and divide that by the number of nucleotide possibilties in the
sequence (see _iupac_match_score). Example:

sequence k has as possible nucleotides g and t.
primer g has as only possibility g.

#Sequence nucleotides possible: 2 (g is one, t is another)

Sequence nucleotides covered by primer: g is covered, t is not. Therefore,
#Sequence nucleotides covered = 1

The formula is: #Sequence nucleotides covered / #Sequence nucleotides possible
Substitute values: 1 / 2
"""


import dataclasses
import enum
import functools
import itertools
import logging
import re

from Bio import pairwise2

from sequences import fasta
from sequences import iupac
from . import exceptions


_log = logging.getLogger(__name__)


# These parameters are the same as MAFFT --localpair defaults.
# Experimentation showed that these values work better than alternatives.
_GAP_OPEN_PENALTY = -2.0
_GAP_EXT_PENALTY = -0.1

_IUPAC_TO_RE = {
    iupaccode: "[" + "".join(nucleotides) + "]" \
               if len(nucleotides) > 1
               else nucleotides[0]
    for iupaccode, nucleotides in iupac.IUPAC_TO_NUCS.items()
}
_IUPAC_CHARS = list(iupac.IUPAC_TO_NUCS.keys())
_IUPAC_ALPHA = "".join(_IUPAC_CHARS)

_NONGAP_RE = re.compile(r"[^-]")
_NONIUPAC_RE = re.compile("[^" + _IUPAC_ALPHA + "]")


def _iupac_match_score(iupac_primer, iupac_seq):
    """Ratio of how many options in sequence are covered by primer."""
    # Ambiguous nucleotide codes in primer indicate certainty of different
    # primers in mixture, whereas ambiguous nucleotide codes in sequence
    # indicate uncertainty as to what nucleotide is there.

    # Therefore, treat uncertainty as uniformly distributed random variate
    # and score for expected match
    nucs_p, nucs_s = [iupac.IUPAC_TO_NUCS[c] for c in [iupac_primer, iupac_seq]]
    n_seq_nucs_covered = sum(1 if nuc in nucs_p else 0 for nuc in nucs_s)

    return n_seq_nucs_covered / len(nucs_s)


# IUPAC similarity matrix, where primer goes left and sequence goes right.
_IUPAC_SIM = {
    (iupac_primer, iupac_seq): _iupac_match_score(iupac_primer, iupac_seq)
    for iupac_primer, iupac_seq in itertools.product(_IUPAC_CHARS, _IUPAC_CHARS)
}


@dataclasses.dataclass(order=True)
class _Alignment:
    """Representation of an alignment between two sequences."""
    score: float
    _span: int
    start: int = dataclasses.field(compare=False)
    seq_a: str = dataclasses.field(compare=False)
    seq_b: str = dataclasses.field(compare=False)

    def __post_init__(self):
        """Negate span so that Alignments order correctly."""
        self._span = -self._span

    def __repr__(self):
        """String representation of alignment."""
        offset = self.start - self.starting_pos()
        positions = [p for p in range(offset, offset+len(self.seq_a))]
        ticks = [t for t in positions if (t%10) == 0]
        indent = 4
        out = "\n"
        out += indent*" " + self.seq_a + "\n"
        out += indent*" " + self.seq_b + "\n"

        out += indent*" "
        out += "".join(["|" if i in ticks else " " for i in positions])
        out += "\n"

        pos = offset - indent - 1
        for tick in ticks:
            tick_s = str(tick)
            width = len(tick_s)
            nspaces = tick - pos - width
            out += nspaces*" " + tick_s
            pos = tick
        out += "\n"

        return out

    @property
    def span(self):
        """Unnegate span for the outside world."""
        return -self._span

    def starting_pos(self):
        """Find the position in the aligned sequences corresponding to start."""
        match_a, match_b = [
            _NONGAP_RE.search(s) for s in [self.seq_a, self.seq_b]
        ]
        # neither should be None
        assert match_a and match_b
        return max(match_a.start(), match_b.start())


class End(enum.Enum):
    """Which end are we processing?"""
    HEAD = 0
    TAIL = 1


def _sliding_window(sequence, window_size):
    """Yield each window of fixed size from sequence."""
    lseq = len(sequence)
    assert window_size <= lseq

    window_starts = range(lseq - window_size)
    window_ends = range(window_size, lseq)

    for start, end in zip(window_starts, window_ends):
        yield sequence[start:end]


def _kmer_to_regex(k_mer):
    """Translate the K-mer to a regex."""
    regex_s = "".join(_IUPAC_TO_RE[iupaccode] for iupaccode in k_mer.lower())

    return re.compile(regex_s, flags=re.IGNORECASE)


def _primer_to_regexes(primer, k_mer_len):
    """Return the regexes that find each K-mer in the primer."""
    return (_kmer_to_regex(k_mer)
            for k_mer in _sliding_window(primer, k_mer_len))


def _search_each(regex, haystack):
    """Regex search repeatedly to find all matches of regex in haystack."""
    haystack_remaining = haystack

    while True:
        match = regex.search(haystack_remaining)
        if not match:
            break

        span = match.span()
        yield span

        haystack_remaining = haystack_remaining[1 + span[0]:]


def _exact_kmer_matches(primer, bio_seq, k_mer_len):
    """Yield inferred primer locations from exact K-mer matches in sequence."""
    locations = set()
    for index, regex in enumerate(_primer_to_regexes(primer, k_mer_len)):
        for start, _ in _search_each(regex, bio_seq):
            matched_primer_location = start - index

            aln = _Alignment(0, k_mer_len, start,
                             primer[index:k_mer_len + index],
                             bio_seq[start:k_mer_len + start])
            _log.debug("Exact K-mer match with primer %r:\n%r", primer, aln)

            if matched_primer_location not in locations:
                locations.add(matched_primer_location)
                yield matched_primer_location


def _tup_to_aln(offset, tup):
    """Convert a tuple from pairwise2 alignment to an Alignment object."""
    seq_a, seq_b, score, start, end = tup
    span = end - start

    return _Alignment(score, span, offset + start, seq_a, seq_b)


def _min_overlap_alignments(primer, bio_seq, min_overlap):
    """Return all alignments between primer and sequence with enough overlap."""
    align_padding = len(primer) // 2

    for match_loc in _exact_kmer_matches(primer, bio_seq, min_overlap):
        start = match_loc - align_padding
        end = match_loc + len(primer) + align_padding

        # primer goes left: see _IUPAC_SIM and _iupac_match_score
        loprimer = primer.lower()
        losegment = bio_seq[start:end].lower()
        try:
            alns_t = pairwise2.align.localds(loprimer, losegment,
                                             _IUPAC_SIM,
                                             _GAP_OPEN_PENALTY,
                                             _GAP_EXT_PENALTY)
        except SystemError:
            localpos = _NONIUPAC_RE.search(losegment).start()
            globalp = start + localpos

            raise exceptions.BadCharacter(None, bio_seq, globalp, _IUPAC_ALPHA)

        alignify = functools.partial(_tup_to_aln, start)
        yield from map(alignify, alns_t)


def _better_alignment(what_end, aln_a, aln_b):
    """Out of the two alignments, choose the preferred one."""
    if aln_a < aln_b:
        return aln_b

    if aln_b < aln_a:
        return aln_a

    if what_end == End.HEAD:
        if aln_a.start < aln_b.start:
            return aln_b
        return aln_a

    assert what_end == End.TAIL

    if aln_a.start < aln_b.start:
        return aln_a
    return aln_b


def _best_primer_alignment(primers, bio_seq, min_overlap, what_end):
    """Select best alignment out of all alignments for all primers."""
    compare = functools.partial(_better_alignment, what_end)

    best_alignments = []
    for primer in primers:
        alignments = _min_overlap_alignments(primer, bio_seq, min_overlap)

        try:
            best = functools.reduce(compare, alignments)
        except TypeError:
            continue

        best_alignments.append(best)

    if not best_alignments:
        return None

    # Out of best alignments within primers, choose whatever primer
    # is most inward. We cannot simply use best global score here, as longer
    # primers will score better and there is no simple way to remedy that.
    if what_end == End.HEAD:
        return max(best_alignments, key=lambda aln: aln.start)

    assert what_end == End.TAIL
    return min(best_alignments, key=lambda aln: aln.start)


def cut_primers(record_cut_pairs, what_end, padding=0):
    """Yield records from which the primers are cut."""
    for record, cut_idx in record_cut_pairs:
        if what_end == End.HEAD:

            start_idx = max(cut_idx - padding, 0)
            cut_sequence = record.sequence[start_idx:]

        elif what_end == End.TAIL:
            end_idx = min(cut_idx + padding, len(record.sequence))
            cut_sequence = record.sequence[:end_idx]

        else:
            raise exceptions.Unreachable()

        yield fasta.FastaRecord(record.identifier, cut_sequence)


def primer_loci(primers, seq_records, min_overlap, what_end, revcomp=False):
    """Yield primer end points with each sequence record."""
    for record in seq_records:
        seq = record.sequence
        try:
            aln = _best_primer_alignment(primers, seq, min_overlap, what_end)
        except exceptions.BadCharacter as exc:
            # Add the identifier, which the callee did not know.
            exc.identifier = record.identifier
            raise exc

        if revcomp:
            rev_record = record.revcomp()
            rev_seq = rev_record.sequence
            rev_aln = _best_primer_alignment(primers, rev_seq,
                                             min_overlap, what_end)

            if rev_aln and \
                (not aln or \
                 rev_aln == _better_alignment(aln, rev_aln, what_end)):
                # continue as if our record was already reverse complemented
                record = rev_record
                seq = rev_seq
                aln = rev_aln

        _log.info("Identified best primer location for %r:\n%r",
                  record.identifier, aln)

        if not aln:
            continue

        if what_end == End.HEAD:
            cut_idx = aln.start + aln.span

        elif what_end == End.TAIL:
            cut_idx = aln.start

        else:
            raise exceptions.Unreachable()

        yield record, cut_idx
