#!/usr/bin/python

import sys
import re
import operator
import predicthaplo
import functools
import fasta

def enum_ph_fas_to_fasta(sample_name, enumerated_record):
    """Make header, sequence pair from a PH fas record with index."""
    index, ph_fas_record = enumerated_record
    identifier = f"{sample_name}_{index}_freq:{ph_fas_record.frequency}"

    return fasta.FastaRecord(identifier, ph_fas_record.sequence)

def main():
    input = sys.argv[1:len(sys.argv)-2]
    output = sys.argv[-2]

    fas_path = max(input, key=predicthaplo._reconstruction_window)

    with open(fas_path, "r") as f_in:
        ph_fas_records = predicthaplo.fas_records(f_in)
        sorted_records = sorted(ph_fas_records,
                                reverse=True,
                                key=operator.attrgetter("frequency"))

        sample_name = sys.argv[-1]
        transmogrify = functools.partial(enum_ph_fas_to_fasta, sample_name)

        fasta_records = map(transmogrify, enumerate(sorted_records))
        with open(output, "w") as f_out:
            fasta.write_records(f_out, fasta_records)

if __name__ == '__main__':
    main()