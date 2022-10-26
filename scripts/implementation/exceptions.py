
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


"""This module is part of the HaploHIV pipeline and exports
exceptions. One of these exceptions is BadCharacter(...), which is
raised when a character is encountered in a context where it makes no
sense (or handling that character was not previously
implemented). Another is the the Unreachable() exception, for indicating
that there should not be any circumstance in which the codepath
leading to raising Unreachable is reached.

"""


import inspect


class BadCharacter(Exception):
    """For when we encounter a character that we cannot handle."""
    def __init__(self, identifier, sequence, position, alphabet):
        super().__init__(self)

        self.identifier = identifier
        self.sequence = sequence
        self.position = position
        self.alphabet = alphabet

    def __str__(self):
        padding = 15
        start = max(0, self.position - padding)
        end = min(len(self.sequence) - 1, self.position + padding)

        return "Encountered character outside of supported alphabet.\n" \
               f"Supported alphabet: {self.alphabet}\n" \
               f"Illegal character: {self.sequence[self.position]}\n" \
               f"Sequence identifier: {self.identifier}\n" \
               f"Sequence excerpt: {self.sequence[start:end]}\n"


class Unreachable(Exception):
    """For codepaths that logically cannot be reached."""
    def __init__(self):
        super().__init__(self)

        self.calling_frame = inspect.stack()[1]
        line = self.calling_frame.lineno
        name = self.calling_frame.frame.f_code.co_name
        filename = self.calling_frame.frame.f_code.co_filename

        msg = f"Line {line} of {filename}, in {name}, should be unreachable"
        msg += ", but was reached nonetheless."
        self.message = msg

    def __str__(self):
        return self.message
