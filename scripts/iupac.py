
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


"""This module is part of the HaploHIV pipeline and contains useful
facts about IUPAC nucleotide codes. This assumes DNA only.

"""


IUPAC_TO_NUCS = {
    "a": ["a"],
    "c": ["c"],
    "g": ["g"],
    "t": ["t"],
    "k": ["g", "t"],
    "m": ["a", "c"],
    "r": ["a", "g"],
    "s": ["g", "c"],
    "w": ["a", "t"],
    "y": ["c", "t"],
    "b": ["c", "g", "t"],
    "d": ["a", "g", "t"],
    "h": ["a", "c", "t"],
    "v": ["a", "c", "g"],
    "n": ["a", "c", "g", "t"],
    "x": ["a", "c", "g", "t"]
}


NUCS_TO_IUPAC = {
    **{tuple(sorted(nucs)): iupac for iupac, nucs in IUPAC_TO_NUCS.items()},
    **{("n",): "n", ("x",): "x"}
}


def _generate_iupac_complements():
    """Generate full IUPAC to complement table."""
    base_complements = {
        "a": "t",
        "t": "a",

        "c": "g",
        "g": "c"
    }
    result = {}

    for iupac, nucs in IUPAC_TO_NUCS.items():
        complement_nucs = tuple(sorted(base_complements[x] for x in nucs))
        result[iupac] = NUCS_TO_IUPAC[complement_nucs]

    return result


COMPLEMENT = _generate_iupac_complements()
