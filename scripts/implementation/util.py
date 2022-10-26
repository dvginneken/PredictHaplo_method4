
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


"""This module contains utility functions for use throughout the
HaploHIV pipeline.

"""


import contextlib
import inspect
import itertools
import logging
import os


class Once:
    """Yield an object once, then always raise StopIteration."""
    def __init__(self, obj):
        self._obj = obj
        self._exhausted = False

    def __iter__(self):
        return self

    def __next__(self):
        if self._exhausted:
            raise StopIteration

        self._exhausted = True
        result = self._obj
        self._obj = None

        return result


def absdirname(path):
    """Absolute path to directory of passed path."""
    dirname = os.path.dirname(path)
    return os.path.abspath(dirname)


def adjacent(iterable, n=2):
    """Return iterator that zips adjacent elements in passed iterable."""
    iterator = iter(iterable)
    iterators = itertools.tee(iterable, n)

    # advance each iterator by the same amount as their index
    for index, iterator in enumerate(iterators):
        for _ in range(index):
            next(iterator)

    # zip the iterators, so that each is advanced at the same time
    # and so that their results are collected into a tuple
    return zip(*iterators)


def can_execute(path):
    """Whether we can execute the file located at path."""
    if not os.path.isfile(path):
        return False
    stat_result = os.stat(path)
    mode_bits = filemode_bitgroups(stat_result.st_mode)

    if user_owns(stat_result.st_uid) and mode_bits[0] & 1 != 0:
        return True

    if group_owns(stat_result.st_gid) and mode_bits[1] & 1 != 0:
        return True

    return mode_bits[2] & 1 != 0


@contextlib.contextmanager
def changed_working_dir(directory):
    """Use in a "with" statement to change working directory temporarily."""
    prev_cwd = os.getcwd()
    mkdir_p(directory)
    try:
        os.chdir(directory)
        yield os.getcwd()
    finally:
        os.chdir(prev_cwd)


def chunks(sequence, chunk_size):
    """Split sequence up into chunks of some size."""
    end = 0

    for start, end in adjacent(range(0, len(sequence), chunk_size)):
        yield sequence[start:end]

    final_chunk = sequence[end:len(sequence)]
    if final_chunk:
        yield final_chunk


def configure_logging(logfile_path, filemode="w", level=logging.INFO,
                      **kwargs):
    """Utility function for configuring the logging module."""
    if logfile_path is not None:
        logging.basicConfig(
            filename=logfile_path,
            filemode=filemode,
            level=level,
            **kwargs
        )
    else:
        logging.basicConfig(
            filemode=filemode,
            level=level,
            **kwargs
        )


def filemode_bitgroups(filemode):
    """Split a filemode integer into separate integers for each bit group."""
    all_bits = filemode & 7
    group_bits = (filemode >> 3) & 7
    user_bits = (filemode >> 6) & 7

    return user_bits, group_bits, all_bits


def find(predicate, iterable):
    """Make generator over index, element pairs for which predicate holds."""
    return filter(lambda tup: predicate(tup[1]), enumerate(iterable))


def find_first(predicate, iterable):
    """Find first index, element pair for which predicate holds."""
    generator = find(predicate, iterable)
    try:
        return next(generator)
    except StopIteration:
        return None


def group_owns(stat_gid):
    """Whether this process' group owns a file of which stat has given UID."""
    return any(x == stat_gid for x in [os.getgid(), os.getegid()])


def have_snakemake_object():
    """Does the caller have the snakemake object in scope?"""
    return is_name_in_scope("snakemake", nth_caller=1)


def is_name_in_scope(obj_name, nth_caller=0):
    """Is the passed name in the caller's scope?"""
    stack = inspect.stack()
    calling_frame = stack[1 + nth_caller].frame

    # This is the recommended way of ensuring that frame is deleted.
    try:
        result = obj_name in calling_frame.f_locals \
                 or obj_name in calling_frame.f_globals
    finally:
        del calling_frame

    return result


def iterator_exhausted(iterator):
    """Advance the passed iterator and report whether unsuccesful."""
    try:
        next(iterator)
    except StopIteration:
        return True
    else:
        return False


def mkdir_p(path):
    """Ensure that directory exists at path, creating parents as necessary."""
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise FileExistsError(path)
        # In this branch, if we did not error, then there is nothing left to do.

    else:
        # First recursively ensure that parent directory exists, then make path.
        parent_path = parent_directory(path)
        mkdir_p(parent_path)
        os.mkdir(path)


def parent_directory(path, n=1):
    """Get parent directory of path."""
    for _ in range(n):
        parent_relative = os.path.join(path, "..")
        path = os.path.normpath(parent_relative)

    return path


def user_owns(stat_uid):
    """Whether this process' user owns a file of which stat has given UID."""
    return any(x == stat_uid for x in [os.getuid(), os.geteuid()])


def unzip(nested_iterator):
    """Transpose outer two nested iterators."""
    return zip(*nested_iterator)


def write_formatted_template(out_file, template_iter, replace_by):
    """Iterate over the template, writing formatted elements to file."""
    for line in template_iter:
        formatted = line.format(**replace_by)
        out_file.write(formatted)
