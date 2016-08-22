# -*- coding: utf-8 -*-

"""Provides utility functions used in the construction of DNAmod.

-------------------------------------------------------------------------------
Copyright (C) 2016  Ankur Jai Sood and Coby Viner

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
-------------------------------------------------------------------------------
"""

from __future__ import with_statement, division, print_function

import csv
import os
import subprocess

_SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))

UNMOD_ALPH = ['A', 'C', 'G', 'T']


def get_constant(constant_name):
    CONSTANTS_SCRIPT = os.path.join(_SOURCE_DIR, 'constants.sh')

    process = subprocess.Popen([CONSTANTS_SCRIPT, constant_name],
                               stdout=subprocess.PIPE)
    return process.communicate()[0]


def _get_list_data(list_name):
    result = []

    with open(get_constant(list_name), 'rb') as targetlist:
        reader = csv.reader(targetlist, dialect='excel-tab')
        for line in reader:
            if line == []:
                continue
            if len(line) > 1:
                stringline = line[1]
                stringline = stringline.strip()
                result.append(stringline)
    return result


def get_whitelist():
    return _get_list_data('whitelist')


def get_blacklist():
    return _get_list_data('blacklist')
