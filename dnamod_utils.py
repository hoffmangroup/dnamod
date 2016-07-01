# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

import csv
import os
import subprocess

SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))


def get_constant(constant_name):
    CONSTANTS_SCRIPT = os.path.join(SOURCE_DIR, 'constants.sh')

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
