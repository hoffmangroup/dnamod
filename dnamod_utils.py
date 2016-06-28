# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

import csv
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))


def get_list(listname):
    returnlist = []
    with open(FILE_PATH + '/DNA_mod_site/static/whitelist/'+listname,
              'r') as targetlist:
        reader = csv.reader(targetlist, dialect='excel-tab')
        for line in reader:
            if line == []:
                continue
            if len(line) > 1:
                stringline = line[1]
                stringline = stringline.strip()
                returnlist.append(stringline)
    return returnlist
