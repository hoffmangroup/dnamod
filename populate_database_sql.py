#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

"""Creates the DNAmod database.

-------------------------------------------------------------------------------
Copyright (C) 2016  Ankur Jai Sood

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

from contextlib import contextmanager
import json
from more_itertools import peekable
from pysqlite2 import dbapi2 as sqlite3  # needed for latest SQLite
import sys
import unicodecsv as csv
import time
from datetime import datetime, timedelta
import socket
import argparse
import os

from Bio import Entrez
from suds.client import Client  # Using Suds web services client for soap

import dnamod_utils

Entrez.email = "DNAmod-L-request@listserv.utoronto.ca"
Entrez.tool = "DNAmod"

EXP_ALPH_TABLE_NAME = dnamod_utils.get_constant('exp_alph_table')
SEQ_TABLE_NAME = dnamod_utils.get_constant('seq_annot_table')
NATURE_TABLE_NAME = dnamod_utils.get_constant('nature_annot_table')

# Search Variables made up of CHEBI object attributes
SYNONYM_SEARCH_STR = 'Synonyms'
IUPAC_SEARCH_STR = 'IupacNames'
SMILES_SEARCH_STR = 'smiles'
FORMULA_SEARCH_STR = 'Formulae'
CHARGE_SEARCH_STR = 'charge'
MASS_SEARCH_STR = 'mass'
CITATION_SEARCH_STR = 'Citations'

ONTOLOGY_SEARCH_STR = "OntologyParents"
ONTOLOGY_HAS_ROLE = "has role"
ONTOLOGY_FP = "has functional parent"
ONTOLOGY_IS_A = "is a"
ONTOLOGY_IS_TAUTOMER = "is tautomer of"

RESET_TABLES = False
ChEBI_ID_PREFIX = "CHEBI:"

BLACK_LIST = []
WHITE_LIST = []

DNA_BASES = ['cytosine', 'thymine', 'adenine', 'guanine', 'uracil']

DATABASE_FILE_FULLPATH = dnamod_utils.get_constant('database')
ALPHABET_FILE_FULLPATH = dnamod_utils.get_constant('annot_exp_alph')
SEQ_REF_ANNOTS_FULLPATH = dnamod_utils.get_constant('annot_seq')
NATURE_REF_ANNOTS_FULLPATH = dnamod_utils.get_constant('annot_nature')
JSON_INDEX_FILE_FULLPATH = dnamod_utils.get_constant('json')
MANUALADD_FILE_FULLPATH = dnamod_utils.get_constant('manual_additions')

OTHER_BASE_ID = "CHEBI:other"
OTHER_BASE_NAME = "other"
OTHER_BASE_DEF = "A DNA base category added to encompass entities which do not fit into ATCG in the CHEBI ontology."

url = 'https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
client = Client(url)


class RequestMonitor:
    def __init__(self):
        self.requestCount = 0
        self.timeStart = time.time()

    def add_request(self):
        if self.requestCount < 3:
            self.requestCount += 1
        else:
            self.timeCurrent = time.time()
            self.waitTime = 1 - (self.timeCurrent - self.timeStart)
            if self.waitTime <= 1 and self.waitTime > 0:
                time.sleep(self.waitTime)
            self.timeStart = time.time()
            self.requestCount = 1


class CustomMod:
    def __init__(self, id, par, typ):
        self.chebiId = id
        self.parent = par
        self.type = typ


class CustomBase:
    def __init__(self, id, asciiname, definition):
        self.chebId = id
        self.chebiAsciiName = asciiname
        self.definition = definition


def check_time():
    utc = datetime.utcnow()
    time = utc.hour + utc.minute / 60. + utc.second / 3600.
    if time < 2 or time > 10:
        print("NCBI E-utilities restricts scripts to the off peak hours of "
              "9PM to 5AM EST. Please try running this script again during "
              "unrestricted hours. Current time (UTC): " + str(utc))
        sys.exit(0)


def field_not_found(field):
    print("WARNING: unable to find {} for a reference. "
          "Field left empty.".format(field), file=sys.stderr)


def search_exact_match(searchKey, requestMonitor, client):
    # Parameters for getLiteEntity() search
    searchCategory = 'ALL'
    maximumResults = sys.maxint
    starsCategory = 'THREE ONLY'

    # Save results from query into list
    requestMonitor.add_request()

    results = client.service.getLiteEntity(searchKey, searchCategory,
                                           maximumResults, starsCategory)
    if not results:
        result = 'Invalid Input'
        return result

    # Copy results.ListElement[] is messy to deal with
    # so I copy it over again into results
    results = results.ListElement

    # Initialize result as DNE in the case of no result
    result = 'DNE'

    # Iterate through all of the results and find the actual base
    for entities in results:
        if entities.chebiAsciiName == searchKey:
            result = entities

    return result


def search_for_bases(client, requestMonitor):
    # Initialize empy list
    result = []

    # Search CHEBI for bases and return lite entries
    for base in DNA_BASES:
        # print elements # Output for debugging
        result.append(search_exact_match(base, requestMonitor, client))

    return result


def filter_and_build_from_ontology(content, stars, client, requestMonitor):
    result = []

    for entity in content:
        temphold = entity
        entity = get_complete_entity(entity.chebiId, requestMonitor, client)

        if entity == 'Invalid Input' or entity == 'DNE':
            continue
        elif (entity.entityStar == stars and
              (temphold.type == ONTOLOGY_FP
               or temphold.type == ONTOLOGY_IS_A)):
            result.append(entity)

    return result


def find_manual_additions(base):
    additions = []
    with _read_csv_ignore_comments(MANUALADD_FILE_FULLPATH, True) as reader:
        for num, line in enumerate(reader):
            assert len(line) == 2
            if num == 0:  # header
                line.pop(0)
            else:
                id = "CHEBI:" + line[0]
                parent = line[1]
                if parent == base:
                    add = CustomMod(id, parent, ONTOLOGY_FP)
                    additions.append(add)
                    print("----- Manual Addition: {}".format(add.chebiId))
    return additions


def get_children(bases, requestMonitor, client):
    modDictionary = {}
    for base in bases:
        requestMonitor.add_request()

        result = client.service.getOntologyChildren(base.chebiId)
        print("----- BASE: {}".format(base.chebiAsciiName))
        result = result.ListElement
        manualadds = find_manual_additions(base.chebiAsciiName)
        for add in manualadds:
            result.append(add)
        result = filter_and_build_from_ontology(result, 3, client,
                                                requestMonitor)
        additionalChildren = get_further_children(result, client,
                                                  requestMonitor)

        for child in additionalChildren:
            result.append(child)

        result = set(result)
        result = list(result)

        modDictionary[base.chebiAsciiName] = result
        
    requestMonitor.add_request()
    print("----- BASE: {}".format("other"))
    manualadds = find_manual_additions("other")
    result = []
    for add in manualadds:
        result.append(add)
    result = filter_and_build_from_ontology(result, 3, client,
                                            requestMonitor)
    additionalChildren = get_further_children(result, client,
                                              requestMonitor)

    for child in additionalChildren:
        result.append(child)

    result = set(result)
    result = list(result)

    modDictionary["other"] = result

    return modDictionary


def get_recursive_children(entity, client, childrenverified, requestMonitor,
                           additionalChildren):
    if entity.chebiAsciiName not in BLACK_LIST:
        if childrenverified or [child for child in additionalChildren
                                if child['chebiId'] == entity.chebiId and
                                child['verifiedstatus']]:
            # TODO improve the way duplicate verified ontology entries are
            # handled
            entity['verifiedstatus'] = 1
        else:
            entity['verifiedstatus'] = 0

        print("---------- CHILD of BASE: {0:100} Verified: {1} "
              "".format(entity.chebiAsciiName, childrenverified))

        requestMonitor.add_request()

        result = client.service.getOntologyChildren(entity.chebiId)

        if result:
            result = result.ListElement
            result = filter_and_build_from_ontology(result, 3, client,
                                                    requestMonitor)

            for child in result:
                recursivestep = get_recursive_children(child, client,
                                                       childrenverified,
                                                       requestMonitor,
                                                       additionalChildren)
                result = result + recursivestep

            result = set(result)
            result = list(result)

            additionalChildren.extend(result)
    else:
        entity['verifiedstatus'] = 0

    return additionalChildren


def get_further_children(entities, client, requestMonitor):
    childrenverified = False
    additionalChildren = []

    for entity in entities:
        if entity.chebiAsciiName in WHITE_LIST:
            entity['verifiedstatus'] = 1
            childrenverified = True
        else:
            entity['verifiedstatus'] = 0
            childrenverified = False

        additionalChildren = get_recursive_children(entity, client,
                                                    childrenverified,
                                                    requestMonitor,
                                                    additionalChildren)

    return additionalChildren


def get_complete_entity(CHEBIid, requestMonitor, client):
    requestMonitor.add_request()

    result = client.service.getCompleteEntity(CHEBIid)
    return result


def get_complete_bases(bases, requestMonitor, client):
    result = []
    for base in bases:
        result.append(get_complete_entity(base.chebiId, requestMonitor,
                      client))
    otherBase = CustomBase(OTHER_BASE_ID, OTHER_BASE_NAME, OTHER_BASE_DEF)
    result.append(otherBase)
    return result


def create_base_table(conn, sql_conn_cursor, bases):
    if RESET_TABLES is True:
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS base''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS base
                            (baseid TEXT PRIMARY KEY NOT NULL,
                             commonname text,
                             basedefinition text)''')
    conn.commit()

    for base in bases:
        sql_conn_cursor.execute('''UPDATE base SET baseid=?, commonname=?,
                                basedefinition=? WHERE commonname=?''',
                                (base.chebiAsciiName[0].capitalize(),
                                 base.chebiAsciiName,
                                 base.chebiAsciiName, base.chebiAsciiName[0]))
        conn.commit()
        sql_conn_cursor.execute("INSERT OR IGNORE INTO base VALUES(?,?,?)",
                                (base.chebiAsciiName[0].capitalize(),
                                 base.chebiAsciiName, base.definition))
    conn.commit()


def concatenate_list(child, attribute):
    if attribute in dir(child):
        return [synonym_data.data for synonym_data in
                getattr(child, attribute)]
    else:
        return []


def get_entity(child, attribute):
    if attribute in dir(child):
        return ([getattr(child, attribute)])
    else:
        return []


def get_ontology_data(child, attribute, selectors):
    return [ontologyitem for ontologyitem in
            getattr(child, attribute) if ontologyitem.type in selectors]


def get_full_citation(PMID):
    if PMID[0] == 'I':
        return

    print("Adding Citation: {}".format(PMID))
    result = []
    # isbook = False # Unused at the moment
    isarticle = False

    handle = Entrez.efetch("pubmed", id=PMID, retmode="xml")
    records = Entrez.read(handle)

    articleTitle = []
    publicationDate = []
    authors = []

    journalName = None
    Volume = None
    Issue = None
    publisherName = None
    publisherLocation = None

    for record in records['PubmedArticle']:
        if 'MedlineCitation' in record:
            isarticle = True
            article = record['MedlineCitation']['Article']
            journalrecord = record['MedlineCitation']['Article']['Journal']
            daterecord = record['PubmedData']['History']
            if'ArticleTitle' in article.keys():
                articleTitle = article['ArticleTitle']
            else:
                articleTitle = None
            if'AuthorList' in article.keys():
                authors = article['AuthorList']
            else:
                authors = None
            '''if'PubDate' in article.keys():
                publicationDate = article['PubDate']
            elif'ArticleDate' in article.keys():
                publicationDate = article['ArticleDate']'''

            publicationDate = [date for date in daterecord if
                               date.attributes['PubStatus'] == "pubmed"]

            if 'Title' in journalrecord.keys():
                journalName = journalrecord['Title']
            else:
                journalName = None

            if 'JournalIssue' in journalrecord.keys():
                journalIssue = journalrecord['JournalIssue']
                if 'Volume' in journalIssue.keys():
                    Volume = journalIssue['Volume']
                else:
                    Volume = None
                if 'Issue' in journalIssue.keys():
                    Issue = journalIssue['Issue']
                else:
                    Issue = None
            else:
                Volume = None
                Issue = None

        else:
            # XXX TODO refactor
            # isbook = True # Unused at the moment
            article = record['BookDocument']['Book']

            if'BookTitle' in article.keys():
                articleTitle = article['BookTitle']
            else:
                articleTitle = None
            if'AuthorList' in article.keys():
                authors = article['AuthorList']
            else:
                authors = None
            if'PubDate' in article.keys():
                publicationDate = article['PubDate']
            if 'Publisher' in article.keys():
                publisherName = article['Publisher']['PublisherName']
                publisherLocation = article['Publisher']['PublisherLocation']
            else:
                publisherName = None

    handle.close()

    # XXX TODO refactor not found instances

    if articleTitle:
        result.append(articleTitle.encode('utf-8'))
    else:
        result.append('')
        field_not_found('title')

    if isarticle:
        if publicationDate:
            publicationDate = publicationDate[0]

            date = publicationDate['Year']
            result.append(date)
        else:
            result.append('')
            field_not_found('date')
    else:
        if publicationDate:
            date = publicationDate['Year']
            result.append(date)
        else:
            result.append('')
            field_not_found('date (non-article entry)')
    if authors:
        result.append("{0}, {1}, et al.".format(
                      authors[0]['LastName'].encode("utf-8"),
                      authors[0]['Initials'].encode("utf-8")))
    else:
        result.append('')
        field_not_found('author(s)')

    if journalName:
        result.append(journalName.encode('utf-8'))
    else:
        field_not_found('journal name')
        result.append('')

    if Volume:
        result.append(Volume)
    else:
        field_not_found('volume')
        result.append('')

    if Issue:
        result.append(Issue)
    else:
        field_not_found('issue')
        result.append('')

    if publisherName:
        result.append(publisherName.encode('utf-8'))
    elif journalName is None:
        field_not_found('publisher name')
        result.append('')

    if publisherLocation:
        result.append(publisherLocation)
    elif Volume is None:
        field_not_found('publisher location')
        result.append('')

    return result


def create_other_tables(conn, sql_conn_cursor, children, bases):
    # Reset Tables
    if RESET_TABLES is True:
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS covmod''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS names''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS baseprops''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS citation_lookup''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS roles_lookup''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS citations''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS roles''')
        sql_conn_cursor.execute('''DROP TABLE IF EXISTS modbase''')
        conn.commit()

    # Create Tables
    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS modbase
                            (nameid text PRIMARY KEY NOT NULL,
                             position text,
                             baseid text,
                             formulaid text,
                             cmodid integer,
                             verifiedstatus integer,
                             FOREIGN KEY(baseid) REFERENCES base(baseid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             FOREIGN KEY(nameid) REFERENCES names(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             FOREIGN KEY(formulaid) REFERENCES
                             baseprops(formulaid) ON DELETE CASCADE
                             ON UPDATE CASCADE,
                             FOREIGN KEY(cmodid) REFERENCES
                             covmod(cmodid) ON DELETE CASCADE
                             ON UPDATE CASCADE)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS covmod
                            (cmodid integer PRIMARY KEY NOT NULL,
                             symbol text,
                             netcharge text,
                             definition text)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS names
                            (chebiname text,
                             nameid text PRIMARY KEY NOT NULL,
                             iupacname text,
                             othernames text,
                             smiles text)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS baseprops
                            (formulaid text PRIMARY KEY NOT NULL,
                             avgmass text)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS citations
                            (citationid text PRIMARY KEY NOT NULL,
                             title text,
                             pubdate text,
                             authors text,
                             journalnameorpublishername text,
                             volumeorpublisherlocation text,
                             issue text)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS roles
                            (roleid text PRIMARY KEY NOT NULL,
                             role text)''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS citation_lookup
                            (nameid text,
                             citationid text,
                             FOREIGN KEY(nameid) REFERENCES modbase(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             FOREIGN KEY(citationid) REFERENCES
                             citations(citationid) ON DELETE CASCADE
                             ON UPDATE CASCADE,
                             PRIMARY KEY(nameid, citationid))''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS roles_lookup
                            (nameid text,
                             roleid text,
                             FOREIGN KEY(nameid) REFERENCES modbase(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             FOREIGN KEY(roleid) REFERENCES roles(roleid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             PRIMARY KEY(nameid, roleid))''')

    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS modbase_parents
                            (nameid text,
                             parentid text,
                             UNIQUE (nameid, parentid) ON CONFLICT IGNORE,
                             FOREIGN KEY(nameid) REFERENCES modbase(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE,
                             FOREIGN KEY(parentid) REFERENCES modbase(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE)''')

    conn.commit()

    # Populate Tables
    added_entry_IDs = []
    formulafill = 0

    for base in bases:
        # Retreive and process list of childs of current base
        for child in children[base.chebiAsciiName]:
            added_entry_IDs.append(child.chebiId)

            # Parse CHEBI datastructure for relevant info
            synonyms = concatenate_list(child, SYNONYM_SEARCH_STR)
            newnames = []
            for name in synonyms:
                name = str(name.encode('utf-8')).lower()
                if name != child.chebiAsciiName:
                    newnames.append(name)
            synonyms = list(set(newnames))

            iupac = concatenate_list(child, IUPAC_SEARCH_STR)
            smiles = get_entity(child, SMILES_SEARCH_STR)

            formula = concatenate_list(child, FORMULA_SEARCH_STR)
            if not formula:
                formulafill = formulafill + 1
                formula = formulafill

            charge = get_entity(child, CHARGE_SEARCH_STR)
            mass = get_entity(child, MASS_SEARCH_STR)
            citations = concatenate_list(child, CITATION_SEARCH_STR)

            roles = get_ontology_data(child, ONTOLOGY_SEARCH_STR,
                                      [ONTOLOGY_HAS_ROLE])
            role_names = [role.chebiName for role in roles]
            role_ids = [role.chebiId for role in roles]

            # Populate roles and citations tables with unique data
            for role in range(len(roles)):
                sql_conn_cursor.execute('''SELECT count(*) FROM roles
                                        WHERE roleid = ?''', (role_ids[role],))
                data = sql_conn_cursor.fetchone()[0]
                if data == 0:
                    sql_conn_cursor.execute('''UPDATE roles SET role=?
                                            WHERE roleid=?''',
                                            (role_names[role], role_ids[role]))
                    conn.commit()
                    sql_conn_cursor.execute('''INSERT OR IGNORE
                                            INTO roles VALUES(?,?)''',
                                            (role_ids[role], role_names[role]))
                    conn.commit()

            for citation in citations:
                sql_conn_cursor.execute('''SELECT * FROM citations
                                        WHERE citationid = ?''',
                                        [citation])
                data = sql_conn_cursor.fetchone()
                if True:
                    citationinfo = get_full_citation(citation)
                    if citationinfo:
                        citationinfo_uni = [info.decode('utf-8') for info in
                                            citationinfo]
                    else:
                        citationinfo_uni = [" ", " ", " ", " ", " ", " "]

                    sql_conn_cursor.execute('''UPDATE citations SET title=?,
                                            pubdate=?, authors=?,
                                            journalnameorpublishername=?,
                                            volumeorpublisherlocation=?,
                                            issue=? WHERE citationid=?''',
                                            (citationinfo_uni[0],
                                             citationinfo_uni[1],
                                             citationinfo_uni[2],
                                             citationinfo_uni[3],
                                             citationinfo_uni[4],
                                             citationinfo_uni[5],
                                             citation))
                    conn.commit()

                    sql_conn_cursor.execute('''INSERT OR IGNORE
                                            INTO citations
                                            VALUES(?,?,?,?,?,?,?)''',
                                            (citation, citationinfo_uni[0],
                                             citationinfo_uni[1],
                                             citationinfo_uni[2],
                                             citationinfo_uni[3],
                                             citationinfo_uni[4],
                                             citationinfo_uni[5]))
                    conn.commit()

            sql_conn_cursor.execute('''INSERT OR IGNORE
                                    INTO baseprops VALUES(?,?)''',
                                    (str(formula), str(mass)))
            sql_conn_cursor.execute('''INSERT OR IGNORE
                                    INTO names VALUES(?,?,?,?,?)''',
                                    (child.chebiAsciiName, child.chebiId,
                                     str(iupac), str(synonyms), str(smiles)))
            sql_conn_cursor.execute('''INSERT OR IGNORE
                                    INTO covmod VALUES(NULL,?,?,?)''',
                                    ('0', str(charge), child.definition))
            conn.commit()
            rowid = sql_conn_cursor.lastrowid

            sql_conn_cursor.execute('''INSERT OR IGNORE
                                    INTO modbase VALUES(?,?,?,?,?,?)''',
                                    (child.chebiId, '0',
                                     base.chebiAsciiName[0].capitalize(),
                                     str(formula),
                                     rowid, child['verifiedstatus']))
            conn.commit()

            for role in range(len(roles)):
                sql_conn_cursor.execute('''INSERT OR IGNORE
                                        INTO roles_lookup VALUES(?,?)''',
                                        (child.chebiId, role_ids[role]))
            conn.commit()

            for citation in citations:
                sql_conn_cursor.execute('''INSERT OR IGNORE
                                        INTO citation_lookup VALUES(?,?)''',
                                        (child.chebiId, citation))
            conn.commit()

        # for each base, go through children again and annotate parents within
        # database
        # NB: only bases with other modified bases are annotated as parents
        #     the non-modified (A/C/G/T/U) base is annotated elsewhere
        for child in children[base.chebiAsciiName]:
            parent_entires = get_ontology_data(child, ONTOLOGY_SEARCH_STR,
                                               [ONTOLOGY_FP, ONTOLOGY_IS_A])

            parents_to_annot = [parent.chebiId for parent in parent_entires
                                if parent.chebiId in added_entry_IDs]

            for parent_to_annot in parents_to_annot:
                sql_conn_cursor.execute('''INSERT OR IGNORE
                                        INTO modbase_parents VALUES(?,?)''',
                                        (child.chebiId, parent_to_annot))
            conn.commit()


@contextmanager
def _read_csv_ignore_comments(file_name, enforce_header_num_cols=None):
    try:
        with open(file_name, 'rb') as file:
            # use a generator expression to ignore comments
            reader = peekable(csv.reader((row for row in file if not
                              row.startswith('#')), delimiter="\t"))
            if enforce_header_num_cols:
                # look at header, leaving it in gen.
                num_cols = len(reader.peek())
                yield (row + ([None] * (num_cols - len(row))) for row in
                       reader)
            else:
                yield reader
    finally:
        pass


def _create_modbase_annot_table(conn, sql_conn_cursor, header, table_name):
    col_names_create_spec = ''
    for num, col in enumerate(header):
        col_names_create_spec += '[{}] text'.format(col)
        if (num < len(header)):
            col_names_create_spec += ','

    # TODO fix below to perform stringent validation
    # (see: http://stackoverflow.com/questions/25387537/
    # sqlite3-operationalerror-near-syntax-error)
    # maybe ignore, since never "user" sourced (since this is static)
    sql_conn_cursor.execute('''CREATE TABLE IF NOT EXISTS {}
                            (nameid text,
                             {}
                             FOREIGN KEY(nameid) REFERENCES modbase(nameid)
                             ON DELETE CASCADE ON UPDATE CASCADE)
                             '''.format(table_name, col_names_create_spec))
    conn.commit()


# Returns the id actually used, since it may be modified (i.e. prefix added)
# The id returned may be None, if it did not exist in the database
def _add_annots_for_id(id, sql_conn_cursor, line, last_col_num, table_name):
    # add ID prefix, if missing
    if not id.startswith(ChEBI_ID_PREFIX):
        id = ChEBI_ID_PREFIX + id

    if (sql_conn_cursor.execute("SELECT * FROM modbase WHERE nameid=?",
                                (id,)).fetchone()):
        # allow for one column per column in the header, plus the ID column
        col_names_wildcards = '?,' * (1 + len(line[:last_col_num]))
        col_names_wildcards = col_names_wildcards[:-1]  # remove final comma

        sql_conn_cursor.execute("INSERT OR IGNORE INTO {} VALUES({})"
                                "".format(table_name, col_names_wildcards),
                                tuple([id] + line[:last_col_num]))

        return id  # id exists and is now valid
    else:
        print("WARNING: an annotation for {} was skipped, "
              "since it is not in the database.".format(id), file=sys.stderr)
        return None  # indicate that id does not exist


def create_exp_alph_table(conn, sql_conn_cursor, exp_alph_file_name):
    print("---------- Adding Nomenclature ----------")

    with _read_csv_ignore_comments(exp_alph_file_name, True) as reader:
        for num, line in enumerate(reader):
            assert len(line) > 2  # min. of three columns

            if num == 0:  # header
                # remove the first column, since it is our foreign key and not
                # displayed
                line.pop(0)
                _create_modbase_annot_table(conn, sql_conn_cursor, line,
                                            EXP_ALPH_TABLE_NAME)
            else:
                ids = line[0].split(",")
                line.pop(0)  # remove the ID, since it is processed above
                for id in ids:
                    id = _add_annots_for_id(id, sql_conn_cursor,
                                            line, len(line),
                                            EXP_ALPH_TABLE_NAME)
    conn.commit()


def create_annot_citation_tables(conn, sql_conn_cursor, ref_annots_file_name,
                                 table_name):
    print("---------- Adding Custom Annotations ----------")

    with _read_csv_ignore_comments(ref_annots_file_name) as reader:
        for num, line in enumerate(reader):
            if num == 0:  # header
                # remove the first column, since it is our foreign key and
                # not displayed
                line.pop(0)
                _create_modbase_annot_table(conn, sql_conn_cursor, line,
                                            table_name)
            else:
                # TODO refactor to check or at least output a descriptive error
                # if PMID is not last col
                references = line[-1].split(",")

                for reference in references:
                    if not reference:
                        continue

                    citationinfo = get_full_citation(reference)
                    citationinfo_uni = [info.decode('utf-8') for info in
                                        citationinfo]

                    sql_conn_cursor.execute('''UPDATE citations
                                            SET title=?, pubdate=?, authors=?
                                            WHERE citationid=?''',
                                            (citationinfo_uni[0],
                                             citationinfo_uni[1],
                                             citationinfo_uni[2], reference))
                    conn.commit()

                    sql_conn_cursor.execute('''INSERT OR IGNORE INTO citations
                                            VALUES(?,?,?,?,?,?,?)''',
                                            (reference,
                                             citationinfo_uni[0],
                                             citationinfo_uni[1],
                                             citationinfo_uni[2],
                                             citationinfo_uni[3],
                                             citationinfo_uni[4],
                                             citationinfo_uni[5]))

                ids = line[0].split(",")
                line.pop(0)  # remove the ID, since it is processed above
                for id in ids:
                    id = _add_annots_for_id(id, sql_conn_cursor, line,
                                            len(line), table_name)

                    # id may now be None, if it is not (yet) in the database
                    if id:
                        for reference in references:
                            # some annotations may lack a reference
                            if not reference:
                                continue

                            sql_conn_cursor.execute('''INSERT OR IGNORE
                                                    INTO citation_lookup
                                                    VALUES(?,?)''',
                                                    (id, reference))

        conn.commit()


def create_search_index(conn, sql_conn_cursor, JSON_fullpath):
    # create or re-create the JSON
    with open(JSON_fullpath, 'w+') as JSON:
        feeds = []
        sql_conn_cursor.execute("SELECT nameid FROM modbase")
        result = sql_conn_cursor.fetchall()

        for modification in result:
            sql_conn_cursor.execute('''SELECT Name FROM expanded_alphabet
                                    WHERE nameid =?''', (modification))
            nomenclature = sql_conn_cursor.fetchone()
            name = nomenclature
            sql_conn_cursor.execute('''SELECT * FROM names WHERE nameid = ?''',
                                    (modification))
            data = sql_conn_cursor.fetchone()
            chebiname = data[0]
            chebiid = data[1]
            iupacname = data[2]
            synonyms = data[3]

            sql_conn_cursor.execute('''SELECT formulaid, verifiedstatus
                                    FROM modbase WHERE nameid = ?''',
                                    (modification))
            furtherdata = sql_conn_cursor.fetchone()
            formula = furtherdata[0]
            verified = furtherdata[1]

            sql_conn_cursor.execute('''SELECT Abbreviation
                                    FROM expanded_alphabet WHERE nameid = ?''',
                                    (modification))
            abbrevdata = sql_conn_cursor.fetchone()
            abbreviation = abbrevdata

            sql_conn_cursor.execute('''SELECT Symbol FROM expanded_alphabet
                                    WHERE nameid = ?''', (modification))
            symbol = sql_conn_cursor.fetchone()

            common_name = name if name else chebiname

            writedata = {
                'CommonName': common_name,
                'ChEBIId': chebiid,
                'IUPACName': iupacname,
                'Synonyms': synonyms,
                'ChemicalFormula': formula,
                'Abbreviation': abbreviation,
                'Verified': verified,
                'Symbol': symbol,
                'Refname': chebiname
            }

            if (writedata['ChEBIId'] is not None):
                feeds.append(writedata.copy())

        # write the JSON
        json.dump(feeds, JSON, sort_keys=True, indent=4)

    with open(JSON_fullpath, 'r') as JSON:
        # check the JSON
        try:
            json.load(JSON)
        except ValueError as JSON_err:
            print("WARNING: invalid JSON ({}). Search is unlikely"
                  "to function.".format(JSON_err),
                  file=sys.stderr)


def populate_tables(conn, sql_conn_cursor, bases, children, client):
    create_base_table(conn, sql_conn_cursor, bases)
    create_other_tables(conn, sql_conn_cursor, children, bases)
    create_exp_alph_table(conn, sql_conn_cursor, ALPHABET_FILE_FULLPATH)
    create_annot_citation_tables(conn, sql_conn_cursor,
                                 SEQ_REF_ANNOTS_FULLPATH, SEQ_TABLE_NAME)
    create_annot_citation_tables(conn, sql_conn_cursor,
                                 NATURE_REF_ANNOTS_FULLPATH, NATURE_TABLE_NAME)
    print("4/5 Creating Search Index...")
    create_search_index(conn, sql_conn_cursor, JSON_INDEX_FILE_FULLPATH)


def fix_verified_status(conn, sql_conn_cursor, client, requestMonitor):
    sql_conn_cursor.execute('''SELECT nameid FROM modbase''')
    result = sql_conn_cursor.fetchall()
    ids2unverify = []
    for nameid in result:
        sql_conn_cursor.execute('''SELECT verifiedstatus FROM modbase
                                WHERE nameid = ?''', nameid)
        verified = sql_conn_cursor.fetchone()
        if verified[0] == 1:
            entity = get_complete_entity(nameid, requestMonitor, client)
            for ontologyItem in entity.OntologyParents:
                if ontologyItem.type == ONTOLOGY_IS_TAUTOMER:
                    tautId = ontologyItem.chebiId
                    ids2unverify.append(tautId)

        sql_conn_cursor.execute('''SELECT othernames FROM names
                                WHERE nameid = ?''', nameid)
        synonyms = sql_conn_cursor.fetchone()
        if verified[0] == 1 and synonyms[0]:
            synonyms = synonyms[0]
            synonyms = synonyms[1:-1]
            synonyms = synonyms.split(', ')
            for name in synonyms:
                name = name[1:-1]
                sql_conn_cursor.execute('''SELECT chebiname FROM names
                                        WHERE chebiname = ?''', (name,))
                matchedName = sql_conn_cursor.fetchone()
                if matchedName:
                    ids2unverify.append(nameid[0])

    uniqueAbbreviations = []
    with _read_csv_ignore_comments(ALPHABET_FILE_FULLPATH, True) as reader:
            for num, line in enumerate(reader):
                if num == 0:
                    line.pop(0)
                else:
                    id = line[0].split(",")
                    abbreviation = line[1]
                    if abbreviation in uniqueAbbreviations:
                        if id[0] not in ids2unverify:
                            ids2unverify.append(id[0])
                    else:
                        uniqueAbbreviations.append(abbreviation)

    print('Verified status revoked for: ', ids2unverify)
    for id in ids2unverify:
        sql_conn_cursor.execute('''UPDATE modbase SET verifiedstatus = 0
                                WHERE nameid = ?''', (id,))
    conn.commit()

WHITE_LIST = dnamod_utils.get_whitelist()
BLACK_LIST = dnamod_utils.get_blacklist()

requestMonitor = RequestMonitor()
socket.setdefaulttimeout(300)
#check_time()

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reset', action='store_true')
args = parser.parse_args()

if args.reset == True:
    RESET_TABLES = True

if RESET_TABLES == True:
    if os.path.exists(DATABASE_FILE_FULLPATH):
        os.remove(DATABASE_FILE_FULLPATH)
    print("Reseting Database...")

conn = sqlite3.connect(DATABASE_FILE_FULLPATH)

sql_conn_cursor = conn.cursor()

if RESET_TABLES == False:
    sql_conn_cursor.execute('''PRAGMA foreign_keys = ON''')

conn.commit()
os.chmod(DATABASE_FILE_FULLPATH, 0755)

print("1/5 Searching for bases...")
bases = search_for_bases(client, requestMonitor)

print("2/5 Searching for children...")
children = get_children(bases, requestMonitor, client)
bases = get_complete_bases(bases, requestMonitor, client)

print("3/5 Creating tables...")
populate_tables(conn, sql_conn_cursor, bases, children, client)

print("5/5 Finishing up...")
fix_verified_status(conn, sql_conn_cursor, client, requestMonitor)

conn.close()

print("Done!")
