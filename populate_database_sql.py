#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

'''
Ankur Jai Sood
29/5/2015

Populate_database.py
Function:
1. Performs search of CHEBI database for DNA bases
2. Returns chebiId and chebiAsciiName of bases
3. Searches CHEBI database for all entities of which
   the DNA bases are functional parents
4. Using above results, populates DNA post-transciptional modification table
'''

# Using Suds web services client for soap
import csv
import os
import sqlite3
import sys
from suds.client import Client
from Bio import Entrez
from sys import maxint

import dnamod_utils

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
RESET_TABLES = True

# Program Constants
BLACK_LIST = []
WHITE_LIST = []

DNA_BASES = ['cytosine', 'thymine', 'adenine', 'guanine', 'uracil']

FILE_PATH = os.path.dirname(os.path.abspath(__file__))

print("Operating from: {}".format(FILE_PATH))

DATABASE_FILE_FULLPATH = os.path.join(FILE_PATH, "DNA_mod_database.db")

REF_ANNOTS_FULLPATH = os.path.join(FILE_PATH, "ref_annots_sequencing.txt")

url = 'https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
client = Client(url)


def search_exact_match(searchKey, client):

    # Parameters for getLiteEntity() search
    searchCategory = 'ALL'
    maximumResults = maxint
    starsCategory = 'THREE ONLY'

    # Save results from query into list
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


def search_for_bases(client):
    # Initialize empy list
    result = []

    # Search CHEBI for bases and return lite entries
    for base in DNA_BASES:
        # print elements # Output for debugging
        result.append(search_exact_match(base, client))

    return result


def filter_stars(content, stars, client):
    result = []
    for entity in content:
        # print content[elements] # For Debugging
        # content[elements] = search_exact_match(content[elements].chebiName,
        #                                        client)
        temphold = entity
        # print temphold # For debugging
        #time.sleep(1)
        entity = get_complete_entity(entity.chebiId, client)
        # print content[elements] # For Debugging
        if entity == 'Invalid Input' or entity == 'DNE':
            continue
        elif (entity.entityStar == stars and
              (temphold.type == 'has functional parent' 
               or temphold.type == 'is a')):
            result.append(entity)
    return result


def get_children(bases, client):
    modDictionary = {}
    for base in bases:
        result = client.service.getOntologyChildren(base.chebiId)
        print("----- BASE: {}".format(base.chebiAsciiName))
        result = result.ListElement
        result = filter_stars(result, 3, client)
        additionalChildren = get_further_children(result, client)

        for child in additionalChildren:
            result.append(child)

        result = set(result)
        result = list(result)

        modDictionary[base.chebiAsciiName] = result

    return modDictionary


def get_recursive_children(entity, client, childrenverified, additionalChildren):
    if entity.chebiAsciiName not in BLACK_LIST:
        if childrenverified:
            entity['verifiedstatus'] = 1;
        else:
            entity['verifiedstatus'] = 0;
        print("---------- CHILD of BASE: {0:100} Verified: {1} ".format(entity.chebiAsciiName, childrenverified))
        result = client.service.getOntologyChildren(entity.chebiId)
        if result:
            result = result.ListElement
            result = filter_stars(result, 3, client)
            for child in result:
                recursivestep = get_recursive_children(child, client, childrenverified, additionalChildren)
                result = result + recursivestep
            result = set(result)
            result = list(result)
            additionalChildren.extend(result)
    return additionalChildren

def get_further_children(entities, client):
    childrenverified = False;
    additionalChildren = []
    for entity in entities:
        if entity.chebiAsciiName in WHITE_LIST:
            entity['verifiedstatus'] = 1;
            childrenverified = True;
        else:
            entity['verifiedstatus'] = 0;
            childrenverified = False;
        additionalChildren = get_recursive_children(entity, client, childrenverified, additionalChildren)
    return additionalChildren


def get_complete_entity(CHEBIid, client):
    result = client.service.getCompleteEntity(CHEBIid)
    return result


def get_complete_bases(bases, client):
    result = []
    for base in bases:
        result.append(get_complete_entity(base.chebiId, client))
    return result


def create_base_table(conn, bases):
    c = conn.cursor()

    c.execute('''PRAGMA foreign_keys = ON''')
    conn.commit()

    if RESET_TABLES is True:
        c.execute('''DROP TABLE IF EXISTS base''')

    c.execute('''CREATE TABLE IF NOT EXISTS base
                (baseid TEXT PRIMARY KEY NOT NULL,
                 commonname text,
                 basedefinition text)''')
    conn.commit()

    for base in bases:
        c.execute("UPDATE base SET commonname=?, basedefinition=? WHERE commonname=?",(base.chebiAsciiName, base.chebiAsciiName, base.chebiAsciiName[0]))
        conn.commit()
        c.execute("INSERT OR IGNORE INTO base VALUES(?,?,?)",
                  (base.chebiAsciiName[0],
                   base.chebiAsciiName, base.definition))
    conn.commit()
    conn.close()


def concatenate_list(child, attribute):
    if attribute in dir(child):
        return [synonym_data.data for synonym_data in getattr(child, attribute)]
    else:
        return []


def get_entity(child, attribute):
    if attribute in dir(child):
        return ([getattr(child, attribute)])
    else:
        return []

def get_roles(child, attribute, selector):
    return [ontologyitem for ontologyitem in
            getattr(child, attribute) if ontologyitem.type == selector] 


def get_full_citation(PMID):
    print("Adding Citation: ",PMID)
    result = []
    # isbook = False # Unsed at the moment
    isarticle = False
    Entrez.email = "jai.sood@hotmail.com"
    handle = Entrez.efetch("pubmed", id=PMID, retmode="xml")
    records = Entrez.parse(handle)

    articleTitle = []
    publicationDate = []
    authors = []

    for record in records:
        if 'MedlineCitation' in record.keys():
            isarticle = True
            article = record['MedlineCitation']['Article']
            if'ArticleTitle' in article.keys():
                articleTitle = article['ArticleTitle']
            else:
                articleTitle = None
            if'AuthorList' in article.keys():
                authors = article['AuthorList']
            else:
                authors = None
            if'ArticleDate' in article.keys():
                publicationDate = article['ArticleDate']
            else:
                publicationDate = None
        else:
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
            else:
                publicationDate = None

    handle.close()

    if articleTitle:
        result.append(articleTitle.encode('utf-8'))
    else:
        result.append('Title Not Found')

    if isarticle:
        if publicationDate:
            date = (publicationDate[0]['Month'] + '-' +
                    publicationDate[0]['Day'] + '-' +
                    publicationDate[0]['Year'])
            result.append(date)
        else:
            result.append('Publication Date Not Found')
    else:
        if publicationDate:
            date = (publicationDate['Month'] + '-' + publicationDate['Day']
                    + '-' + publicationDate['Year'])
            result.append(date)
        else:
            result.append('Publication Date Not Found')
    if authors:
        result.append("{0}, {1}, et al.".format(authors[0]['LastName'].encode("utf-8"),
                                              authors[0]['Initials'].encode("utf-8")))
    else:
        result.append('Authors Not Found')

    return result


def create_other_tables(conn, children, bases):
    c = conn.cursor()

    c.execute('''PRAGMA foreign_keys = ON''')
    conn.commit()

    # Reset Tables
    if RESET_TABLES is True:
        c.execute('''DROP TABLE IF EXISTS covmod''')
        c.execute('''DROP TABLE IF EXISTS names''')
        c.execute('''DROP TABLE IF EXISTS baseprops''')
        c.execute('''DROP TABLE IF EXISTS citation_lookup''')
        c.execute('''DROP TABLE IF EXISTS roles_lookup''')
        c.execute('''DROP TABLE IF EXISTS citations''')
        c.execute('''DROP TABLE IF EXISTS roles''')
        c.execute('''DROP TABLE IF EXISTS modbase''')
        conn.commit()

    # Create Tables
    c.execute('''CREATE TABLE IF NOT EXISTS modbase
                (nameid text PRIMARY KEY NOT NULL,
                 position text,
                 baseid text,
                 formulaid text,
                 cmodid integer,
                 verifiedstatus integer,
                 FOREIGN KEY(baseid) REFERENCES base(baseid) ON DELETE CASCADE ON UPDATE CASCADE,
                 FOREIGN KEY(nameid) REFERENCES names(nameid) ON DELETE CASCADE ON UPDATE CASCADE,
                 FOREIGN KEY(formulaid) REFERENCES baseprops(formulaid) ON DELETE CASCADE ON UPDATE CASCADE,
                 FOREIGN KEY(cmodid) REFERENCES covmod(cmodid) ON DELETE CASCADE ON UPDATE CASCADE)''')

    c.execute('''CREATE TABLE IF NOT EXISTS covmod
                (cmodid integer PRIMARY KEY NOT NULL,
                 symbol text,
                 netcharge text,
                 definition text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS names
                (chebiname text,
                 nameid text PRIMARY KEY NOT NULL,
                 iupacname text,
                 othernames text,
                 smiles text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS baseprops
                (formulaid text PRIMARY KEY NOT NULL,
                 avgmass text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS citations
                (citationid text PRIMARY KEY NOT NULL,
                 title text,
                 pubdate text,
                 authors text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS roles
                (roleid text PRIMARY KEY NOT NULL,
                 role text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS citation_lookup
                (nameid text,
                 citationid text,
                 FOREIGN KEY(nameid) REFERENCES modbase(nameid) ON DELETE CASCADE ON UPDATE CASCADE,
                 FOREIGN KEY(citationid) REFERENCES citations(citationid) ON DELETE CASCADE ON UPDATE CASCADE,
                 PRIMARY KEY(nameid, citationid))''')

    c.execute('''CREATE TABLE IF NOT EXISTS roles_lookup
                (nameid text,
                 roleid text,
                 FOREIGN KEY(nameid) REFERENCES modbase(nameid) ON DELETE CASCADE ON UPDATE CASCADE,
                 FOREIGN KEY(roleid) REFERENCES roles(roleid) ON DELETE CASCADE ON UPDATE CASCADE,
                 PRIMARY KEY(nameid, roleid))''')
    conn.commit()

    # Populate Tables
    formulafill = 0;
    for base in bases:
        childlist = children[base.chebiAsciiName]
        # Retreive childlist of current base
        for child in childlist:
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
                formulafill = formulafill + 1;
                formula = formulafill;

            charge = get_entity(child, CHARGE_SEARCH_STR)
            mass = get_entity(child, MASS_SEARCH_STR)
            citations = concatenate_list(child, CITATION_SEARCH_STR)

            roles = get_roles(child, ONTOLOGY_SEARCH_STR, ONTOLOGY_HAS_ROLE)
            role_names = [role.chebiName for role in roles]
            role_ids = [role.chebiId for role in roles]

            # Populate roles and citations tables with unique data
            for role in range(len(roles)):
                c.execute("SELECT count(*) FROM roles WHERE roleid = ?",
                          (role_ids[role],))
                data = c.fetchone()[0]
                if data == 0:
                    c.execute("UPDATE roles SET role=? WHERE roleid=?",(role_names[role], role_ids[role]))
                    conn.commit()
                    c.execute("INSERT OR IGNORE INTO roles VALUES(?,?)",
                              (role_ids[role], role_names[role]))
                    conn.commit()

            for citation in citations:
                c.execute("SELECT * FROM citations WHERE citationid = ?",
                          [citation])
                data = c.fetchone()
                if data is None:
                    citationinfo = get_full_citation(citation)
                    citationinfo_uni = [info.decode('utf-8') for info in citationinfo]

                    c.execute("UPDATE citations SET title=?, pubdate=?, authors=? WHERE citationid=?",(citationinfo_uni[0], citationinfo_uni[1], citationinfo_uni[2], citation))
                    conn.commit()
                    c.execute("INSERT OR IGNORE INTO citations VALUES(?,?,?,?)",
                              (citation, citationinfo_uni[0], citationinfo_uni[1],
                               citationinfo_uni[2]))
                    conn.commit()

            c.execute("INSERT OR IGNORE INTO baseprops VALUES(?,?)",
                      (str(formula), str(mass)))
            c.execute("INSERT OR IGNORE INTO names VALUES(?,?,?,?,?)",
                      (child.chebiAsciiName, child.chebiId,
                       str(iupac), str(synonyms), str(smiles)))
            c.execute("INSERT OR IGNORE INTO covmod VALUES(NULL,?,?,?)",
                      ('0', str(charge), child.definition))
            conn.commit()
            rowid = c.lastrowid

            c.execute("INSERT OR IGNORE INTO modbase VALUES(?,?,?,?,?,?)",
                      (child.chebiId, '0', base.chebiAsciiName[0],
                       str(formula), rowid, child['verifiedstatus']))
            conn.commit()

            for role in range(len(roles)):
                c.execute("INSERT OR IGNORE INTO roles_lookup VALUES(?,?)",
                          (child.chebiId, role_ids[role]))
            conn.commit()

            for citation in citations:
                c.execute("INSERT OR IGNORE INTO citation_lookup VALUES(?,?)",
                          (child.chebiId, citation))
            conn.commit()

    conn.close()

def create_custom_citations(conn, ref_annots_file_name):
    c = conn.cursor()

    c.execute('''PRAGMA foreign_keys = ON''')

    c.execute('''CREATE TABLE IF NOT EXISTS sequencing_citations
                (nameid text,
                 PMid text,
                 seq_technology text,
                 resolution text,
                 enrichment_method text,
                 FOREIGN KEY(nameid) REFERENCES modbase(nameid) ON DELETE CASCADE ON UPDATE CASCADE)''')
    conn.commit()
    
    with open(ref_annots_file_name) as file:
        reader = csv.reader(file, dialect='excel-tab')
        for line in reader:
            if line[0] == '#':
                continue
            if len(line) > 1:
                stringline = line
                print(stringline)
                c.execute("INSERT OR IGNORE INTO sequencing_citations VALUES(?,?,?,?,?)",
                          (stringline[0], stringline[1], stringline[2], stringline[3], stringline[4]))

    conn.commit()
    conn.close()


def populate_tables(conn, bases, children, client):
    create_base_table(conn, bases)
    create_other_tables(conn, children, bases)
    create_custom_citations(conn, REF_ANNOTS_FULLPATH)


def check_for_duplicates(conn):
    c = conn.cursor()

    c.execute('''PRAGMA foreign_keys = ON''')
    conn.commit()

    for nameid in c.execute('''SELECT nameid FROM modbase'''):
        synonyms = c.execute('''SELECT othernames FROM names WHERE nameid = ?''', nameid)
        for name in synonyms:
            c.execute('''SELECT nameid FROM names WHERE chebiname = ?''', name)
            matchname = c.fetchone()
            print(name, matchname)
            if name == matchname:
                print("Match!")

WHITE_LIST = dnamod_utils.get_list('whitelist')
BLACK_LIST = dnamod_utils.get_list('blacklist')

print("1/4 Searching for bases...")
bases = search_for_bases(client)

print("2/4 Searching for children...")
children = get_children(bases, client)
bases = get_complete_bases(bases, client)

conn = sqlite3.connect(DATABASE_FILE_FULLPATH)

print("3/4 Creating tables...")
populate_tables(conn, bases, children, client)

print("4/4 Finishing up...")
check_for_duplicates(conn)

print("Done!")
