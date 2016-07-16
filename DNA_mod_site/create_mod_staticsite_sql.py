#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

'''
Ankur Jai Sood
09/6/2015

Create_mod_staticsite.py
Function:
1. Opens FSDB database files
2. Returns a modification
3. Creates a html page for the modifcation based on created jinja2 template
4. Creates a home page with links to html modification pages
'''
import codecs
from collections import OrderedDict
from itertools import izip
import os
import pybel
import sqlite3
import sys

# Using Jinja2 as templating engine
from jinja2 import Environment
from jinja2 import FileSystemLoader

# permit import from parent directory
sys.path.append(os.path.join(sys.path[0], '..'))
import dnamod_utils

# Program Constants
ENCODING = 'utf8'
BASES = 'Adenine', 'Cytosine', 'Guanine', 'Thymine', 'Uracil'
VERIFIED_BASES = 'Adenine', 'Cytosine', 'Guanine', 'Thymine'
UNVERIFIED_BASES = ('UnverifiedAdenine', 'UnverifiedThymine',
                    'UnverifiedCytosine', 'UnverifiedGuanine',
                    'UnverifiedUracil')
BASE_DICT = {'adenine': 'CHEBI:16708', 'thymine': 'CHEBI:17821',
             'cytosine': 'CHEBI:16040', 'guanine': 'CHEBI:16235',
             'uracil': 'CHEBI:17568'}
REF_COL_NAMES = ['citationid', 'title', 'pubdate', 'authors']

HTML_FILES_DIR = dnamod_utils.get_constant('site_html_dir')
TEMPLATE_DIR = dnamod_utils.get_constant('site_template_dir')

# TODO move to constants.sh ...
SEQ_ANNOT_TABLE = 'sequencing_citations'
NATURE_ANNOT_TABLE = 'nucleobase_nature_info'
REFERENCES_TABLE = 'citations'

IMAGE_FORMAT = 'svg'


def render_image(smiles, name):
    """Creates images of chemical structures from their SMILES, with PyBel.
       The SVG image is saved. The path of the image is returned, if it exists.
       Otherwise, None is returned.
    """

    img_dir = dnamod_utils.get_constant('site_image_dir')

    if not os.path.exists(img_dir):
        os.makedirs(img_dir)
    if not smiles:
        return

    image_path = os.path.join(img_dir, "{}.{}".format(name, IMAGE_FORMAT))

    mol = pybel.readstring('smi', smiles)
    mol.write(IMAGE_FORMAT, image_path, overwrite=True)

    if not os.path.isfile(image_path):
        print("Warning: failed to render image"
              "for {}".format(name), file=sys.stderr)
        image_path = None

    return image_path


def get_citations(lookup_key, cursor):
    c = cursor.cursor()
    citationList = []
    c.execute("SELECT * FROM citation_lookup WHERE nameid = ?", (lookup_key,))
    results = c.fetchall()
    for row in results:
        citationid = row[1]
        c.execute("SELECT * FROM citations WHERE citationid = ?",
                  (citationid,))
        query = c.fetchone()

        citationList.append(dict(izip(REF_COL_NAMES,
                            [item.encode(ENCODING) for item in query])))
    return citationList


def get_table_headers(cursor, table_name):
    c = cursor.cursor()
    c.execute("PRAGMA table_info({})".format(table_name))
    header = [result[1] for result in c.fetchall()]
    return list(header)


def get_expanded_alphabet(id, cursor):
    c = cursor.cursor()

    exp_alph_dict = {}

    exp_alph_headers = get_table_headers(cursor, 'expanded_alphabet')

    # TODO refactor to dynamically select columns from header length
    # seq_headers contains the header for this table
    # overall orders first by date, but still grouped by method
    c.execute('''SELECT DISTINCT nameid,
                 [{1}], [{2}], [{3}], [{4}], [{5}]
                 FROM expanded_alphabet
                 WHERE nameid = ?
                 '''.format(*exp_alph_headers),
              (id,))
    results = c.fetchall()

    for result in results:
        exp_alph_dict.update(dict(izip(exp_alph_headers, result)))

    return exp_alph_dict


def get_mod_base_ref_annot_data(id, cursor, table):
    """Get data for mod. base annotations with references.

    Keyword arguments:
    id -- the entries ChEBI ID
    cursor -- the SQLite cursor
    table -- the table containing the annotations

    Returns:
    A dictionary containing the headers as keys and
    a list of each data row per header as values.
    """

    c = cursor.cursor()

    annot_dict_list = []

    table_header = get_table_headers(cursor, table)

    # use all columns except first (id) and last (ref.)
    table_header.pop(0)
    reference_col_name = table_header.pop(-1)

    subquery_alias = 'modbase_refs_ord'

    sel_cols_str = ''
    for num, col in enumerate(table_header, 1):
        sel_cols_str += "{}.[{}]".format(subquery_alias, col)

        if num < len(table_header):
            sel_cols_str += ", "
        else:
            break

    # Overall orders first by reference's date,
    # both within groups and to determine group order,
    # but still grouped by the second column of data.
    c.execute('''SELECT DISTINCT
                     GROUP_CONCAT({5}.citationid, ';'),
                     GROUP_CONCAT({5}.title, ';'),
                     GROUP_CONCAT({5}.pubdate, ';'),
                     GROUP_CONCAT({5}.authors, ';'),
                     {2}
                 FROM (
                        SELECT * FROM {0}
                        JOIN {1} AS ref ON {0}.[{3}]
                            LIKE '%' || ref.citationid || '%'
                        WHERE nameid = ?
                        ORDER BY COALESCE(date(ref.pubdate),
                                         ref.authors, 1)
                 ) {5}
                 GROUP BY {2}
                 ORDER BY COALESCE({5}.[{4}],
                                   date({5}.pubdate),
                                   {5}.authors, 1)
                 '''.format(table, REFERENCES_TABLE, sel_cols_str,
                            reference_col_name, table_header[0],
                            subquery_alias),
              (id,))

    results = c.fetchall()

    for result in results:
        annot_dict_list += [OrderedDict(izip(REF_COL_NAMES +
                                        table_header, result))]
    return annot_dict_list


def create_html_pages():
    # Load in SQLite database
    conn = sqlite3.connect(dnamod_utils.get_constant('database'))
    c = conn.cursor()

    # Create a Jinja 2 environment object and load in templates
    env = Environment()
    env.loader = FileSystemLoader(TEMPLATE_DIR)

    page_template = env.get_template('modification.html')

    # Dictionary to store links for hompage
    homepageLinks = {}
    links = []
    blacklist = []

    c.execute('''DROP TABLE IF EXISTS temp''')
    c.execute('''CREATE TABLE temp AS SELECT * FROM
                    (SELECT * from
                        (SELECT * FROM
                            (SELECT * FROM modbase
                             AS MB NATURAL JOIN baseprops)
                        AS MB_BP NATURAL JOIN covmod)
                    AS MB_CV NATURAL JOIN names)
                AS MB_B NATURAL JOIN base''')
    conn.commit()

    if not os.path.exists(HTML_FILES_DIR):
        os.makedirs(HTML_FILES_DIR)

    for BASE in BASES:
        links = []
        blacklist = []

        baseid = BASE[0].lower()
        mods = c.execute("SELECT * FROM temp WHERE baseid = ?", baseid)
        conn.commit()

        for mod in mods:
            # Read data:
            formula = mod[3]
            netcharge = mod[8]
            avgmass = mod[6]
            definition = mod[9]
            chebiname = mod[10].encode('ascii')
            chebiid = mod[0]
            iupacname = mod[11]
            synonyms = mod[12]
            smiles = mod[13].encode('ascii')
            commonname = mod[14]

            citation_lookup = mod[0]

            # XXX TODO cleanup commented-out code
            # roles_lookup = mod[7] # unused as roles are not on site

            citations = get_citations(citation_lookup, conn)

            roles = []
            roles_ids = []

            print("Creating page: " + chebiname + ".html")
            # Process SMILES to render image
            smiles = smiles[1:-1]

            image_path = render_image(smiles, chebiname)

            # Process synonyms for list
            synonyms = synonyms[1:-1]
            synonyms = synonyms.split(', ')
            if synonyms == ['']:
                synonyms = []
            result = []
            for name in synonyms:
                name = name[1:-1]
                result.append(name)
            synonyms = result

            smiles = smiles.decode('ascii')
            chebiname = chebiname.decode('ascii')

            # Formatting
            formula = formula[1:-1]
            netcharge = netcharge[1:-1]
            iupacname = iupacname[1:-1]
            avgmass = avgmass[1:-1]

            # Write html page
            writefile = os.path.join(HTML_FILES_DIR, chebiname + '.html')
            f = codecs.open(writefile, 'w+', encoding=ENCODING)

            for key in REF_COL_NAMES:  # decode encoded values
                for citation in citations:
                    citation[key] = citation[key].decode(ENCODING)

            seq_annot = get_mod_base_ref_annot_data(citation_lookup,
                                                    conn, SEQ_ANNOT_TABLE)

            nature_annot = get_mod_base_ref_annot_data(citation_lookup,
                                                       conn,
                                                       NATURE_ANNOT_TABLE)

            # TODO revise second name
            ref_annot_tab_names = ['Mapping Techniques', 'Nature']

            ref_annots = [seq_annot, nature_annot]

            expandedalpha = get_expanded_alphabet(chebiid, conn)

            render = page_template.render(ChebiName=chebiname,
                                          Definition=definition,
                                          Formula=formula,
                                          NetCharge=netcharge,
                                          AverageMass=avgmass,
                                          IupacName=iupacname,
                                          Smiles=smiles,
                                          Synonyms=synonyms,
                                          ChebiId=chebiid,
                                          CommonName=commonname,
                                          Image=image_path,
                                          Citations=citations,
                                          ParentLink=BASE_DICT[commonname],
                                          Roles=roles,
                                          RolesChebi=roles_ids,
                                          RefAnnotTabNames=ref_annot_tab_names,
                                          RefAnnots=ref_annots,
                                          RefAnnotsRefColNames=REF_COL_NAMES,
                                          # pass ExpandedAlpha=None to disable
                                          ExpandedAlpha=expandedalpha)
            f.write(render)
            f.close()

            # Check if mod is on whitelist
            link = chebiname
            if mod[5]:
                links.append(link)
            else:
                blacklist.append(link)

        links = sorted(links, key=lambda s: s.lower())
        homepageLinks[BASE] = links
        blacklist = sorted(blacklist, key=lambda s: s.lower())
        blacklistBase = 'Unverified' + BASE
        homepageLinks[blacklistBase] = blacklist
    return homepageLinks


def create_homepage(homepageLinks):
    verifiedBases = {}
    unverifiedBases = {}
    for base in BASES:
        verifiedBases[base] = homepageLinks[base]
        unvername = 'Unverified' + base
        unverifiedBases[base] = homepageLinks[unvername]

    verifiedBases['Thymine'] = (verifiedBases['Thymine'] +
                                verifiedBases['Uracil'])
    del(verifiedBases['Uracil'])

    env = Environment()
    env.loader = FileSystemLoader(TEMPLATE_DIR)

    home_template = env.get_template('homepage.html')

    writefile = os.path.join(HTML_FILES_DIR, 'index.html')

    f = codecs.open(writefile, 'w+', encoding=ENCODING)

    render = home_template.render(bases=VERIFIED_BASES,
                                  modifications=verifiedBases,
                                  unverifiedbases=BASES,
                                  unverifiedmodifications=unverifiedBases)
    f.write(render)
    f.close()

print("Generating Static Site....")
links = create_html_pages()
create_homepage(links)
print("Static Site Generated")
