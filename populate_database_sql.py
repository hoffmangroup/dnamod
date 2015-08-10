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
import os
import sqlite3
from suds.client import Client
from sys import maxint
from Bio import Entrez

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
DNA_BASES = ['cytosine', 'thymine', 'adenine', 'guanine', 'uracil']
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
print FILE_PATH

# Url to CHEBI wsdl schema
url = 'https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'

# Create a client object from above url
client = Client(url)

# This line prints out all methods and data types available from webservice
# print client


def search_exact_match(searchKey, client):

    # Parameters for getLiteEntity() search
    searchCategory = 'ALL'
    maximumResults = maxint
    starsCategory = 'THREE ONLY'

    # Query CHEBI using above key and store in list results
    # getLiteEntity(xs:string search, SearchCategory searchCategory,
    #               xs:int maximumResults, StarsCategory stars)

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
    for elements in range(len(results)):
        if results[elements].chebiAsciiName == searchKey:
            result = results[elements]

    return result


def search_for_bases(client):
    # Initialize DNA bases
    dnaBases = DNA_BASES

    # Code for debugging
    # print dnaBases[0]
    # print len(dnaBases)
    # print range(len(dnaBases))

    # Initialize empy list
    result = []

    # Search CHEBI for bases and return lite entries
    for elements in range(len(dnaBases)):
        # print elements # Output for debugging
        result.append(search_exact_match(dnaBases[elements], client))

    return result


def filter_stars(content, stars, client):
    result = []
    for elements in range(len(content)):
        # print content[elements] # For Debugging
        # content[elements] = search_exact_match(content[elements].chebiName,
        #                                        client)
        temphold = content[elements]
        # print temphold # For debugging
        content[elements] = get_complete_entity(content[elements].chebiId,
                                                client)
        # print content[elements] # For Debugging
        if content[elements] == 'Invalid Input' or content[elements] == 'DNE':
            continue
        elif (content[elements].entityStar == stars and
              temphold.type == 'has functional parent'):
            result.append(content[elements])
    return result


def get_children(bases, client):
    modDictionary = {}
    for elements in range(len(bases)):
        result = client.service.getOntologyChildren(bases[elements].chebiId)
        result = result.ListElement
        result = filter_stars(result, 3, client)
        # print bases[elements] # Debugging output
        modDictionary[bases[elements].chebiAsciiName] = result
    return modDictionary


def get_complete_entity(CHEBIid, client):
    result = client.service.getCompleteEntity(CHEBIid)
    return result


def get_complete_bases(bases, client):
    result = []
    for elements in range(len(bases)):
        result.append(get_complete_entity(bases[elements].chebiId, client))
    return result


def create_base_table(bases):
    conn = sqlite3.connect('DNA_mod_database.db')
    c = conn.cursor()

    c.execute('''DROP TABLE IF EXISTS base''')
    c.execute('''CREATE TABLE IF NOT EXISTS base
                (baseid TEXT PRIMARY KEY NOT NULL,
                 commonname text,
                 basedefinition text)''')
    conn.commit()

    for base in bases:
        c.execute("INSERT INTO base VALUES(?,?,?)",
                  (base.chebiAsciiName[0],
                   base.chebiAsciiName, base.definition))
    conn.commit()
    conn.close()


def concatenate_list(child, attribute):
    return ([synonym_data.data
             for synonym_data in getattr(child, attribute)]
            if attribute in dir(child) else [])


def get_entity(child, attribute):
    return ([getattr(child, attribute)]
            if attribute in dir(child) else [])


def get_roles(child, attribute, selector):
    return ([ontologyObject
             for ontologyObject in getattr(child, attribute)
            if ontologyObject.type == ONTOLOGY_HAS_ROLE])


def get_full_citation(PMID):
    result = []
    # isbook = False # Unsed at the moment
    isarticle = False
    Entrez.email = "jai.sood@hotmail.com"
    handle = Entrez.efetch("pubmed", id=PMID, retmode="xml")
    records = Entrez.parse(handle)

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

    if articleTitle is not None and articleTitle != []:
        result.append(articleTitle)
    else:
        result.append('Title Not Found')

    if isarticle:
        if publicationDate is not None and publicationDate != []:
            date = (publicationDate[0]['Month'] + '-' +
                    publicationDate[0]['Day'] + '-' +
                    publicationDate[0]['Year'])
            result.append(str(date))
        else:
            result.append('Publication Date Not Found')
    else:
        if publicationDate is not None and publicationDate != []:
            date = (publicationDate['Month'] + '-' + publicationDate['Day']
                    + '-' + publicationDate['Year'])
            result.append(date)
        else:
            result.append('Publication Date Not Found')

    if authors is not None and authors != []:
        name = (authors[0]['LastName'] + ', ' + authors[0]['Initials'] +
                '. et al.')
        result.append(name)
    else:
        result.append('Author Not Found')

    return result


def create_other_tables(children, bases):
    conn = sqlite3.connect('DNA_mod_database.db')
    c = conn.cursor()

    # Reset Tables
    if RESET_TABLES is True:
        c.execute('''DROP TABLE IF EXISTS modbase''')
        c.execute('''DROP TABLE IF EXISTS covmod''')
        c.execute('''DROP TABLE IF EXISTS names''')
        c.execute('''DROP TABLE IF EXISTS baseprops''')
        c.execute('''DROP TABLE IF EXISTS citations''')
        c.execute('''DROP TABLE IF EXISTS roles''')
        c.execute('''DROP TABLE IF EXISTS citation_lookup''')
        c.execute('''DROP TABLE IF EXISTS roles_lookup''')
        conn.commit()

    # Create Tables
    c.execute('''CREATE TABLE IF NOT EXISTS modbase
                (modbaseid INTEGER PRIMARY KEY NOT NULL,
                 position text,
                 baseid text,
                 nameid integer,
                 propertyid integer,
                 citationid text,
                 cmodid integer,
                 roleid integer,
                 FOREIGN KEY(nameid) REFERENCES names(nameid),
                 FOREIGN KEY(propertyid) REFERENCES baseprops(propertyid),
                 FOREIGN KEY(cmodid) REFERENCES covmod(cmodid),
                 FOREIGN KEY(citationid) REFERENCES citations(citationid),
                 FOREIGN KEY(roleid) REFERENCES roles(roleid))''')

    c.execute('''CREATE TABLE IF NOT EXISTS covmod
                (cmodid INTEGER PRIMARY KEY NOT NULL,
                 symbol text,
                 definition text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS names
                (nameid INTEGER PRIMARY KEY NOT NULL,
                 chebiname text,
                 chebiid text,
                 iupacname text,
                 othernames text,
                 smiles text)''')

    c.execute('''CREATE TABLE IF NOT EXISTS baseprops
                (propertyid INTEGER PRIMARY KEY NOT NULL,
                 formula text,
                 netcharge text,
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
                (modid int,
                 citationid text,
                 FOREIGN KEY(modid) REFERENCES modbase(modbaseid),
                 FOREIGN KEY(citationid) REFERENCES citations(citationid),
                 PRIMARY KEY(modid, citationid))''')

    c.execute('''CREATE TABLE IF NOT EXISTS roles_lookup
                (modid int,
                 roleid text,
                 FOREIGN KEY(modid) REFERENCES modbase(modbaseid),
                 FOREIGN KEY(roleid) REFERENCES roles(roleid),
                 PRIMARY KEY(modid, roleid))''')
    conn.commit()

    # Populate Tables
    for base in bases:
        childlist = children[base.chebiAsciiName]
        # Retreive childlist of current base
        for child in childlist:
            # Parse CHEBI datastructure for relevant info
            synonyms = concatenate_list(child, SYNONYM_SEARCH_STR)
            iupac = concatenate_list(child, IUPAC_SEARCH_STR)
            smiles = get_entity(child, SMILES_SEARCH_STR)
            formula = concatenate_list(child, FORMULA_SEARCH_STR)
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
                    c.execute("INSERT INTO roles VALUES(?,?)",
                              (role_ids[role], role_names[role]))

            for citation in citations:
                c.execute("SELECT * FROM citations WHERE citationid = ?",
                          [citation])
                data = c.fetchone()
                if data is None:
                    citationinfo = get_full_citation(citation)
                    c.execute("INSERT INTO citations VALUES(?,?,?,?)",
                              (citation, citationinfo[0], citationinfo[1],
                               citationinfo[2]))

            c.execute("INSERT INTO baseprops VALUES(NULL,?,?,?)",
                      (str(formula), str(charge), str(mass)))
            c.execute("INSERT INTO names VALUES(NULL,?,?,?,?,?)",
                      (child.chebiAsciiName, child.chebiId,
                       str(iupac), str(synonyms), str(smiles)))
            c.execute("INSERT INTO covmod VALUES(NULL,?,?)",
                      ('0', child.definition))
            rowid = c.lastrowid

            for role in range(len(roles)):
                c.execute("INSERT INTO roles_lookup VALUES(?,?)",
                          (rowid, role_ids[role]))
            for citation in citations:
                c.execute("INSERT INTO citation_lookup VALUES(?,?)",
                          (rowid, citation))

            c.execute("INSERT INTO modbase VALUES(NULL,?,?,?,?,?,?,?)",
                      ('0', base.chebiAsciiName[0], rowid,
                       rowid, rowid, rowid, rowid))
            conn.commit()

    conn.close()


def populate_tables(bases, children, client):
    create_base_table(bases)
    create_other_tables(children, bases)


# 1. Performs search of CHEBI database for DNA bases
# 2. Returns chebiId and chebiAsciiName of bases
print "1/4 Searching for bases..."
bases = search_for_bases(client)
# print bases # Debugging Output

# 3. Searches CHEBI database for all entities of which
#    the DNA bases are functional parents and are rated three or more stars
print "2/4 Searching for children..."
children = get_children(bases, client)
bases = get_complete_bases(bases, client)
# print children.keys() # Debugging output
# print bases[0] # Debugging output
# print children[bases[0].chebiAsciiName][0]  # Debugging output
# print children[bases[0].chebiAsciiName]  # Debugging output
# print children.items()

# 4. Using above results, populates DNA post-transciptional modification table
print "3/4 Creating tables..."
populate_tables(bases, children, client)
print "4/4 Finishing up..."
print "Done!"
