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

SYNONYM_SEARCH_STR = 'Synonyms'
IUPAC_SEARCH_STR = 'IupacNames'
SMILES_SEARCH_STR = 'smiles'
FORMULA_SEARCH_STR = 'Formulae'
CHARGE_SEARCH_STR = 'charge'
MASS_SEARCH_STR = 'mass'

# Using Suds web services client for soap
from suds.client import Client
from sys import maxint

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
    dnaBases = ['cytosine', 'thymine', 'adenine', 'guanine', 'uracil']

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
    with open("base.fsdb", "w+") as basefile:
        basefile.write('#fsdb -F S base_id common_name description\n')
        for base in bases:
            basefile.write(base.chebiAsciiName[0] + '  ' + base.chebiAsciiName
                           + '  ' + base.definition + '\n')


def concatenate_name(child, attribute):
    return ([synonym_data.data
             for synonym_data in getattr(child, attribute)]
            if attribute in dir(child) else [])


def get_entity(child, attribute):
    return ([getattr(child, attribute)]
            if attribute in dir(child) else [])


def create_other_tables(children, bases):
    modbase = open("mod_base.fsdb", "w+")
    modbase.write("#fsdb -F S modbase_id position_location base_id name_id\
                  property_id role_id cmod_id\n")
    covmod = open("cov_modification.fsdb", "w+")
    covmod.write("#fsdb -F S cmod_id symbol definition\n")
    name = open("name.fsdb", "w+")
    name.write("#fsdb -F S name_id chebi_name chebi_id iupac_name other_names\
               smiles\n")
    baseprop = open("base_properties.fsdb", "w+")
    baseprop.write("#fsdb -F S property_id formula net_charge avg_mass\n")

    id_counter = 0
    for base in bases:
        childlist = children[base.chebiAsciiName]
        for child in childlist:
            modbase.write(('modbase'+str(id_counter)) + '  ' + '0' + '  ' +
                          base.chebiAsciiName[0] + '  ' +
                          ('name'+str(id_counter)) + '  ' +
                          ('prop'+str(id_counter)) + '  ' + '0' + '  ' +
                          ('covmod'+str(id_counter)) + '\n')

            covmod.write(('covmod'+str(id_counter)) + '  ' + '0' + '  ' +
                         str(child.definition) + '\n')

            synonyms = concatenate_name(child, SYNONYM_SEARCH_STR)
            iupac = concatenate_name(child, IUPAC_SEARCH_STR)
            smiles = get_entity(child, SMILES_SEARCH_STR)
            formula = concatenate_name(child, FORMULA_SEARCH_STR)
            charge = get_entity(child, CHARGE_SEARCH_STR)
            mass = get_entity(child, MASS_SEARCH_STR)

            name.write('name'+str(id_counter) + '  ' + child.chebiAsciiName
                       + '  ' + child.chebiId + '  ' + str(iupac) + '  ' +
                       str(synonyms) + '  ' + str(smiles) + '\n')

            baseprop.write('property'+str(id_counter) + '  ' +
                           str(formula) + '  ' + str(charge) + '  ' +
                           str(mass) + '\n')

            id_counter = id_counter + 1

    modbase.close()
    covmod.close()
    name.close()
    baseprop.close()


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
