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
# Using Jinja2 as templating engine
import os
import csv
import openbabel
import pybel
import subprocess
from jinja2 import Environment
from jinja2 import FileSystemLoader

BASES = 'Adenine','Thymine','Cytosine','Guanine','Uracil', 'Other Links'
BASE_DICT = {'adenine': 'CHEBI:16708','thymine': 'CHEBI:17821','cytosine': 'CHEBI:16040','guanine': 'CHEBI:16235','uracil': 'CHEBI:17568'}

def render_image(smiles, name):
    pwd = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/images'
    if not os.path.exists(pwd):
        os.makedirs(pwd)
    if not smiles:
        return
    path = pwd + '/' + name + '.svg'
    mol = pybel.readstring('smi', smiles)
    mol.write('svg', path, overwrite=True)


def get_modifications():
    # FSDB command designed to get database info:
    # dbrow '_base_id =~ /c/' -i mod_base.fsdb | dbjoin -i - -i base.fsdb -R base_id | dbjoin -i - -i name.fsdb -R name_id | dbjoin -i - -i cov_modification.fsdb -R cmod_id | dbjoin -i - -i base_properties.fsdb -R property_id > result.fsdb

    # Need to point to this directory: /mnt/work1/users/asood/mordor_home/DNA_Base_Database/
    ans = subprocess.call("dbrow '_base_id =~ /c/' -i mod_base.fsdb | dbjoin -i - -i base.fsdb -R base_id | dbjoin -i - -i name.fsdb -R name_id | dbjoin -i - -i cov_modification.fsdb -R cmod_id | dbjoin -i - -i base_properties.fsdb -R property_id > result2.fsdb",cwd = '/mnt/work1/users/home2/asood/DNA_Base_Database', shell = True)

    if ans == 1:
        print "Error: Reevaluate FSDB pipeline"
        sys.exit()

    # Attempt to set up subprocess pipeline. Not working for some reason :( 
    '''
    ans = subprocess.call(["pwd"])
    print ans
    p1 = subprocess.Popen(["dbrow","'_base_id =~ /c/'", "-i", "mod_base.fsdb"], cwd = '/mnt/work1/users/home2/asood/DNA_Base_Database', stdout = subprocess.PIPE)
    p2 = subprocess.Popen(["dbjoin","-i","-","-i","base.fsdb","-R","base_id"], stdin=p1.stdout, stdout=subprocess.PIPE) # XXX subprocess.PIPE here?
    p1.stdout.close()
    output = p2.communicate()[0]
    print output
    '''


def check_whitelist(molName):
    with open('/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/whitelist/whitelist', 'r') as whitelist:
        reader = csv.reader(whitelist, dialect='excel-tab')
        for line in reader:
            if line == []:
                continue
            if len(line) > 1:
                #print line
                stringline = line[0]
                stringline = stringline.strip()
                flag = line[1]
                #print stringline
                #print flag
                if len(line) > 1:
                    if stringline == molName and flag == '*':
                        return True
        return False


def create_html_pages():
    # Create a Jinja 2 environment object and load in templates
    env = Environment()
    env.loader=FileSystemLoader('/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/Templates')

    page_template = env.get_template('modification.html')
    home_template = env.get_template('homepage.html')

    # Dictionary to store links for hompage
    homepageLinks = {}
    links = []
    blacklist = []

    # result file header:
    # #fsdb -F t property_id cmod_id name_id base_id modbase_id position_location role_id common_name description chebi_name chebi_id iupac_name other_names smiles symbol definition formula net_charge avg_mass
    for BASE in BASES:
        links = []
        with open('/mnt/work1/users/home2/asood/DNA_Base_Database/result.fsdb', 'r') as f1:
            reader = csv.reader(f1, dialect='excel-tab')
            for line in reader:
                if line[0][0] == '#':
                    continue
                #data = line.split('\t')
                data = line
                if data[5] == BASE[0].lower():
                    pwd = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/'
                    if not os.path.exists(pwd):
                        os.makedirs(pwd)

                    writefile = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/'+ data[10] +'.html'
                    f = open(writefile, 'w+')

                    temp = data[13]
                    temp = temp[1:-1]
                    synonyms = temp.split(', ')
                    if synonyms == ['']: synonyms = []

                    temp = data[20]
                    temp = temp[1:-1]
                    citations = temp.split(', ')
                    if citations == ['']: citations = []

                    iupac = data[12]
                    iupac = iupac[1:-1]

                    temp = data[21]
                    temp = temp[1:-1]
                    roles = temp.split(', ')
                    if roles == ['']: roles = []

                    temp = data[22]
                    temp = temp[1:-1]
                    roles_ids = temp.split(', ')
                    if roles_ids == ['']: roles_ids = []

                    render = page_template.render(ChebiName=data[10], Definition=data[16], Formula=data[17], NetCharge=data[18], AverageMass=data[19], IupacName=iupac, Smiles=data[14], Synonyms=synonyms, ChebiId=data[11], CommonName=data[8],Description=data[9], Citations = citations, ParentLink = BASE_DICT[data[8]], Roles = roles, RolesChebi = roles_ids)
                    f.write(render)
                    f.close()

                    link = data[10]
                    #print link
                    if check_whitelist(link):
                        links.append(link)
                    else:
                        blacklist.append(link)

                    smiles  = data[14]
                    smiles = smiles[1:-1]
                    render_image(smiles, data[10])

        homepageLinks[BASE] = links
    homepageLinks['Other Links'] = blacklist
    return homepageLinks

def create_homepage(homepageLinks):
    env = Environment()
    env.loader=FileSystemLoader('/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/Templates')

    page_template = env.get_template('modification.html')
    home_template = env.get_template('homepage.html')

    writefile = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/index.html'
    f = open(writefile, 'w+')
    render = home_template.render(bases=BASES, modifications = homepageLinks)
    f.write(render)
    f.close()


#get_modifications()
links = create_html_pages()
#print links['Adenine']
#print links['Cytosine']
create_homepage(links)
print "Static Site Generated"
