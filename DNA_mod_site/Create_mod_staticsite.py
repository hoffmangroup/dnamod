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
import subprocess
from jinja2 import Environment
from jinja2 import FileSystemLoader

BASES = 'Adenine','Thymine','Cytosine','Guanine','Uracil'

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


def create_html_pages():
    # Create a Jinja 2 environment object and load in templates
    env = Environment()
    env.loader=FileSystemLoader('/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/Templates')

    page_template = env.get_template('modification.html')
    home_template = env.get_template('homepage.html')

    # Dictionary to store links for hompage
    homepageLinks = {}
    links = []

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
                if data[3] == BASE[0].lower():
                    pwd = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/'
                    if not os.path.exists(pwd):
                        os.makedirs(pwd)

                    writefile = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/'+ data[9] +'.html'
                    f = open(writefile, 'w+')
                    temp = data[12]
                    temp = temp[1:-1]
                    synonyms = temp.split(', ')
                    iupac = data[11]
                    iupac = iupac[1:-1]
                    render = page_template.render(ChebiName=data[9], Definition=data[15], Formula=data[16], NetCharge=data[17], AverageMass=data[18], IupacName=iupac, Smiles=data[13], Synonyms=synonyms, ChebiId=data[10], CommonName=data[7],Description=data[8])
                    f.write(render)
                    f.close()
                    link = data[9]
                    links.append(link)

        homepageLinks[BASE] = links
    return homepageLinks

def create_homepage(homepageLinks):
    env = Environment()
    env.loader=FileSystemLoader('/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/Templates')

    page_template = env.get_template('modification.html')
    home_template = env.get_template('homepage.html')

    writefile = '/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static/homepage.html'
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
