# DNAmod #

The DNAmod webpage is: [https://DNAmod.hoffmanlab.org](http://DNAmod.hoffmanlab.org) (which redirects to [https://www.pmgenomics.ca/hoffmanlab/proj/dnamod](https://www.pmgenomics.ca/hoffmanlab/proj/dnamod)).

This outlines the steps needed to set up a Python environment to reproduce or run a local DNAmod instance.

## Summary ##

Covalent DNA modifications have been found in numerous organisms and more are continually being discovered and characterized, as detection methods improve. Many of these modifications can affect the conformation of the DNA double helix, often resulting in downstream effects upon transcription factor binding. Some of these modifications have been demonstrated to be stable, while others are viewed as merely transient.

DNAmod catalogs information on known DNA modifications. It aims to profile modifications' properties, with diagrams and links to databases such as [ChEBI](http://www.ebi.ac.uk/chebi). It also captures citations to published literature on these modifications. This makes it easier for those new to the field to explore the array of potential modifications.

DNAmod is comprised of this static website and a backing [SQLite](https://www.sqlite.org/) database. The database is created using Python, including the SOAP client [suds](https://fedorahosted.org/suds/) and [Biopython](http://biopython.org/wiki/Main_Page). The website is also created using Python, makes use of [Open Babel](http://openbabel.org/) via its Python wrapper, [Pybel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html), and uses the [Jinja2](http://jinja.pocoo.org/) templating engine.

This project is open source and updated regularly. We also manually curate the modifications verified to occur *in vivo*.

## Installation ##

To set up ensure that Python is installed on your machine along with the dependencies listed below.

```DNA_mod_database.db``` contains the SQLite database for DNAmod.

```populate_database_sql.py``` pulls data from ChEBI and imports the data into the SQLite database.

```create_mod_staticsite_sql.py``` takes the SQLite database and pushes it through the Jinja2 templates located in the templates folder to create the static site. The location of the static site is the ```static``` folder.

The whitelist and blacklist located in the ```whitelist``` folder control the modifications which show up as verified on the finished site. These can be reconfigured depending on your preferences.

The ```sequencing``` file allows the user to add custom citations to the database and associate these citations to particular modifications. Changing the format of this file may also warrant modification of the Jinja2 template ```modification.html```.

To update and recreate the database one would usually run ```populate_database_sql.py``` followed by ```create_mod_staticsite_sql.py``` however running the shell script ```update_dnamod.sh``` will complete those steps. The shell script ```sync_live_site.sh``` can be reconfigured to allow you to quickly update your live site location for DNAmod. 

### Dependencies ###

The following dependencies, in addition to Python (≥ version 2.7.6 necessary for some systems), are needed to build DNAmod.

1. [OpenBabel](http://openbabel.org/wiki/Category:Installation)

2. [SQLite](https://www.sqlite.org/) ≥ version 3.8.3

#### PyPI packages ####

The following packages can be obtained from [The Python Package Index, PyPI](https://pypi.python.org/pypi) and installed with ```pip```.

1. [Biopython](http://biopython.org/wiki/Main_Page): ```biopython``` on PyPI

2. [```Jinja2```](http://jinja.pocoo.org/)

3. [Pybel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html): ```openbabel``` on PyPI

4. [```suds```](https://fedorahosted.org/suds/) (a Python SOAP client)

5. [```unicodecsv```](https://pypi.python.org/pypi/unicodecsv)

6. [```pysqlite2```](https://pypi.python.org/pypi/pysqlite)

## License ##

Released under the [GNU General Purpose License, version 2 (GPLv2)](http://www.gnu.org/licenses/gpl-2.0-standalone.html).

## Contacts ##

This project was created at the [Hoffman Lab](https://www.pmgenomics.ca/hoffmanlab/) by Ankur Jai Sood, Coby Viner, and Michael M. Hoffman.

There is a moderated [dnamod-announce](https://listserv.utoronto.ca/cgi-bin/wa?A0=DNAMOD-ANNOUNCE-L&X=E5FDFD12D6CD9E97CC&Y) mailing list that you can subscribe to for information on new releases of DNAmod.

There is also a [dnamod-users](https://listserv.utoronto.ca/cgi-bin/wa?A0=DNAMOD-L&X=E5FDFD12D6CD9E97CC&Y) mailing list for general discussion and questions about the use of DNAmod.