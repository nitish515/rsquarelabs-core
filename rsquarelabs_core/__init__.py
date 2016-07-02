__author__ = 'rrmerugu'


__title__ = 'R2-Core Framework'
__version__ = '0.0.1'
__author__ = 'Ravi RT Merugu'
__license__ = 'Apache License, Version 2.0'
__copyright__ = 'Copyright 2012-2016 www.rsquarelabs.com'

# Version synonym
VERSION = __version__





import os, sys, subprocess, pip, platform
from .config import RSQ_PROJECTS_HOME, RSQ_HOME, RSQ_DB_PATH, RSQ_PROJECTS_CONFIG, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH, RSQ_EXPORT_PATH, RSQ_IMPORT_PATH

if not os.path.exists(RSQ_PROJECTS_HOME):
    os.mkdir(RSQ_PROJECTS_HOME,0755)

if not os.path.exists(RSQ_HOME):
    os.mkdir(RSQ_HOME,0755)

if not os.path.exists(RSQ_SCRIPT_PATH):
    os.mkdir(RSQ_SCRIPT_PATH,0755)

if not os.path.exists(RSQ_EXPORT_PATH):
    os.mkdir(RSQ_EXPORT_PATH,0755)

if not os.path.exists(RSQ_IMPORT_PATH):
    os.mkdir(RSQ_IMPORT_PATH,0755)

if not os.path.exists(RSQ_BACKUP_PATH):
    os.mkdir(RSQ_BACKUP_PATH,0755)

if not os.path.exists(RSQ_PROJECTS_CONFIG): # not very much needed
    os.mkdir(RSQ_PROJECTS_CONFIG, 0755)

# Checking the python version for gromacs.
if not sys.version_info[:2] == (2, 7):
    print ("Sorry! This tool works only with python version 2.7.")
    exit()

# Checking the installation of gromacs.
if not subprocess.call("gmx", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
    print "Please! Install gromacs version 5.x"
    print  "For more info: https://github.com/rsquarelabs/rsquarelabs-core/wiki"
    exit()

# Checking if the requirements are installed.
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version)
     for i in installed_packages])

required_packages = ['requests>=2.10.0', 'termcolor>=1.1.0', 'pip>=6.1.2']

for package in installed_packages_list:
    if package in required_packages:
        required_packages.remove(package)

# if not required_packages==[]:
#     print "Please! Reinstall the gromacs."
#     exit()

# Checking the operating system for gromacs.
if platform.system()=='Windows':
    print "Installation of gromacs is denied in this operating system"
    exit()


