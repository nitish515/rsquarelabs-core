__author__ = 'rrmerugu'

import os, sys, subprocess, pip, platform
from .config import RSQ_PROJECTS_HOME, RSQ_HOME, RSQ_DB_PATH, RSQ_PROJECTS_CONFIG

if not os.path.exists(RSQ_PROJECTS_HOME):
    os.mkdir(RSQ_PROJECTS_HOME,0755)

if not os.path.exists(RSQ_HOME):
    os.mkdir(RSQ_HOME,0755)

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

required_packages = ['requests==2.10.0', 'termcolor==1.1.0', 'pip==8.1.2']

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


