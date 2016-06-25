#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'rrmerugu'

import os, re, sys, shutil
from distutils.core import setup
from setuptools import find_packages

# readme_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.md')
# print readme_file
# readme = open(readme_file).read()






def get_version(package):
    """
    Return package version as listed in `__version__` in `init.py`.
    """
    init_py = open(os.path.join(package, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)


version = get_version('rsquarelabs_core')

if sys.argv[-1] == 'publish':
    try:
        import pypandoc
    except ImportError:
        print("pypandoc not installed.\nUse `pip install pypandoc`.\nExiting.")
    if os.system("pip freeze | grep wheel"):
        print("wheel not installed.\nUse `pip install wheel`.\nExiting.")
        sys.exit()
    if os.system("pip freeze | grep twine"):
        print("twine not installed.\nUse `pip install twine`.\nExiting.")
        sys.exit()
    os.system("python setup.py sdist bdist_wheel")
    os.system("twine upload dist/*")
    print("You probably want to also tag the version now:")
    print("  git tag -a %s -m 'version %s'" % (version, version))
    print("  git push --tags")
    shutil.rmtree('dist')
    shutil.rmtree('build')
    shutil.rmtree('djangorestframework.egg-info')
    sys.exit()

readme = "This is the library of automation pipeline modules developed at RSQUARE LABS."

github_url = 'http://github.com/rsquarelabs/rsquarelabs-core'
version = "0.0.7dev"

setup(name='rsquarelabs-core',
version= version,
description='This is the library of automation pipeline modules developed at RSQUARE LABS.',
long_description= readme,
author='Ravi RT Merugu',
author_email='rrmerugu@gmail.com',
url = github_url,
packages = find_packages(),
package_data={'rsquarelabs_core' : ['*']},
install_requires=['bottle','termcolor','requests'],
download_url='%s/tarball/%s' %(github_url,version ),
keywords = ['Computational Biology', 'Molecular Modelling', 'Bioinformatics', 'Automation'])