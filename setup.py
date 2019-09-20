#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup

with open(os.path.join('cdenrichrgenestoterm', '__init__.py')) as ver_file:
    for line in ver_file:
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'gseapy'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='cdenrichrgenestoterm',
    version=version,
    description="Python Boilerplate contains all the boilerplate you need to create a Python package.",
    long_description=readme + '\n\n' + history,
    author="Christopher Churas",
    author_email='churas.camera@gmail.com',
    url='https://github.com/ndexbio/cdenrichrgenestoterm',
    packages=[
        'cdenrichrgenestoterm',
    ],
    package_dir={'cdenrichrgenestoterm':
                 'cdenrichrgenestoterm'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='cdenrichrgenestoterm',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    scripts=['cdenrichrgenestoterm/cdenrichrgenestoterm.py'],
    test_suite='tests',
    tests_require=test_requirements
)
