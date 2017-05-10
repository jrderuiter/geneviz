#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    README = readme_file.read()

with open('HISTORY.rst') as history_file:
    HISTORY = history_file.read()

REQUIREMENTS = [
    'future', 'pandas', 'matplotlib', 'pysam', 'toolz', 'intervaltree',
    'seaborn', 'pybiomart'
]

EXTRAS_REQUIRE = {
    'dev': [
        'sphinx', 'sphinx-autobuild', 'sphinx-rtd-theme', 'bumpversion',
        'pytest>=2.7', 'pytest-mock', 'pytest-helpers-namespace', 'pytest-cov',
        'python-coveralls'
    ]
}

setup(
    name='geneviz',
    version='0.1.0',
    description="Python library for creating visualizations of genomic data.",
    long_description=README + '\n\n' + HISTORY,
    author="Julian de Ruiter",
    author_email='julianderuiter@gmail.com',
    url='https://github.com/jrderuiter/geneviz',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=REQUIREMENTS,
    extras_require=EXTRAS_REQUIRE,
    license="MIT license",
    zip_safe=False,
    keywords='geneviz',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ])
