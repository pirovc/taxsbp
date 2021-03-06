import io
import os
import re

from setuptools import find_packages
from setuptools import setup

def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())

setup(
    name="taxsbp",
    version="1.1.1",
    url="https://www.github.com/pirovc/taxsbp",
    license='MIT',

    author="Vitor C. Piro",
    author_email="pirovc@posteo.net",

    description="Implementation of the approximation algorithm for the hierarchically structured bin packing problem adapted for the NCBI Taxonomy database.",
    long_description=read("README.md"),

    packages=['taxsbp'],
    install_requires=['binpacking==1.4.3'],

    entry_points = {'console_scripts': ['taxsbp=taxsbp.taxsbp:main_cli']},

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
