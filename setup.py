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
    version="0.1.2",
    url="https://www.github.com/pirovc/taxsbp",
    license='MIT',

    author="Vitor C. Piro",
    author_email="vitorpiro@gmail.com",

    description="Implementation of the approximation algorithm for the hierarchically structured bin packing problem based on the NCBI Taxonomy database.",
    long_description=read("README.md"),

    packages=['taxsbp'],
    install_requires=['binpacking==1.4.1'],

    entry_points = {'console_scripts': ['taxsbp=taxsbp.TaxSBP:main']},

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)