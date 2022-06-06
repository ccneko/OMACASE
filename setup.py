#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20211111


"""
OMACASE, a tool to analyze optical mapping data
with _de novo_ repeat motif detection and analysis capabilities
distributed under license GPL-3.0-or-later
By Claire Chung
"""


from setuptools import setup, find_packages

from omacase import __author__, __version__, __email__

with open('README.md', 'r') as f:
    long_description = f.read()

with open('LICENSE.md', 'r') as f:
    license_text = f.read()

with open('requirements.txt', 'r') as f:
    install_reqs = f.readlines()
print(install_reqs)

setup(
    name="omacase",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="a tool to analyze optical mapping data, with _de novo_ repeat motif detection and analysis capabilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license=license_text,
    url="https://github.com/ccneko/OMACASE",
    install_requires=install_reqs,
    packages=find_packages(),
    include_package_data=True,
    package_data={"omacase":["assets/*/*", "config/*", "test_data/*"]},
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    platforms='any',
)