# Installation

## How to install OMACASE
OMACASE can be installed via one of the following ways:
- From source:
    1. Download the code release from GitHub
    2. Run `python3 -m pip install .` in the code directory; OR
    3. run `docker build .` in the code directory
- From package / container image repository (to be available at manuscript submission):
    - `python3 -m pip install omacase`; OR
    - `docker pull omacase`

## Prerequisites
- [Python](https://www.python.org/) 3.8+
- [Docker](https://www.docker.com/) (only if running inside a Docker container)
- [wkhtmltopdf](https://github.com/JazzCore/python-pdfkit/wiki/Installing-wkhtmltopdf) for PDF print
    - Make sure X11 server has been started for PDF printing
