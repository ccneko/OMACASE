# OMACASE
an open-source standalone Python package to analyze optical mapping data, with _de novo_ repeat motif detection and analysis capabilities

## License
- GPL-3.0-or-later

## Written in & for
- Python 3.8+

## Requirements
- Python 3.8+
- Docker (only if running inside a Docker container)
- [wkhtmltopdf](https://github.com/JazzCore/python-pdfkit/wiki/Installing-wkhtmltopdf) for PDF print
    - Make sure X11 server has been started for PDF printing

## Installation
- `python3 -m pip install --user omacase`

## (Alternatively) Runninng inside a Docker container
- `docker build .`
- `docker run {image ID}`
OR
- `docker run ccneko/omacase` (pending registrations)

## Getting Started
Below are some simple instructions to help you start using OMACASE

### Importing the package in Python
- `import omacase`

### Reading BNX / CMAP input in custom Python script
- `your_bnx = opmap.BNX(your_bnx_path)`
- `your_cmap = opmap.CMAP(your_cmap_path)`

## Publication
- Chung CYL, Chan TF. OMACASE: Optical Mapping Quality Control and Repeat Detection Software. (Manuscript in preparation) 

