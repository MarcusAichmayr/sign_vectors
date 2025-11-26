# Sign vectors

## Description

SageMath package for sign vectors, oriented matroids and chirotopes

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 3 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

Sage 10.0 or later is recommended.

## Installation

### Install from GitHub (recommended)

To install the latest development version on a system where Sage
was built from source or installed from official packages, run:

    sage -pip install git+https://github.com/MarcusAichmayr/sign_vectors.git

### Install from PyPI (recommended)

To install the package for Python directly, run:

    pip install sign-vectors

### Local install from source

Download the source from the git repository:

    git clone https://github.com/MarcusAichmayr/sign_vectors.git

Change to the root directory of the repository and run:

    make install

### Local install from source (no Sage installation required)

Download the source from the git repository:

    git clone https://github.com/MarcusAichmayr/sign_vectors.git

Change to the root directory of the repository and run:

    python3 -m venv venv
    . venv/bin/activate
    pip install -v -e ".[passagemath]"

## Documentation

The documentation of this package is available on GitHub:

https://marcusaichmayr.github.io/sign_vectors/

To generate it, run

    make doc

or

    make doc-pdf

at the root directory of the repository.

## Testing

To run the test suite, install the package and run the command

    make test

at the root directory of the repository.
