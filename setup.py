from setuptools import setup


# Get information from separate files (README)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


setup(
    name="sign_vectors",
    version="1.2",
    description="SageMath package for sign vectors, oriented matroids and chirotopes",
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/MarcusAichmayr/sign_vectors",
    author="Marcus S. Aichmayr",
    author_email="aichmayr@mathematik.uni-kassel.de",
    license="GPL-3.0-or-later",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords=[
        "sign vectors",
        "oriented matroids",
        "chirotopes",
    ],
    packages=["sign_vectors",],
    extras_require={
        "passagemath": [
            "passagemath-symbolics",
            "passagemath-flint",
            "passagemath-graphs",
            "passagemath-repl",
        ],
    },
)
