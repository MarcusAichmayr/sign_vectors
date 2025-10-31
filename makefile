.PHONY: install test

install:
	sage -pip install --upgrade .

test:
	sage -t sign_vectors/

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
