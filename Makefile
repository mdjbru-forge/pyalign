### * Variables

PYTHON_MODULE=pyalign
PYTHON_MODULE_EGG=$(PYTHON_MODULE).egg-info

COVERED_PACKAGES=$(PYTHON_MODULE)
SPHINX_DOC_FOLDER=docs/

### * Help

help:
	@echo "Makefile for the $(PYTHON_MODULE) Python module                "
	@echo "                                                               "
	@echo "Type \"make <target>\" where <target> is one of the following: "
	@echo "                                                               "
	@echo "  install     Install the module                               "
	@echo "  uninstall   Uninstall the module                             "
	@echo "  test        Run the tests with coverage output               "
	@echo "  clean       Clean everything (coverage, tests, pyc files)    "


### * Main targets

### ** test
tests: test
test:
	nosetests tests/ --with-coverage --cover-package=$(COVERED_PACKAGES) --cover-html \
          --with-html --html-file=tests/nosetests.html
	@echo -e "\nThe coverage results are accessible from cover/index.html"
	@echo "The html version of the test results are accessible from tests/nosetests.html"

### ** clean
clean:
	rm -f .coverage
	rm -fr cover
	rm -f tests/nosetests.html
	@# http://superuser.com/questions/112078/delete-matching-files-in-all-subdirectories
	find . -name \*.pyc -type f -delete

### ** install
install:
	rm -fr $(PYTHON_MODULE_EGG)
	pip install --upgrade --user .

### ** uninstall
uninstall:
	pip uninstall -y $(PYTHON_MODULE)
