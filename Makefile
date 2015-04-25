PYTHON=`which python`
DESTDIR=/
BUILDIR=$(CURDIR)/debian/python-scattnlay
PROJECT=python-scattnlay
VERSION=0.3.1

all:
	@echo "make source - Create source package for Python extension"
	@echo "make cython - Convert Cython code to c++"
	@echo "make python_ext - Create Python extension using C++ code"
	@echo "make cython_ext - Create Python extension using Cython code"
	@echo "make install - Install Python extension on local system"
	@echo "make buildrpm - Generate a rpm package for Python extension"
	@echo "make builddeb - Generate a deb package for Python extension"
	@echo "make standalone - Create a standalone programs (scattnlay and fieldnlay)"
	@echo "make clean - Delete temporal files"
#	make standalone

source:
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../

cython: scattnlay.pyx
	cython --cplus scattnlay.pyx

python_ext: nmie.cc  py_nmie.cc scattnlay.cpp
	export CFLAGS='-std=c++11' && python setup.py build_ext --inplace

cython_ext: nmie.cc  py_nmie.cc scattnlay.pyx
	export CFLAGS='-std=c++11' && python setup_cython.py build_ext --inplace

install:
	$(PYTHON) setup.py install --root $(DESTDIR) $(COMPILE)

buildrpm:
	#$(PYTHON) setup.py bdist_rpm --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall
	$(PYTHON) setup.py bdist_rpm --dist-dir=../

builddeb:
	# build the source package in the parent directory
	# then rename it to project_version.orig.tar.gz
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../ --prune
	rename -f 's/$(PROJECT)-(.*)\.tar\.gz/$(PROJECT)_$$1\.orig\.tar\.gz/' ../*
	# build the package
	export CFLAGS='-std=c++11' && dpkg-buildpackage -i -I -rfakeroot

standalone: standalone.cc nmie.cc
	export CFLAGS='-std=c++11' && c++ -DNDEBUG -O2 -Wall -std=c++11 farfield.cc nmie.cc  -lm -o scattnlay
	mv scattnlay ../
	export CFLAGS='-std=c++11' && c++ -DNDEBUG -O2 -Wall -std=c++11 nearfield.cc nmie.cc  -lm -o fieldnlay
	mv fieldnlay ../

clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
