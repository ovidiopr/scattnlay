PYTHON=`which python`
CYTHON=`which cython`
DESTDIR=/
PROJECT=python-scattnlay
VERSION=2.2
BUILDIR=$(CURDIR)/debian/$(PROJECT)
SRCDIR=$(CURDIR)/src
MULTIPREC=100

all:
	@echo "make source - Create source package for Python extension"
	@echo "make cython - Convert Cython code to C++"
	@echo "make python_ext - Create Python extension using C++ code"
	@echo "make cython_ext - Create Python extension using Cython code"
	@echo "make install - Install Python extension on local system"
	@echo "make buildrpm - Generate a rpm package for Python extension"
	@echo "make builddeb - Generate a deb package for Python extension"
	@echo "make standalone - Create standalone programs (scattnlay and fieldnlay)"
	@echo "make clean - Delete temporal files"
#	make standalone

source:
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../

cython: $(SRCDIR)/scattnlay.pyx
    # create c++ code for double precision module
	$(CYTHON) --cplus $(SRCDIR)/scattnlay.pyx -o $(SRCDIR)/scattnlay.cpp
	# create c++ code for MP module
	ln -s $(SRCDIR)/scattnlay.pyx $(SRCDIR)/scattnlay_mp.pyx
	$(CYTHON) --cplus $(SRCDIR)/scattnlay_mp.pyx -o $(SRCDIR)/scattnlay_mp.cpp
	rm $(SRCDIR)/scattnlay_mp.pyx

python_ext: $(SRCDIR)/nmie.cc $(SRCDIR)/py_nmie.cc $(SRCDIR)/scattnlay.cpp $(SRCDIR)/scattnlay_mp.cpp
	$(PYTHON) setup.py build_ext --inplace

cython_ext: cython python_ext

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
	dpkg-buildpackage -i -I -rfakeroot

standalone: $(SRCDIR)/farfield.cc $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc
    # create standalone programs with DP
	c++ -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o ../scattnlay
	c++ -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o ../fieldnlay
	# create standalone programs with MP
	c++ -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o ../scattnlay-mp
	c++ -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o ../fieldnlay-mp

clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
