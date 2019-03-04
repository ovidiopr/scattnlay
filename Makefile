PYTHON=`which python`
CYTHON=`which cython`
DESTDIR=/
PROJECT=python-scattnlay
VERSION=2.2
BUILDIR=$(CURDIR)/debian/$(PROJECT)
SRCDIR=$(CURDIR)/src
MULTIPREC=100
CXX_NMIE_HEADERS=$(SRCDIR)/nmie.hpp $(SRCDIR)/nmie-impl.hpp $(SRCDIR)/nmie-precision.hpp

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

standalone: scattnlay fieldnlay scattnlay-mp fieldnlay-mp

# standalone programs with DP
scattnlay: $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o scattnlay $(CXXFLAGS) $(LDFLAGS)
fieldnlay: $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o fieldnlay $(CXXFLAGS) $(LDFLAGS)
# standalone programs with MP
scattnlay-mp: $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o scattnlay-mp $(CXXFLAGS) $(LDFLAGS)
fieldnlay-mp: $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o fieldnlay-mp $(CXXFLAGS) $(LDFLAGS)


lib:
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` $(SRCDIR)/nmie.cc $(SRCDIR)/nmie-pybind11.cc  -lm -o example`python3-config --extension-suffix`
clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
