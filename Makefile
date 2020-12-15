PYTHON=`which python3`
DESTDIR=/
PROJECT=python-scattnlay
VERSION=2.3
BUILDIR=$(CURDIR)/debian/$(PROJECT)
SRCDIR=$(CURDIR)/src
MULTIPREC=100
CXX_NMIE_HEADERS=$(SRCDIR)/nmie.hpp $(SRCDIR)/nmie-impl.cc $(SRCDIR)/nmie-precision.hpp

all:
	@echo "make source - Create source package for Python extension"
	@echo "make ext - Create Python extension in place"
	@echo "make install - Install Python extension on local system"
	@echo "make rpm - Generate a rpm package for Python extension"
	@echo "make deb - Generate a deb package for Python extension"
	@echo "make standalone - Create standalone programs (scattnlay and fieldnlay)"
	@echo "make clean - Delete temporal files"
#	make standalone

source:
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../

ext: $(SRCDIR)/nmie.cc $(SRCDIR)/nmie-pybind11.cc $(SRCDIR)/pb11_wrapper.cc
	$(PYTHON) setup.py build_ext --inplace

install:
	$(PYTHON) setup.py install --root $(DESTDIR) $(COMPILE)

rpm:
	#$(PYTHON) setup.py bdist_rpm --post-install=rpm/postinstall --pre-uninstall=rpm/preuninstall
	$(PYTHON) setup.py bdist_rpm --dist-dir=../

deb: source
	# build the source package in the parent directory
	# then rename it to project_version.orig.tar.gz
	rename -f 's/$(PROJECT)-(.*)\.tar\.gz/$(PROJECT)_$$1\.orig\.tar\.gz/' ../*
	# build the package
	dpkg-buildpackage -i -I -rfakeroot

standalone: scattnlay-dp fieldnlay-dp scattnlay-mp fieldnlay-mp

# standalone programs with DP
scattnlay-dp: $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o scattnlay-dp $(CXXFLAGS) $(LDFLAGS)

fieldnlay-dp: $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o fieldnlay-dp $(CXXFLAGS) $(LDFLAGS)

# standalone programs with MP
scattnlay-mp: $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/farfield.cc $(SRCDIR)/nmie.cc  -lm -o scattnlay-mp $(CXXFLAGS) $(LDFLAGS)

fieldnlay-mp: $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc $(CXX_NMIE_HEADERS)
	$(CXX) -DNDEBUG -DMULTI_PRECISION=$(MULTIPREC) -O2 -Wall -std=c++11 $(SRCDIR)/nearfield.cc $(SRCDIR)/nmie.cc  -lm -o fieldnlay-mp $(CXXFLAGS) $(LDFLAGS)

clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
