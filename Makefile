PYTHON=`which python3`
DESTDIR=/
PROJECT=python-scattnlay
VERSION=2.3
BUILDIR=$(CURDIR)/debian/$(PROJECT)
SRCDIR=$(CURDIR)/src
MULTIPREC=100
CXX_NMIE_HEADERS=$(SRCDIR)/nmie.hpp $(SRCDIR)/nmie-impl.hpp $(SRCDIR)/nmie-precision.hpp

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

ext: $(SRCDIR)/nmie.cc $(SRCDIR)/nmie-pybind11.cc $(SRCDIR)/pb11_wrapper.cc $(CXX_NMIE_HEADERS)
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

wasm: $(SRCDIR)/nmie-js-wrapper.cc $(CXX_NMIE_HEADERS)
#    emcc -lm -Wall -O2 -std=c++11 --bind -s ENVIRONMENT="web" -s MODULARIZE=1 -s SINGLE_FILE=1 -s WASM=1 -o nmie.js $(SRCDIR)/nmie-js-wrapper.cc
#	emcc --bind -lm -Wall -O2 -std=c++11 -s WASM=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['addOnPostRun']" -o nmiejs.js $(SRCDIR)/nmie-js-wrapper.cc

# 	emcc --bind -lm -Wall -O2 -std=c++11 -s MODULARIZE=1 -s WASM=1 -o nmie.js $(SRCDIR)/nmie-js-wrapper.cc
	emcc --bind -lm -Wall -Oz -std=c++11 -s MODULARIZE=1 -s ASSERTIONS=1 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXPORT_NAME="nmiejs" -s ENVIRONMENT="web" -o nmiejs.js $(SRCDIR)/nmie-js-wrapper.cc
#	emcc --bind -lm -Wall -O3 -std=c++11 -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='["cwrap"]' -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s 'EXPORT_NAME="nmiejs"' -o ./nmiejs.js $(SRCDIR)/nmie-js-wrapper.cc

# 	emcc -g --bind -lm -Wall -std=c++11 -s WASM=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['addOnPostRun']" -o nmie.js $(SRCDIR)/nmie-js-wrapper.cc
# 	emcc --bind -lm -Wall -O2 -std=c++11 -s MODULARIZE=1 -s EXPORT_ES6=1 -s WASM=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['addOnPostRun']" -o nmie.js $(SRCDIR)/nmie-js-wrapper.cc

#     "build:codec": "emcc -O3 -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"cwrap\"]' -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s 'EXPORT_NAME=\"fibonacci\"' -o ./fibonacci.js fibonacci.c",

# emcc -O2 \
#     oniguruma/src/.libs/libonig.so \
#     src/onigasm.cc\
#     -Isrc -Ioniguruma/src \
#     -s ENVIRONMENT=shell \
#     -s NO_EXIT_RUNTIME=1 \
#     -s NO_FILESYSTEM=1 \
#     -s TOTAL_MEMORY=157286400 \
#     -s ALLOW_MEMORY_GROWTH=1 \
#     -s DEMANGLE_SUPPORT=1 \
#     -s EXTRA_EXPORTED_RUNTIME_METHODS="['ccall']" \
#     -s MODULARIZE=1 \
#     -s EXPORT_NAME="'Onigasm'" \
#     -o lib/onigasm.js


	@cp -f nmiejs.js vue-cli3-webapp/src/
	@cp -f nmiejs.wasm vue-cli3-webapp/public/

clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
