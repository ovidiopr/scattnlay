![output example](/doc/OutputExample.png)
**Output example:** Field distribution inside layered Si\Ag\Si sphere
and Poynting vector distribution in Ag sphere with powerflow lines
calculated with Scattnlay (scripts  field-SiAgSi-flow.py and
field-Ag-flow.py from example section as [revision](https://github.com/ovidiopr/scattnlay/commit/57c7261705a5776f78420c1f486e929517d5f584) ).

Discuss:
--------

Try to join our Gitter chat: [![Join the chat at https://gitter.im/scattnlay/Lobby](https://badges.gitter.im/scattnlay/Lobby.svg)](https://gitter.im/scattnlay/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
 
Fill the issue here: [Issues](https://github.com/ovidiopr/scattnlay/issues).

Stable releases
===============

- Version 2.0.1 (Jan 17, 2017). [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.248729.svg)](https://doi.org/10.5281/zenodo.248729)
- Version 2.0.0 (Apr 1, 2016).
- Version 1.0.0 (Nov 22, 2014).

How to use scattnlay
====================

**Table of contents:**
- [Compile code](#compile-code)
- [Binary install](#binary-install)
- [Use](#use)
- [Papers](#papers)
- [Acknowledgment](#acknowledgment)
- [License](#license)

Compile Code:
-------------
To compile the source you will need a C++11 capable compiler. To use
MultiPrecision feature you need to install Boost.Multiprecision
library:

 - **libboost-all-dev (>= 1.58.0)**

To compile the Python extension you need [NumPy](http://www.numpy.org/):

 - **python-numpy (>= 1.0)**
 - **python-all-dev (any version)**
 - **python-numpy-dev (any version)**

And to compile the Debian package you need some tools:

 - **debhelper (>=7.0.0)**
 - **dh-python (any version)**
 - **cdbs (>= 0.4.49)**

Compilation options

 - **make src** - Create source package for Python extension
 - **make ext** - Create Python extension using C++ code
 - **make install** - Install Python extension on local system
 - **make rpm** - Generate a rpm package for Python extension
 - **make deb** - Generate a deb package for Python extension
 - **make standalone** - Create standalone programs (scattnlay and fieldnlay)
 - **make clean** - Delete temporal files

Binary install:
--------------

Binary files for Ubuntu and derivative distributions can be found at
[Launchpad](https://launchpad.net/~ovidio/+archive/ubuntu/scattering/+packages)
To install it you must configure the repository:
``` bash
sudo add-apt-repository ppa:ovidio/scattering
sudo apt-get update
```
and then you simply install the package:
``` bash
sudo apt-get install python-scattnlay
```
You can also install it from PyPi via
```bash
sudo pip install python-scattnlay
```

Use:
----

1. Python library
  * Use scattnlay directly
  
  ```python
from scattnlay import scattnlay, fieldnlay
...
x = ...
m = ...
coords = ...
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
terms, E, H = fieldnlay(x, m, coords)
...
  ```
  
  * Execute some of the test scripts (located in the folder 'tests/python')
          Example:
		  
  ```bash
./test01.py
  ```
  
2. Standalone program
  * Execute scattnlay directly:

  ```bash
scattnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]
  ```

  * Execute fieldnlay directly:

  ```bash
fieldnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] -p xi xf nx yi yf ny zi zf nz [-c comment]
  ```

  * Execute some of the test scripts (located in the folder 'tests/shell'):

  ```bash
./test01.sh > test01.csv
  ```
  
3. C++ library

Scattnlay "Hello world!" example:

```C++
    try {
      nmie::MultiLayerMieApplied<double> multi_layer_mie; 
      multi_layer_mie.AddTargetLayer(core_width, index_Si);
      multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
      multi_layer_mie.AddTargetLayer(outer_width, index_Si);
      multi_layer_mie.SetWavelength(WL);
      multi_layer_mie.RunMieCalculation();
      double Qabs = multi_layer_mie.GetQabs();
      printf("Qabs = %g\n", Qabs);
    } catch( const std::invalid_argument& ia ) {
      // Will catch if  multi_layer_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      return -1;
    }
```

The complete `example-minimal.cc` and a bit more complicated
`example-get-Mie.cc` can be found in example directory along with
`go-cc-examples.sh` script with build commands.

`example-get-Mie.cc` can be compiled using double precision or
multiple precision (just include `-DMULTI_PRECISION=200` to use 200
digits for calculations). 

Related papers
--------------

1. O. Peña and U. Pal, "Scattering of electromagnetic radiation by a
   multilayered sphere," Comput. Phys. Commun. 180, 2348-2354 (2009).
   http://dx.doi.org/10.1016/j.cpc.2009.07.010

2. K. Ladutenko, O. Peña-Rodríguez, I. Melchakova, I. Yagupov and P. Belov,
   "Reduction of scattering using thin all-dielectric shells designed by
   stochastic optimizer," J. Appl. Phys. 116, 184508 (2014).
   http://dx.doi.org/10.1063/1.4900529

3. K. Ladutenko, P. Belov, O. Peña-Rodríguez, A. Mirzaei, A. Miroshnichenko
   and I. Shadrivov, "Superabsorption of light by nanoparticles,"
   Nanoscale 7, 18897-18901 (2015).
   http://dx.doi.org/10.1039/C5NR05468K

4. K. Ladutenko, U. Pal, A. Rivera, and O. Peña-Rodríguez, "Mie calculation
   of electromagnetic near-field for a multilayered sphere,"
   Comp. Phys. Comm. 214, 225-230 (2017).
   http://dx.doi.org/j.cpc.2017.01.017

Acknowledgment
--------------

We expect that all publications describing work using this software,
or all commercial products using it, cite at least one of the following references:
> [1] O. Peña and U. Pal, "Scattering of electromagnetic radiation
>     by a multilayered sphere," Computer Physics Communications,
>     vol. 180, Nov. 2009, pp. 2348-2354.
>
> [2] K. Ladutenko, U. Pal, A. Rivera and O. Peña-Rodríguez, "Mie calculation
>     of electromagnetic near-field for a multilayered sphere,"
>     Computer Physics Communications, vol. 214, May 2017, pp. 225-230.


License
-------

GPL v3+
