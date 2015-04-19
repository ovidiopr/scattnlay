![output example](/doc/SiAgSi.png)
**Output example:** Field distribution inside layered Si\Ag\Si sphere
calculated with Scattnlay.

How to use scattnlay
====================

**Table of contents:**
- [Compile code](#compile-code)
- [Use](#use)
- [Papers](#papers)
- [Acknowledgment](#acknowledgment)
- [License](#license)

Compile Code:
-------------

Compilation options

 - **make source** - Create source package (python library)
 - **make install** - Install on local system (python library)
 - **make buildrpm** - Generate a rpm package (python library)
 - **make builddeb** - Generate a deb package (python library)
 - **make standalone** - Create a standalone program
 - **make clean** - Delete temporal files

Use:
----

1. Python library
  * Use scattnlay directly
  
  ```python
from scattnlay import scattnlay
...
x = ...
m = ...
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
...
  ```
  
  * Execute some of the test scripts (located in the folder 'tests/python')
          Example:
		  
  ```bash
./test01.py
  ```
  
2. Standalone program
  * Execute scattnlay directly
          Usage:
		  
  ```bash
scattnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-c comment]
  ```
  * Execute some of the test scripts (located in the folder 'tests/shell')
      Example:

  ```bash
./test01.sh > test01.csv
  ```
3. C++ library

```C++
    try {
      MultiLayerMie multi_layer_mie;
      multi_layer_mie.SetLayersSize(x);
      multi_layer_mie.SetLayersIndex(m);

      multi_layer_mie.RunMieCalculation();

      *Qsca = multi_layer_mie.GetQsca();
      *Qabs = multi_layer_mie.GetQabs();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  multi_layer_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }
```

Papers
------

1. "Scattering of electromagnetic radiation by a multilayered sphere"
   O. Pena and U. Pal,  Computer Physics Communications, vol. 180,
   Nov. 2009, pp. 2348-2354. http://dx.doi.org/10.1016/j.cpc.2009.07.010

2. "Reduction of scattering using thin all-dielectric shells designed by stochastic optimizer"
   Konstantin Ladutenko, Ovidio Peña-Rodríguez, Irina Melchakova, Ilya
   Yagupov, and Pavel Belov  J. Appl. Phys., vol. 116, pp. 184508,
   2014 http://dx.doi.org/10.1063/1.4900529 

Acknowledgment
--------------

We expect that all publications describing work using this software,
or all commercial products using it, cite the following reference:
> O. Pena and U. Pal, "Scattering of electromagnetic radiation
> by a multilayered sphere," Computer Physics Communications,
> vol. 180, Nov. 2009, pp. 2348-2354.

License
-------

GPL v3+
