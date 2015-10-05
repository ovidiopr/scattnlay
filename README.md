![output example](/doc/OutputExample.png)
**Output example:** Field distribution inside layered Si\Ag\Si sphere
and Poynting vector distribution in Ag sphere with poweflow lines
calculated with Scattnlay (scripts  field-SiAgSi-flow.py and
field-Ag-flow.py from example section as
57c7261705a5776f78420c1f486e929517d5f584 ).

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

 - **make source** - Create source package for Python extension
 - **make cython** - Convert Cython code to C++
 - **make python_ext** - Create Python extension using C++ code
 - **make cython_ext** - Create Python extension using Cython code
 - **make install** - Install Python extension on local system
 - **make buildrpm** - Generate a rpm package for Python extension
 - **make builddeb** - Generate a deb package for Python extension
 - **make standalone** - Create standalone programs (scattnlay and fieldnlay)
 - **make clean** - Delete temporal files

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
