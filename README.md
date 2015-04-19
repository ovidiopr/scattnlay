How to use scattnlay
====================

Compile Code:
-------------

Compilation options

 - make source - Create source package (python library)
 - make install - Install on local system (python library)
 - make buildrpm - Generate a rpm package (python library)
 - make builddeb - Generate a deb package (python library)
 - make standalone - Create a standalone program
 - make clean - Delete temporal files

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
          Usage: scattnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-c comment]

      * Execute some of the test scripts (located in the folder 'tests/shell')
          Example:  ./test01.sh > test01.csv

*******************************************************************************

