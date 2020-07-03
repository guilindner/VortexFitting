Installation
============

The VortexFitting code works only on  **Python 3.6** and superior.
The following packages are used:

* NumPy
* SciPy
* Matplotlib
* netCDF4

Using PIP
------------

We recommend creating a new virtual environment for your projects

.. code-block:: bash

   $ python3 -m venv NAME
   
   $ source NAME/bin/activate
   
   $ pip install vortexfitting
   
Simply run "pip install vortexfitting"


Documentation
-------------

In order to generate the documentation you need the **sphinx** package.

To perform changes on the documentation, you should modify directly the *.rst*
files, located in the *docs* folder. After modifying the desired files, run the
script *update_docs.py* (for Linux and Windows):

.. code-block:: bash

   $ python update_docs.py

After this, the html files will be generated in the same folder.
