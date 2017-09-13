Tutorial
========

Requirements
------------
You need the following **Python 3** libraries:

* NumPy
* SciPy
* Matplotlib
* netCDF4

Data input
----------
The data used must be in the netCDF4 format. The axis of the velocity
components are ordered like z, y, x, where z can also be the number of
samples in a 2D field.

The variable names should be changed in the **class.py** file

Using the code
--------------

The code comes with a default test case for the HIT, located in *../data/test_data.nc*.
The default configuration os **class.py** is set to handle the HIT case.

If no input file with the *-i* argument has been specified, the *test_data.nc* will be used.

.. code-block:: bash
   
   $ python3 vortexfitting.py

If we want to define a threshold for the swirling strength, we can specify with
the *-s* argument, like this:

.. code-block:: bash

   $ python3 vortexfitting.py -t 0.5

The differencing scheme can be changed with the *-s* argument:

.. code-block:: bash

   $ python3 vortexfitting.py -s 4

.. note:: Available schemes:
          
          * 2 - Second Order (default)
          * 22 - Least-square filter
          * 4 - Fourth order

We can as well change the detection method with the *-d* argument

.. code-block:: bash

   $ python3 vortexfitting.py -d Q

.. note:: Available methods:
          
          * Q - Q criterion
          * swirling - Swirling Strenght (default)
          * delta - Delta criterion
