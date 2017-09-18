Tutorial
========

Requirements
------------
You need the following **Python 3** libraries:

* NumPy
* SciPy
* Matplotlib
* netCDF4

In order to generate the documentation you need the **sphinx** package.


Using the code
--------------

Data input
``````````

The data used must be in the netCDF4 format. The axis of the velocity
components are ordered like z, y, x, where z can also be the number of
samples in a 2D field.

The variable names should be changed in the **class.py** file

In this example we have the DNS HIT format used in this tutorial:

.. code-block:: python

    self.u = np.array(grp1.variables['velocity_x'][time,:,:])
    self.v = np.array(grp1.variables['velocity_y'][time,:,:])
    self.w = np.array(grp1.variables['velocity_z'][time,:,:])
    self.samples = self.u.shape[1]
    self.dx = np.linspace(0,self.samples,self.samples)
    self.dy = np.linspace(0,self.samples,self.samples)
    self.norm = False
    self.normdir = 1

Here the names of the variables *velocity_x*, *velocity_y* and *velocity_z* are
defined and an equally spaced mesh is created using *np.linspace*. 

.. note:: No normalization is required for this case due to its homogeneity.


Parameters
``````````

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

Data output
```````````

The results will be written to the *../results/* folder with the following files:

* accepted.svg: The location and size of the accepted vortices
* linked.svg: same as *accepted.svg* but can be open on the web browser with
  clickable vortices
* vortex#_1.png: Comparison of the velocity field of the vortex and the model
* vortex#_2.png: Comparison of the velocity field of the vortex and the model,
  subtracting the convection velocity
* vortices.dat: parameters of all the vortices


Generating a custom Vortex
--------------------------

It's possible to generate a custom vortex using the **generateNetCDF.py** module.
It will create a netCDF4 file with the same characteristics as the DNS HIT file.

.. code-block:: bash

   $ python3 generateNetCDF.py

This command will create a file *generatedField.nc* at the data folder.

You can tune the characteristics and position of the vortex by changing the 
following values directly on *generatedField.nc*:

* coreR;
* gamma;
* fxCenter;
* u_conv;
* v_conv.

The size of the domain can also be changed on the *ndim* variable.

Converting NC to ASCII
----------------------

If for any reason you need to convert the NC file to a text format (ASCII), the
module **convertToASCII.py** can do the job. It will open the *infile* and save
all z planes (or time) into separeted files.

.. code-block:: bash

   $ python3 convertToASCII.py -i input.nc -o output

Depending on the file you need to change the variable names like *velocity_x*
and such for the corresponding variable.
