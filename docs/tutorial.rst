Tutorial
========

Using the code
--------------

Data input
``````````

The data used can be of different format.
For NetCDF, the axis of the velocity components are ordered like z, y, x, 
where z can also be the number of samples in a 2D field.

Tecplot format or OpenFoam export are also accepted.

If you want to implement a custom format, change the variable names, their order ...
it should be directly changed in the **class.py** file


Parameters
``````````

The code comes with default test cases, located in *../data/*.

If no input file with the *-i* (*--input*) argument has been specified, the *example_dataHIT.nc* will be used.

.. code-block:: bash
   
   $ vortexfitting

If we want to define a threshold for the swirling strength, we can specify with
the *-t* (*--threshold*) argument, like this:

.. code-block:: bash

   $ python3 vortexfitting.py -t 0.5

The differencing scheme can be changed with the *-s* (*--scheme*) argument:

.. code-block:: bash

   $ python3 vortexfitting.py -s 4

.. note:: Available schemes:
          
          * 2 - Second Order (default)
          * 22 - Least-square filter
          * 4 - Fourth order

We can as well change the detection method with the *-d* (*--detect*) argument

.. code-block:: bash

   $ python3 vortexfitting.py -d Q

.. note:: Available methods:
          
          * Q - Q criterion
          * swirling - Swirling Strenght (default)
          * delta - Delta criterion

An initial guessing radius can be set with *-rmax* argument. 

.. code-block:: bash

   $ python3 vortexfitting.py -rmax 15

An output directory can be specified with the *-o* (*--output*) argument. 

.. code-block:: bash

   $ python3 vortexfitting.py -o ../results/MY_DIRECTORY

Use arguments *-first*, *-last* and *-step* to analyze a set of images. Default for *-step* is 1.

For example, if you want to compute from image #10 to #20, each 2 images, enter:

.. code-block:: bash

   $ python3 vortexfitting.py -first 10 -last 20 -step 2


By default, the correlation threshold to detect a vortex is 0.75. This value may be changed with the
*-ct* (*--corrthreshold*) argument.

.. code-block:: bash

   $ python3 vortexfitting.py -ct 0.85

To avoid vortices overlapping, the box size parameter *-b* (*--boxsize*) can be used. 
It takes an integer distance in mesh units, between two vortex centers.

.. code-block:: bash

   $ python3 vortexfitting.py -b 10


The plot method is chosen with the *-p* (*--plot*) argument

.. note:: Available schemes:
          
          * fit - detection and fitting, saves images (default)
          * detect - Locate the potenial vortices (without fitting)
          * fields - display the velocity fields and vorticity

.. code-block:: bash

   $ python3 vortexfitting.py -p fields



Data output
```````````

The results will be written to the *../results/* folder with the following files:

* accepted.svg: The location and size of the accepted vortices
* linked.svg: same as *accepted.svg* but can be open on the web browser with
  clickable vortices
* vortex#_initial_vfield.png: Comparison of the velocity field of the vortex and the model
* vortex#_advection_field_subtracted.png: Comparison of the velocity field of the vortex and the model,
  subtracting the advection velocity
* vortices.dat: parameters of all the detected vortices

If you want to update the output format of *vortices.dat*, it should be done in the **output.py** file.

The plot results are handled in the **fitting.py** module.

Generating a custom Vortex
--------------------------

It's possible to generate a custom vortex using the **generateNetCDF.py** module.
It will create a NetCDF file with the same characteristics as the DNS HIT file.

.. code-block:: bash

   $ python3 generateNetCDF.py

This command will create a file *generatedField.nc* at the data folder.

You can tune the characteristics and position of the vortex by changing the 
following values directly on *generatedField.nc*:

* core_radius;
* gamma;
* x_center;
* y_center;
* u_advection;
* v_advection.

The size of the domain can also be changed on the *ndim* variable.

You can use the *output* option (*-o*) to specify the name of the created file, 
and *ndim* (*-ndim*) option to change the domain size.
For example:

.. code-block:: bash

   $ python3 generateNetCDF.py -o ./data/testGenerate.nc -ndim 300


will produce a 300x300 mesh, in a file named *testGenerate.nc*.


Converting NetCDF to ASCII (and vice-versa)
-------------------------------------------

If for any reason you need to convert a netCDF file to a text format (ASCII), the
module **convertToASCII.py** can do the job. It will open the *infile* and save
all z planes (or time) into separated files.

.. code-block:: bash

   $ python3 convertToASCII.py -i input.nc -o output.dat

Depending on the file you need to change the variable names like *velocity_x*
and such for the corresponding variable.

The module **convertToNC.py** can convert an ASCII file to a NetCDF4 format.
You can specify the spatial dimensions (*nx*, *ny* respectively for x and y directions),
with the options *-nx* or *-ny*

.. code-block:: bash

   $ python3 convertToNC.py -i input.dat -o output.nc
