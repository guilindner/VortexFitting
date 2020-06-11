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

The data used can be of different format.
For netCDF4, the axis of the velocity components are ordered like z, y, x, 
where z can also be the number of samples in a 2D field.

TecPlot format or OpenFoam export are also accepted.

The variable names should be changed in the **class.py** file

In this example we have the DNS HIT format used in this tutorial:

.. code-block:: python

    self.u_velocity_matrix = np.array(datafile_read.variables['velocity_x'][time_step, :, :])
    self.v_velocity_matrix = np.array(datafile_read.variables['velocity_y'][time_step, :, :])
    self.w_velocity_matrix = np.array(datafile_read.variables['velocity_z'][time_step, :, :])
    self.x_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[1], self.u_velocity_matrix.shape[1])
    self.y_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[0], self.u_velocity_matrix.shape[0])
    self.z_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[0], self.u_velocity_matrix.shape[0])
    self.normalization_flag = False
    self.normalization_direction = False

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

An initial radius can be set with *-rmax* argument. 

.. code-block:: bash

   $ python3 vortexfitting.py -rmax 15

An output directory can be specified with *-o* / *--output* argument. 

.. code-block:: bash

   $ python3 vortexfitting.py -o ../results/MY_DIRECTORY

Use arguments *-first*, *-last* and *-step* to analyze a set of images. Default for *-step* is 1.

For example, if you want to compute from image #10 to #20, each 2 images, enter:

.. code-block:: bash

   $ python3 vortexfitting.py -first 10 -last 20 -step 2


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


Generating a custom Vortex
--------------------------

It's possible to generate a custom vortex using the **generateNetCDF.py** module.
It will create a netCDF4 file with the same characteristics as the DNS HIT file.

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


Converting NC to ASCII
----------------------

If for any reason you need to convert the NC file to a text format (ASCII), the
module **convertToASCII.py** can do the job. It will open the *infile* and save
all z planes (or time) into separated files.

.. code-block:: bash

   $ python3 convertToASCII.py -i input.nc -o output

Depending on the file you need to change the variable names like *velocity_x*
and such for the corresponding variable.



Documentation
-------------

To perform changes on the documentation, you should modify directly the *.rst*
files, located in the *docs* folder. After modifying the desired files, run the
script *update_docs.sh*:

.. code-block:: bash

   $ sh update_docs.sh

After this, the html files will be generated in the same folder. All files in
the doc folder should be commited to github.com to appear online.
