Examples
========

The example data files can be downloaded from the `github project page <https://github.com/guilindner/VortexFitting/tree/master/data>`_ 

DNS case
--------
This velocity field is from a DNS simulation of an Homogeneous Isotropic Turbulence (HIT) case.

Using default input file *example_dataHIT.nc*, default scheme *least square filter* and
default detection *swirling strength* we test here two cases:

* Case 1: threshold on swirling strength = 0.0
* Case 2: threshold on swirling strength = 0.1
* Case 3: threshold on swirling strength = 0.2
* Case 4: threshold on swirling strength = 0.4

If you set *self.normalization_flag = True* on the classes.py module, the value of swirling
strength will be normalized. The values of the threshold must be set to avoid
unnecessary noise.

Case 1
``````
.. code-block:: bash
   
   $ vortexfitting -t 0.0 -i data/example_dataHIT.nc

.. image:: _images/HIT_00.svg
   :width: 80%
   :alt: detected vortices for example_dataHIT.nc, with a threshold of 0.0
   :align: center

361 vortices detected with 84 accepted.

Case 2
``````
.. code-block:: bash
   
   $ vortexfitting -t 0.1 -i data/example_dataHIT.nc

.. image:: _images/HIT_01.svg
   :width: 80%
   :alt: detected vortices for example_dataHIT.nc, with a threshold of 0.1
   :align: center

162 vortices detected with 72 accepted.

Case 3
``````
.. code-block:: bash
   
   $ vortexfitting -t 0.2 -i data/example_dataHIT.nc

.. image:: _images/HIT_02.svg
   :width: 80%
   :alt: detected vortices for example_dataHIT.nc, with a threshold of 0.2
   :align: center

58 vortices detected with 40 accepted.

Case 4
``````

.. code-block:: bash
   
   $ vortexfitting -t 0.4 -i data/example_dataHIT.nc

.. image:: _images/HIT_04.svg
   :width: 80%
   :alt: detected vortices for example_dataHIT.nc, with a threshold of 0.4
   :align: center

9 vortices detected with 8 accepted.

Below two vortices are displayed, where in the left we have the normal field
and to the right we have the advection velocity subtracted.

.. image:: _images/DNSvortex0_1.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model
.. image:: _images/DNSvortex0_2.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model

.. image:: _images/DNSvortex1_1.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model
.. image:: _images/DNSvortex1_2.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model

+----+---------+--------+--------+
|Case|Threshold|Detected|Accepted|
+====+=========+========+========+
|1   |0.0      |361     |84      |
+----+---------+--------+--------+
|2   |0.1      |162     |72      |
+----+---------+--------+--------+
|3   |0.2      |58      |40      |
+----+---------+--------+--------+
|4   |0.4      |9       |8       |
+----+---------+--------+--------+

PIV case - NetCDF file
----------------------

For PIV data we need to update the format, to match NetCDF file.

It is done with the *-ft piv_netcdf* (*file type*) argument.

Here, since we have an advection velocity, we have to set *self.normalization_flag = True*
and *self.normalization_direction = 'y'*.
This is done directly in the **classes.py** module.

The *-rmax* argument leaves the software calculate the initial radius.

.. code-block:: bash
   
   $ vortexfitting -i data/example_dataPIV.nc -ft piv_netcdf -t 1.5 -rmax 0

.. image:: _images/piv_15.svg
   :width: 90 %

203 vortices detected with 24 accepted.

Below two vortices are displayed, where in the left we have the normal field
and to the right we have the advection velocity subtracted.

.. image:: _images/PIVvortex0_1.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model
.. image:: _images/PIVvortex0_2.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model

.. image:: _images/PIVvortex1_1.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model
.. image:: _images/PIVvortex1_2.png
   :width: 45 %
   :alt: detected vortex and its Lamb-Oseen model

PIV case - Tecplot file
-----------------------

For PIV data with Tecplot, we need to update the format, to match Tecplot file.

It is done with the *-ft piv_tecplot* (*file type*) argument.

.. code-block:: bash
   
   $ vortexfitting.py data/adim_vel_{:06d}.dat -first 10 -t 5 -b 10 -ft piv_tecplot 

.. image:: _images/PIV_accepted_10.svg
   :width: 40 %
   :align: center
   :alt: detected vortices for an experimental case

An average field can be subtracted, using *-mf* argument (*mean file*)

If you want to analyze a set of images, use arguments *-first*, *-last* and *-step*.

(please modify data input to format the image number: *dim_vel_{:06d}.dat* with *-first 10* is formatted as *dim_vel_000010.dat*).

.. code-block:: bash
   
   $ vortexfitting.py -i data/dim_vel_{:06d}.dat -mf data/mean.dat -t 50 -first 10 -ft piv_tecplot

As presented in the image, one main vortex is found in the velocity field provided.

Numerical case - OpenFOAM file
------------------------------

A columnar Lamb-Oseen vortex is generated on OpenFOAM. By default, data are extracted in a text file, with a *.raw* extension.

Here, a z-plane is extracted, with a 100x100 mesh. :math:`(u;v;w)` data are exported.

The spatial mesh for this simulation is quite small, so the default initial radius (*rmax = 10*) is too large.

Specify a smaller value (close to the spatial mesh); *-rmax 0* gets an initial radius of :math:`r_{max} =2\sqrt{dx^2+dy^2}`,

with :math:`dx` and :math:`dy` the spatial resolution.

.. code-block:: bash
   
   $ vortexfitting.py -i data/example_Ub_planeZ_0.01.raw -ft openfoam -rmax 0.0

.. image:: _images/openfoam_quiverplot.png
   :width: 45 %
   :alt: detected vortex for an openfoam simulation of an isolated columnar vortex
.. image:: _images/openfoam.png
   :width: 45 %
   :alt: detected vortex for an openfoam simulation of an isolated columnar vortex
