Examples
========

DNS case
--------
DNS simulation of an Homogeneous Isotropic Turbulence (HIT) case.

Using default input file *test_caseHIT.nc*, default scheme *second order* and
default detection *swirling strength* we test here two cases:

* Case 1: Threshold on swirling strength = 0.0
* Case 2: Threshold on swirling strength = 0.1
* Case 3: Threshold on swirling strength = 0.2
* Case 4: Threshold on swirling strength = 0.4

if you set **self.norm = True** on the classes.py module, the value of swirling
strength will be normalized. The values of the threshold must be set to avoid
unnecessary noise.

Case 1
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.0

.. image:: _images/HIT_00.svg
   :width: 90%

361 vortices detected with 88 accepted.

Case 2
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.1

.. image:: _images/HIT_01.svg
   :width: 90%

162 vortices detected with 75 accepted.

Case 3
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.2

.. image:: _images/HIT_02.svg
   :width: 90%

58 vortices detected with 42 accepted.

Case 4
``````

.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.4

.. image:: _images/HIT_04.svg
   :width: 90%

9 vortices detected with 8 accepted.

Below two vortices are displayed, where in the left we have the normal field
and to the right we have the convection velocity subtracted.

.. image:: _images/DNSvortex0_1.png
   :width: 45 %
.. image:: _images/DNSvortex0_2.png
   :width: 45 %

.. image:: _images/DNSvortex1_1.png
   :width: 45 %
.. image:: _images/DNSvortex1_2.png
   :width: 45 %

+----+---------+--------+--------+
|Case|Threshold|Detected|Accepted|
+====+=========+========+========+
|1   |0.0      |361     |88      |
+----+---------+--------+--------+
|2   |0.1      |162     |75      |
+----+---------+--------+--------+
|3   |0.2      |58      |42      |
+----+---------+--------+--------+
|4   |0.4      |9       |8       |
+----+---------+--------+--------+

PIV case - NetCDF file
----------------------

For PIV data we need to change the *classes.py* to match the NetCDF file:

.. code-block:: python

   self.u = np.array(grp1.variables['velocity_n'][time,:,:])
   self.v = np.array(grp1.variables['velocity_s'][time,:,:])
   self.w = np.array(grp1.variables['velocity_z'][time,:,:])
   self.dx = np.array(grp1.variables['grid_n'])
   self.dy = np.array(grp1.variables['grid_z'])
   self.dy = self.dy - self.dy[0] #it does not start at 0
   self.u = self.u - np.mean(self.u,1)[:,None]
   self.v = self.v - np.mean(self.v,1)[:,None]
   self.w = self.w - np.mean(self.w,1)[:,None]
   self.norm = True
   self.normdir = 0
   self.samples = self.u.shape[1]

Here since we have a convection velocity, we have to set the *self.norm = True*
and the *self.normdir = 0* (for the y direction)

.. code-block:: bash
   
   $ python3 vortexfitting.py -i ../data/test_dataPIV -t 1.5

.. image:: _images/piv_15.svg
   :width: 90 %

203 vortices detected with 29 accepted.

Below two vortices are displayed, where in the left we have the normal field
and to the right we have the convection velocity subtracted.

.. image:: _images/PIVvortex0_1.png
   :width: 45 %
.. image:: _images/PIVvortex0_2.png
   :width: 45 %

.. image:: _images/PIVvortex1_1.png
   :width: 45 %
.. image:: _images/PIVvortex1_2.png
   :width: 45 %

PIV case - Tecplot file
-----------------------

For PIV data with Tecplot, we need to update the format, to match Tecplot file.

.. code-block:: bash
   
   $ python3 vortexfitting.py -i ../data/adim_vel_{:06d}.dat -first 10 -last 10 -t 5 -b 10 -ft piv_tecplot 

.. image:: _images/PIV_accepted_10.svg
   :width: 40 %

An average field can be subtracted, using *-mf* argument (*mean file*)

If you want to analyze a set of images, use arguments *-first*, *-last* and *-step*.

(please modify data input to format the image number: *dim_vel_{:06d}.dat* with *-first 10* is formatted as *dim_vel_000010.dat*).

.. code-block:: bash
   
   $ python3 vortexfitting.py -i ../data/dim_vel_{:06d}.dat -mf ../data/mean.dat -t 50 -first 10 -last 10 -ft piv_tecplot

