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

Case 1
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.0

.. image:: _images/HIT_00.svg
   :width: 360px
   :height: 270px
   :align: center
   :alt: HIT case with 0.0 threshold on swirling strength

361 vortices detected with 141 accepted.

case 2
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.1

.. image:: _images/HIT_01.svg
   :width: 360px
   :height: 270px
   :align: center
   :alt: hit case with 0.1 threshold on swirling strength

162 vortices detected with 108 accepted.

case 3
``````
.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.2

.. image:: _images/HIT_02.svg
   :width: 360px
   :height: 270px
   :align: center
   :alt: hit case with 0.2 threshold on swirling strength

58 vortices detected with 51 accepted.

Case 4
``````

.. code-block:: bash
   
   $ python3 vortexfitting.py -t 0.4

.. image:: _images/HIT_04.svg
   :width: 360px
   :height: 270px
   :align: center
   :alt: HIT case with 0.4 threshold on swirling strength

9 vortices detected with 8 accepted.

+----+---------+--------+--------+
|Case|Threshold|Detected|Accepted|
+====+=========+========+========+
|1   |0.0      |361     |141     |
+----+---------+--------+--------+
|2   |0.1      |162     |108     |
+----+---------+--------+--------+
|3   |0.2      |58      |51      |
+----+---------+--------+--------+
|4   |0.4      |9       |8       |
+----+---------+--------+--------+


PIV case
--------
