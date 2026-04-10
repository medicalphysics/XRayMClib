
xraymclib
===================================
xraymclib (Diagnostic X-ray Monte Carlo) is a radiation dose scoring library for diagnostic photon energies written in C++. The main goal for this library is to provide an accurate physics model to describe and model x-ray sources and estimate radiation doses by the Monte-Carlo method.

If you are looking for an application with graphical user interface to perform simulations, OpenXrayMC_ uses this library as Monte Carlo engine and also allows for import CT images and phantoms as scoring volumes. 

.. _OpenXrayMC: https://github.com/medicalphysics/Openxraymc/releases

Example of usage
-----------------
A brief example on how to use xraymclib to simulate a pencilbeam of 60 keV onto a cylindar of aluminum is shown below. 

.. literalinclude:: ../../examples/pencilbeam/pencilbeam.cpp
   :language: c++
   :linenos:

The example can be built by CMake with the following CMakeLists.txt file.

.. literalinclude:: ../../examples/pencilbeam/CMakeLists.txt
   :language: cmake
   :linenos:

See complete examples at:

* Example 1 https://github.com/medicalphysics/xraymcExample1

* Example 2 https://github.com/medicalphysics/xraymcExample2

* Example 3 https://github.com/medicalphysics/xraymcExample3



C++ API
-------------
.. toctree::
   :maxdepth: 1   
   
   api/api_root


Other
-------------
.. toctree::
   :maxdepth: 1   
   
   about
   license
..
   physicsmodel

   

Navigation
==================

* :ref:`genindex`
* :ref:`search`
