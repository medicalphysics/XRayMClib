
DXMClib
===================================
DXMClib (Diagnostic X-ray Monte Carlo) is a radiation dose scoring library for diagnostic photon energies written in C++. The main goal for this library is to provide an accurate enough physics model to describe and model x-ray sources and estimate radiation doses by the Monte-Carlo method.

If you are looking for an application with graphical user interface to perform simulations, OpenDXMC_ uses this library as Monte Carlo engine and also allows for import CT images and phantoms as scoring volumes. 

.. _OpenDXMC: https://github.com/medicalphysics/OpenDXMC/releases

Example of usage
-----------------
A brief example on how to use DXMClib to simulate a pencilbeam of 60 keV onto a cylindar of aluminum is shown below. 

.. literalinclude:: ../../examples/pencilbeam/pencilbeam.cpp
   :language: c++
   :linenos:

The example can be built by CMake with the following CMakeLists.txt file.

.. literalinclude:: ../../examples/pencilbeam/CMakeLists.txt
   :language: cmake
   :linenos:


Check also out
--------------
.. toctree::
   :maxdepth: 2
   
   about
   physicsmodel
   api/api_root
   license
   
Docs
====
.. doxygenclass:: CTDualSource
   :project: DXMClib
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
