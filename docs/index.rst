Welcome to Epidemics on Networks's documentation!
=================================================


**EoN** (Epidemics on Networks)  is a Python module that provides tools to study 
the spread of SIS and SIR diseases in networks.  

.. toctree::
   :maxdepth: 2

   Getting Started <GettingStarted>
   Examples 
   EoN 

Highlights 
----------

**EoN** is based on the book  

  `Mathematics of Epidemics on Networks: from Exact to Approximate Models`_ 

  **EoN** is built on top of NetworkX_.  Its repository_ is on github.   
  EoN's tools fall into two broad categories:

- **Stochastic simulation of SIS and SIR disease**

  - Event-based simulation 
  
    - much faster than traditional Gillespie simulation
    - allows weighted graphs 
    - allows non-Markovian dynamics
  - Gillespie algorithms for Markovian dynamics on unweighted graphs
  - discrete-time (synchronous update) models
- **Numerical solvers for ODE models**

  - pair approximation models
  - effective degree models
  - edge-based compartmental models

There is also support for producing figures and animations and other detailed
investigation of simulations.

.. _repository: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks
.. _Mathematics of epidemics on networks\: from exact to approximate models: http://www.springer.com/us/book/9783319508047
.. _NetworkX: https://networkx.github.io


   
