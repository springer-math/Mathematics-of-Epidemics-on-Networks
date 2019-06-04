Welcome to Epidemics on Networks's documentation!
=================================================


**EoN** (Epidemics on Networks) is a Python module that provides tools to study 
the spread of SIS and SIR diseases in networks. 

**Support EoN**:  
  - The best way to support EoN `is to let me know you're using it`_.
  - This will help my case when applying for grants & promotions and help me justify the time I spend on it. 

**MIT License**:
See :download:`license.txt<../license.txt>` for 
full details.


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
    
  - Gillespie algorithms for Markovian dynamics

    - Through some careful optimization the unweighted SIS/SIR versions are comparable to the event-based simulation.
    - The weighted version is slower, but still reasonably fast.
    - There are methods for generic simple contagions and generic complex contagions.
    
  - discrete-time (synchronous update) models
  - tools for visualizing and animating simulated epidemics.
  
- **Numerical solvers for ODE models**

  - pair approximation models
  - effective degree models
  - edge-based compartmental models

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   Getting Started <GettingStarted>
   Examples 
   EoN 
   Changes


.. _repository: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks
.. _Mathematics of epidemics on networks\: from exact to approximate models: http://www.springer.com/us/book/9783319508047
.. _NetworkX: https://networkx.github.io
.. _is to let me know you're using it: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/31

   
