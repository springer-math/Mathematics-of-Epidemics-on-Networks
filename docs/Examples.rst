EoN Examples
==============


We have collected a number of examples using **EoN** to generate figures.    We 
start with examples from the book 
`Mathematics of Epidemics on Networks: from Exact to Approximate Models`_.  Then
we give a few other examples.  If none of these examples helps with the particular
problem you are facing, `Submit an issue`_ or go to `stackoverflow`_ and use the 'eon' tag.

If you have an example you think would be useful here, please `submit an issue`_
or email me.

Reproducing figures from "Mathematics of Epidemics on Networks"
---------------------------------------------------------------

Here are examples to generate (close approximations to) many of the figures in 
`Mathematics of Epidemics on Networks: from Exact to Approximate Models`_. 

Chapter 1
^^^^^^^^^
Introduction

.. toctree::
   :maxdepth: 2

	      
   Figure 1.2 <examples/fig1p2.rst>

   Figure 1.5 <examples/fig1p5.rst>

Chapter 2
^^^^^^^^^
Top down models

.. toctree::
   :maxdepth: 2

   Figure 2.11 <examples/fig2p11.rst>

Chapter 3
^^^^^^^^^
Bottom up models

.. toctree::
   :maxdepth: 2

   Figure 3.2 <examples/fig3p2.rst>

Chapter 4
^^^^^^^^^
Homogeneous meanfield

.. toctree::
   :maxdepth: 2
   
   Figure 4.1 <examples/fig4p1.rst>
   Figure 4.5 <examples/fig4p5.rst>
   Figure 4.7 <examples/fig4p7.rst>
   Figure 4.8 <examples/fig4p8.rst>
   Figure 4.9 <examples/fig4p9.rst>
   Figure 4.10 <examples/fig4p10.rst>
   Figure 4.11 <examples/fig4p11.rst>
   Figure 4.12 <examples/fig4p12.rst>
   Figure 4.13 <examples/fig4p13.rst>


Chapter 5
^^^^^^^^^
Heterogeneous meanfield

For Chapter 5 figures, these examples use larger populations than the figures in 
the text.

Figures 5.3, 5.4, and 5.5 demonstrate the ease of the X_from_graph versions of 
the analytic equations

.. toctree::
   :maxdepth: 2

   Figure 5.2 <examples/fig5p2.rst>
   Figure 5.3 <examples/fig5p3.rst>
   Figure 5.4 <examples/fig5p4.rst>
   Figure 5.5 <examples/fig5p5.rst>
   
Chapter 6
^^^^^^^^^
Percolation and EBCM

The remainder of these simulations use reduced sizes or numbers of iterations 
compared to the published figure.  This is to save time.

.. toctree::
   :maxdepth: 2
 
   Figure 6.1 (and 6.3) <examples/fig6p1.rst>
   Figure 6.2 <examples/fig6p2.rst>
   Figure 6.3 is done in the same code as 6.1. <examples/fig6p1.rst>
   Figure 6.4 <examples/fig6p4.rst>
   Figure 6.24 <examples/fig6p24.rst>
   

Chapter 7
^^^^^^^^^
Model hierarchies

.. toctree::
   :maxdepth: 2
	      
   Figure 7.2 <examples/fig7p2.rst>
   Figure 7.3 <examples/fig7p3.rst>
   Figure 7.4 <examples/fig7p4.rst>


Chapter 9
^^^^^^^^^
Non-Markovian processes

For Chapter 9 (nonMarkovian) figures, we have not implemented code that solves 
the dynamic equations but we do have code that will do the simulations.  These 
are given here.

.. toctree::
   :maxdepth: 2
 
   Figure 9.2 <examples/fig9p2.rst>
   Figure 9.4 <examples/fig9p4.rst>
   Figure 9.5 <examples/fig9p5.rst>



Additional Examples
-------------------


.. _visualization:


Visualizing or animating disease spread
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can visualize snapshots or animations of disease spread in a network.
For these examples, we'll take a 100x100 grid of nodes [each node is (i,j)]
connected to their 4 nearest neighbors (except the nodes on the edges).  
This isn't the most realistic network, but it is a good example for showing 
the automatic plotting tools.

.. toctree:: 
   :maxdepth: 1
   
   Visualizing a snapshot of SIS disease spread in the network <examples/SIS_display.rst>
   Visualizing an animation of SIR disease spread in the network <examples/SIR_display.rst>
   Visualizing non-SIS/SIR processes (SIRV) spreading in the network <examples/SIRV_display.rst>
   
The last example shows how we might plot things like ``'I'+'R'`` or other combinations of
the data.
   


Non-uniform infectiousness
^^^^^^^^^^^^^^^^^^^^^^^^^^

Perhaps the population is made up of different types of individuals with different
infectiousness.  This example considers a heterosexual population where one group
is more infectious than the other.  This example includes a user-defined transmission
function and uses the simulation investigation tools to look at the dynamics in
each population.

.. toctree::
   :maxdepth: 1

   Bipartite population <examples/bipartite.rst>


Weighted networks
^^^^^^^^^^^^^^^^^

You may have edges (or nodes) with weights affecting transmission or recovery
rates.  (for this be aware that the syntax of edge/node attributes is different
for networkx 2.x and 1.x).

.. toctree::
   :maxdepth: 1

   SIS weighted edges <examples/SIS_weighted.rst>

.. _simple-contagion-section:

Non-SIS/SIR processes with  `Gillespie_simple_contagion <functions/EoN.Gillespie_simple_contagion.html>`_ 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default Gillespie and Event-driven algorithms assume SIS or SIR dynamics.
However, you may want something else, such as SEIR or SIRS or maybe more exotic
things.  If you are willing to assume that events all happen as Poisson processes
(that is, an event happens at a rate), then it's possible to do this with a
Gillespie approach.
 
.. toctree::
   :maxdepth: 1
   
   SEIR <examples/SEIR.rst>
   SIRS <examples/SIRS.rst>
   SIRS with heterogeneity <examples/SIRS_heterogeneous.rst>
   Competing SIR Diseases <examples/Compete.rst>
   Cooperative SIR Diseases <examples/Cooperate.rst>
   


Other
^^^^^

Are you trying to do something but can't figure it out and would like an example?  

`Submit an issue`_ or go to `stackoverflow`_ and use the 'eon' tag.
I'll try to help.

If you have developed something that you think would make a good example and you'd
like to share it, please let me know.

.. _Mathematics of epidemics on networks\: from exact to approximate models: http://www.springer.com/us/book/9783319508047
.. _Submit an issue: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues
.. _stackoverflow: https://stackoverflow.com/questions/tagged/eon?tab=Newest
