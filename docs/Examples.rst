EoN Examples
==============

**EoN** consists of two broad types of algorithms.  One set of algorithms is used to simulate epidemics on networks, and the others solve analytic ODE models that attempt to approximate the disease dynamics.


We have collected a number of examples using **EoN** to generate figures.

Epidemics On Networks figures
-----------------------------

Here are examples to generate (close approximations to) many of the figures in 
`Mathematics of Epidemics on Networks: from Exact to Approximate Models`_. 


* :download:`figure 1.2 <../examples/fig1p2.py>`

* :download:`figure 1.5 <../examples/fig1p5.py>`

* :download:`figure 2.11 <../examples/fig2p11.py>`

* :download:`figure 3.2 <../examples/fig3p2.py>` - 
   - In addition to plots in the book's figure, this also plots the average of 1000 simulations.  
   - For the complete graph, the pair equations run quite slowly (there are N choose 2 edges, and we need equations for each).
   - This code does not include the triangle corrections mentioned after system 3.26.

* :download:`figure 4.1 <../examples/fig4p1.py>`

* :download:`figure 4.5 <../examples/fig4p5.py>`

* :download:`figure 4.7 <../examples/fig4p7.py>`  
   - Note that the book has a typo.  For (c), $\\tau = 1.1\\tau_c$

* :download:`figure 4.8 <../examples/fig4p8.py>`

* :download:`figure 4.9 <../examples/fig4p9.py>`

* :download:`figure 4.10 <../examples/fig4p10.py>`

* :download:`figure 4.11 <../examples/fig4p11.py>`  
   - Note that the book has a typo.  In fact $\\tau = 1.5\\gamma/<K>$

* :download:`figure 4.12 <../examples/fig4p12.py>`

* :download:`figure 4.13 <../examples/fig4p13.py>`

For Chapter 5 figures, these examples use larger populations than the figures in the text.

* :download:`figure 5.2 <../examples/fig5p2.py>`  
   - Note that the book has a typo.  As with fig 4.7, for (c), $\\tau = 1.1\\tau_c$. 
   - It's worth looking at $1.2\\tau_c$ as well.  It's interesting.

* :download:`figure 5.3 <../examples/fig5p3.py>`  
   - This and the next 2 demonstrate the ease of the X_from_graph versions of the analytic equations

* :download:`figure 5.4 <../examples/fig5p4.py>`  

* :download:`figure 5.5 <../examples/fig5p5.py>` 

The remainder of these simulations use reduced sizes or numbers of iterations compared to the published figure.  This is to save time.

* :download:`figure 6.1 <../examples/fig6p1.py>`  
   - This also does figure 6.3

* :download:`figure 6.2 <../examples/fig6p2.py>` 

* figure 6.3 
   - This is done in the same file as figure 6.1.

* :download:`figure 6.4 <../examples/fig6p4.py>` 

* :download:`figure 6.24 <../examples/fig6p24.py>` 

* :download:`figure 7.2 <../examples/fig7p2.py>` 

* :download:`figure 7.3 <../examples/fig7p3.py>` 

* :download:`figure 7.4 <../examples/fig7p4.py>` 

For Chapter 9 (nonMarkovian) figures, we have not implemented code that solves 
the dynamic equations but we do have code that will do the simulations.  These 
are given here.

* :download:`figure 9.2 <../examples/fig9p2.py>` 

* :download:`figure 9.4 <../examples/fig9p4.py>` 

* :download:`figure 9.5 <../examples/fig9p5.py>` 


Additional Examples
-------------------

Sometimes you might want to have the values of parameters change at different 
times.

* :download:`SIS varying tau <../examples/changing_parameters/SIS_change_tau.py>`

* :download:`SIR varying tau <../examples/changing_parameters/SIR_change_tau.py>`

You may have edges (or nodes) with weights affecting transmission or recovery
rates.

* :download:`SIS weighted edges <../examples/weighted_graph/SIS_weighted.py>`

Are you trying to do something but can't figure it out and would like an example?  

`Submit an issue`_ and I'll try to help.


.. _Mathematics of epidemics on networks\: from exact to approximate models: http://www.springer.com/us/book/9783319508047
.. _Submit an issue: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues

