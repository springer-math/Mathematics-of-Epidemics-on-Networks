Changes from v 1.0
==================

New in v 1.0.1
--------------

Returning transmission chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When simulations have `return_full_data=True`, the returned object now includes
information on who infected whom at each time.  This can be accessed through: 

`transmissions <functions/EoN.Simulation_Investigation.transmissions.html>`_
which returns a list of tuples (t,u,v) stating that node u infected node v at 
time t.

`transmission_tree <functions/EoN.Simulation_Investigation.transmission_tree.html>`_
which returns a directed multi graph where an edge from u to v with attribute 'time' 
equal to t means u infected v at time t.

(note that in an SIS epidemic, this "tree" may have cycles and repeated edges)

(addresses issue `21 <https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/21>`_ )

Non-SIS/SIR processes
^^^^^^^^^^^^^^^^^^^^^

It is now possible to run a wide range of non-SIS/SIR processes spreading in
a network.  These processes include competing diseases, SIRS disease, SEIR 
disease, and quite a few other options.  This is done using:

`Gillespie_Arbitrary <functions/EoN.Gillespie_Arbitrary.html>`_.  

Examples are `here <Examples.html#non-sis-sir-processes-with-gillespie-arbitrary>`_.

Currently this does not accept `return_full_data=True`, and it requires that 
the events all occur as Poisson processes (that is, it makes sense to say 
that there is a rate at which things happen, and that rate depends on the 
status of the nodes and perhaps some property of the node or the partnership, 
but nothing else).

(addresses issues 
`13 <https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/13>`_ 
& `17 <https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/17>`_)


.. 
  New in v 1.0.2
.. 
  --------------
  
  No changes (I accidentally made a typo just before uploading v1.0.1 to pypi
  and I can't reupload with the same name).
  
  
..
  New in v 1.0.3
..
  --------------
  No changes to package, but a small change attempting to get readthedocs to
  correctly build.