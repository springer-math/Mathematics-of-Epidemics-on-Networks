Changes from v 1.0
==================

New in v 1.0.8rc4
-----------------
    
    The ``Simulation_Investigation`` code now handles plotting things like
    ``'S+V'`` if we add a time series appropriately.  The last example of
    :ref:`visualization` shows this.
    
New in v 1.0.8rc3
-----------------

    The ``Gillespie_simple_contagion`` and ``Gillespie_complex_contagion`` code 
    can now handle ``return_full_data=True``.
    
    The ``Simulation_Investigation`` objects can now handle arbitrary statuses,
    rather than just SIS and SIR.

    Examples using the ``Simulation_Investigation`` objects have been added or
    updated.
    
New in v 1.0.8rc2
-----------------

    ``Gillespie_simple_contagion`` is now more flexible in how it handles heterogeneity.
    The user can now define a function which will give the 'transmission' rates between
    a pair of nodes and the 'recovery' rates of individual nodes.  So it can be
    more general than the original version.  (a heterogeneous SIRS example is now
    provided)

    Bug fixes in ``basic_discrete_SIS``.

    Changed the discrete SIS and SIR code so that the initial infections occur 
    at t=-1 for the ``simulation_investigation`` objects.
    
    Small change to the default color for Infected nodes (FF2020->FF2000) in 
    simulation_investigation
    


New in v 1.0.7
----------------

   No changes (fixing an error in a tag)

New in v 1.0.6
-----------------

   Documentation for ``Gillespie_complex_contagion`` now includes an example.
   
   Removed print command (left over from debugging) from ``Gillespie_complex_contagion.``
   
New in v 1.0.5
-----------------

   Reintroduced ``Gillespie_Arbitrary`` which just calls ``Gillespie_simple_contagion``
   and provides a warning that it will be discontinued later.
   
   
New in v 1.0.4
-----------------

  
  
  Have added ``Gillespie_complex_contagion`` which can handle complex contagions.
  
  The old ``Gillespie_Arbitrary`` has been renamed ``Gillespie_simple_contagion``.  I 
  have fixed a bug in previous versions that prevented it from handling weighted
  graphs.
  
  
  

  ``Gillespie_Arbitrary`` is now back-compatible to networkx 1.11 (but it has 
  been renamed -- see above). 

  Readthedocs is now providing documentation for each function.
  
  
  

New in v 1.0.3
--------------

  No changes to package, but a small change attempting to get readthedocs to
  correctly build.
    
New in v 1.0.2
--------------
  
  No changes (I accidentally made a typo just before uploading v1.0.1 to pypi
  and I can't reupload with the same name).
  

New in v 1.0.1
--------------

Returning transmission chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When simulations have ``return_full_data=True``, the returned object now includes
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

Currently this does not accept ``return_full_data=True``, and it requires that 
the events all occur as Poisson processes (that is, it makes sense to say 
that there is a rate at which things happen, and that rate depends on the 
status of the nodes and perhaps some property of the node or the partnership, 
but nothing else).

(addresses issues 
`13 <https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/13>`_ 
& `17 <https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/17>`_)


