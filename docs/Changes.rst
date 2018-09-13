Changes from v 1.0
==================

New in v 1.0.1
-------------


Returning transmission chain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When simulations have `return_full_data=True`, the returned object includes
information on who infected whom at each time.  This can be accessed through: 

`transmissions <functions/EoN.Simulation_Investigation.transmissions.html>`_

`transmission_tree <functions/EoN.Simulation_Investigation.transmission_tree.html>`_.
(note that in an SIS epidemic, this may have cycles)

(addresses issue 21) 

Non-SIS/SIR processes
^^^^^^^^^^^^^^^^^^^^^

It is now possible to run a wide range of non-SIS/SIR processes spreading in
a network.  These processes include competing diseases, SIRS disease, SEIR 
disease, and quite a few other options.  This is done using:

`Gillespie_Arbitrary <functions/EoN.Gillespie_Arbitrary.html>`_.  Currently this
does not accept `return_full_data=True`.

(addresses issue 17)