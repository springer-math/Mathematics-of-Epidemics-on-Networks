..  _Simulation_Investigation:

Detailed Simulation Investigation
=================================

When EoN performs a simulation with `return_full_data` set to True, it returns
a Simulation_Investigation object.  At it's core, this has the data about when
each node changed status and what its new status became.  This allows us to 
generate plots of the network at any given instance in time and to produce 
animations.

The basic display produced by a Simulation_Investigation object shows the 
network at a given snapshot in time on the left, and on the right it shows the
time series of S, I, and (if SIR) R.  It has the option to add additional 
curves as might have been calculated by an analytic model, or perhaps
another simulation.

In general, any of the dynamic simulations will produce a Simulation_Investigation
object if we pass it `return_full_data = True`.  

Quick Examples (including animation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For these examples, we'll take a 100x100 grid of nodes [each node is (i,j)]
connected to their nearest neighbors.  This isn't the most realistic network, 
but it is a good example for showing the automatic plotting tools.

::

    import networkx as nx
    import EoN
    import matplotlib.pyplot as plt
    G = nx.grid_2d_graph(100,100) #each node is (u,v) where 0<=u,v<=99
    #we'll initially infect those near the middle 
    initial_infections = [(u,v) for (u,v) in G if 45<u<55 and 45<v<55]
    sim = EoN.fast_SIS(G, 1.0, 1.0, initial_infecteds = initial_infections, 
                   return_full_data=True, tmax = 10)
    pos = {node:node for node in G}
    sim.set_pos(pos)
    sim.display(6, node_size = 4) #display time 6
    plt.savefig('SIS_2dgrid.png')

This produces a snapshot at time 6:

.. image:: SIS_2dgrid.png
    :width: 90 %

::

    import networkx as nx
    import EoN
    import matplotlib.pyplot as plt
    G = nx.grid_2d_graph(100,100) #each node is (u,v) where 0<=u,v<=99
    #we'll initially infect those near the middle 
    initial_infections = [(u,v) for (u,v) in G if 45<u<55 and 45<v<55]
    sim = EoN.fast_SIR(G, 2.0, 1.0, initial_infecteds = initial_infections, 
                   return_full_data=True, tmax = 10)
    pos = {node:node for node in G}
    sim.set_pos(pos)
    ani=sim.animate(ts_plots=['I', 'SIR'], node_size = 4) 
    ani.save('SIR_2dgrid.mp4', fps=5, extra_args=['-vcodec', 'libx264'])

This produces an animation:

.. raw:: html 

   <video controls src="../_static/SIR_2dgrid.mp4", width = 90%></video> 

