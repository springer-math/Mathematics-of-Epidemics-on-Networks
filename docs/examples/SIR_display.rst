SIR Animation Example
---------------------

:download:`Downloadable Source Code <SIR_display.py>` 

The code below produces an animation of an SIR epidemic:

.. raw:: html 

   <video controls src="../_static/SIR_2dgrid.mp4", width = 90%></video> 


You may need to install additional software 
for this to work and modify `extra_args` appropriately.  The commands below 
work on a mac with ffmpeg installed.  The commands below also show an alternate
way to specify the position of nodes passing `pos` through `fast_SIR` to be
used when generating the simulation_investigation object.

::

    import networkx as nx
    import EoN
    import matplotlib.pyplot as plt
    G = nx.grid_2d_graph(100,100) #each node is (u,v) where 0<=u,v<=99
    #we'll initially infect those near the middle 
    initial_infections = [(u,v) for (u,v) in G if 45<u<55 and 45<v<55]
    pos = {node:node for node in G}
    sim_kwargs = {'pos': pos}    
    sim = EoN.fast_SIR(G, 2.0, 1.0, initial_infecteds = initial_infections, 
                   tmax = 40, return_full_data=True, sim_kwargs = sim_kwargs)

    ani=sim.animate(ts_plots=['I', 'SIR'], node_size = 4)  
    ani.save('SIR_2dgrid.mp4', fps=5, extra_args=['-vcodec', 'libx264'])


