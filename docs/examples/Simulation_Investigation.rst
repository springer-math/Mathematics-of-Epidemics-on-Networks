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
                   return_full_data=True, tmax = 40)
    pos = {node:node for node in G}
    sim.set_pos(pos)
    sim.display(6, node_size = 4) #display time 6
    plt.savefig('SIS_2dgrid.png')

This produces a snapshot at time 6:

.. image:: SIS_2dgrid.png
    :width: 90 %

If we changed the `display` command to have `ts_plots=False` or `ts_plots = []` we get

::

    plt.clf()
    sim.display(6, node_size = 4, ts_plots=[]) #display time 6
    plt.savefig('SIS_2dgrid_no_time_series.png')
    
This produces

.. image:: SIS_2dgrid_no_time_series.png
    :width: 50 %


We can also produce animations.  You may need to install additional software 
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
    sim.set_pos(pos)
    ani=sim.animate(ts_plots=['I', 'SIR'], node_size = 4)  
    ani.save('SIR_2dgrid.mp4', fps=5, extra_args=['-vcodec', 'libx264'])

This produces an animation:

.. raw:: html 

   <video controls src="../_static/SIR_2dgrid.mp4", width = 90%></video> 


Non-SIS or SIR processes
------------------------

If we use a model with other states than `'S'`, `'I'`, and `'R'`, the default 
colors aren't specified.  In this case we need to do a little bit more.

Consider a model where the states are `'Sus'`, `'Inf'`, `'Rec'`, or `'Vac'`.  
That is, an SIR model with vaccination.  We will use `Gillespie_simple_contagion`
for this.  I'm choosing the status names to be longer than one character to 
show changes in the argument `ts_plots` stating what the time-series plots 
should show.

In this model, susceptible people have a rate of becoming vaccinated which is
independent of the disease status.  Otherwise, it is just like the SIR disease
above.  So the "spontaneous transitions" are `'Sux'` to `'Vac'` with rate `0.01`
and `'Inf'` to `'Rec'` with rate `1.0`.  The "induced transitions" are 
`('Inf', 'Sus')` to `('Inf', 'Inf')` with rate `2.0`.

::

    import networkx as nx
    import EoN
    import matplotlib.pyplot as plt
    from collections import defaultdict
    G = nx.grid_2d_graph(100,100) #each node is (u,v) where 0<=u,v<=99
    #we'll initially infect those near the middle 
    initial_infections = [(u,v) for (u,v) in G if 45<u<55 and 45<v<55]

    H = nx.DiGraph()  #the spontaneous transitions
    H.add_edge('Sus', 'Vac', rate = 0.01)
    H.add_edge('Inf', 'Rec', rate = 1.0)
    
    J = nx.DiGraph()  #the induced transitions
    J.add_edge(('Inf', 'Sus'), ('Inf', 'Inf'), rate = 2.0)
    
    IC = defaultdict(lambda:'Sus')
    for node in initial_infections:
        IC[node] = 'Inf'
        
    return_statuses = ['Sus', 'Inf', 'Rec', 'Vac']
    
    colordict = {'Sus': '#009a80','Inf':'#ff2000', 'Rec':'gray','Vac': '#5AB3E6'}
    pos = {node:node for node in G}
    tex = False
    sim_kwargs = {'colordict':colordict, 'pos':pos, 'tex':tex}

    sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax=30, return_full_data=True, sim_kwargs=sim_kwargs)

    sim.display(6, node_size = 4, ts_plots=[['Sus', 'Vac'], ['Inf', 'Rec']])
    plt.savefig('SIRV_display.png')

    ani=sim.animate(ts_plots=[['Inf'], ['Sus', 'Inf', 'Rec', 'Vac']], node_size = 4)  
    ani.save('SIRV_animate.mp4', fps=5, extra_args=['-vcodec', 'libx264'])
    
This produces the image

.. image:: SIRV_display.png
    :width: 90 %
    
and the animation

.. raw:: html 

   <video controls src="../_static/SIRV_animate.mp4", width = 90%></video> 

Note that the labels are in plain text rather than math mode (since `tex=False`)


I have not yet included the ability to easily plot things like `'Sus'+'Vac'`.  Please
submit an issue if you want to do this, or ask a question on stackoverflow with
the tag 'EoN'.
