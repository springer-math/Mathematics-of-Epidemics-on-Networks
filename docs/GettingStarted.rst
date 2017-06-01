Getting Started
===============
EoN consists of two broad types of algorithms.  One set of algorithms is used to simulate epidemics on networks, and the others solve analytic ODE models that attempt to approximate the disease dynamics.

Installation
------------
EoN is available at https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks.

To install EoN, you need the folder EoN and its contents to either be in the directory you are using or somewhere in your Python path.

QuickStart Guide
----------------

The code here provides an example of creating a Barabasi-Albert network.  Then it performs several simulations of an SIR epidemic starting with a fraction rho randomly infected initially.  Finally it uses several analytic models to predict the spread of an epidemic in a random network with the given properties.

::

    import networkx as nx
    import matplotlib.pyplot as plt
    import EoN
    import random
    
    N=10**5
    G=nx.barabasi_albert_graph(N, 5) #create a barabasi-albert graph
    
    tmax = 20
    iterations = 5  #run 5 simulations
    tau = 0.1           #transmission rate
    gamma = 1.0    #recovery rate
    rho = 0.005      #random fraction initially infected
    
    for counter in range(iterations): #run simulations
        initial_infections = random.sample(G.nodes(), int(round(rho*N))) 
        t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds=initial_infections, tmax = tmax)
        plt.plot(t, I, color = 'k', alpha=0.3)
            
    t, S, I, R = EoN.SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
    plt.plot(t, I, '-.', label = 'Homogeneous pairwise', linewidth = 5)
    
    t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
    plt.plot(t, I, ':', label = 'Heterogeneous meanfield', linewidth = 5)
    
    t, S, I, R = EoN.EBCM_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
    plt.plot(t, I, '--', label = 'EBCM approximation', linewidth = 5)
    
    plt.legend(loc = 'upper right')
    plt.savefig('SIR_BA_model_vs_sim.pdf')
    
    plt.clf()
   
    #Now run for SIS.   Simulation is much slower so need smaller network
    N=10**4  
    G=nx.barabasi_albert_graph(N, 5) #create a barabasi-albert graph
    print 'got graph'
    for counter in range(iterations):
        initial_infections = random.sample(G.nodes(), int(round(rho*N))) 
        t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds=initial_infections, tmax = tmax)
        plt.plot(t, I, color = 'k', alpha=0.3)
            
    t, S, I = EoN.SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
    plt.plot(t, I, '-.', label = 'Homogeneous pairwise', linewidth = 5)
    
    t, S, I = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
    plt.plot(t, I, ':', label = 'Heterogeneous meanfield', linewidth = 5)
    
    t, S, I = EoN.SIS_compact_pairwise_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
    plt.plot(t, I, '--', label = 'Compact pairwise', linewidth = 5)
    
    plt.legend(loc = 'lower right')
    plt.savefig('SIS_BA_model_vs_sim.pdf')

EoN Examples
------------

We have collected a number of examples using EoN to generate figures.

Epidemics On Networks figures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are examples to generate (close approximations to) many of the figures in 
`Mathematics of Epidemics on Networks: from Exact to Approximate Models`_. 


* :download:`figure 1.2 <../examples/fig1p2.py>`

* :download:`figure 1.5 <../examples/fig1p5.py>`

* :download:`figure 2.11 <../examples/fig2p11.py>`

* :download:`figure 3.2 <../examples/fig3p2.py>` - (This runs quite slowly, and produces slightly different output for the complete graph compared to the figure in the book.  See comments in code for reasons)

* :download:`figure 4.1 <../examples/fig4p1.py>`

* :download:`figure 4.5 <../examples/fig4p5.py>`

* :download:`figure 4.7 <../examples/fig4p7.py>`  - (Note that the book has a typo.  For (c), $\\tau = 1.1\\tau_c$)

* :download:`figure 4.8 <../examples/fig4p8.py>`

* :download:`figure 4.9 <../examples/fig4p9.py>`

* :download:`figure 4.10 <../examples/fig4p10.py>`

* :download:`figure 4.11 <../examples/fig4p11.py>`  - (Note that the book has a typo.  In fact $\\tau = 1.5\\gamma/<K>$)

* :download:`figure 4.12 <../examples/fig4p12.py>`

* :download:`figure 4.13 <../examples/fig4p13.py>`

For Chapter 5 figures, these examples use larger populations than the figures in the text.

* :download:`figure 5.2 <../examples/fig5p2.py>`  - (Note that the book has a typo.  As with fig 4.7, for (c), $\\tau = 1.1\\tau_c$.  It's worth looking at $1.2\\tau_c$ as well.  It's interesting.)

* :download:`figure 5.3 <../examples/fig5p3.py>`  - (Demonstrates the ease of the X_from_graph versions of the analytic equations)

* :download:`figure 5.4 <../examples/fig5p4.py>`  - (Demonstrates the ease of the X_from_graph versions of the analytic equations)

* :download:`figure 5.5 <../examples/fig5p5.py>`  - (Demonstrates the ease of the X_from_graph versions of the analytic equations)

The remainder of these simulations use reduced sizes or numbers of iterations compared to the published figure.  This is to save time.

* :download:`figure 6.1 <../examples/fig6p1.py>`  - (This also does figure 6.3)

* :download:`figure 6.2 <../examples/fig6p2.py>` 

* figure 6.3 - This is done in the same file as figure 6.1.

* :download:`figure 6.4 <../examples/fig6p4.py>` 

* :download:`figure 6.24 <../examples/fig6p24.py>` 

* :download:`figure 7.2 <../examples/fig7p2.py>` 

* :download:`figure 7.3 <../examples/fig7p3.py>` 

* :download:`figure 7.4 <../examples/fig7p4.py>` 



Additional Examples
^^^^^^^^^^^^^^^^^^^

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