SEIR
----

:download:`Downloadable Source Code <arbitrary_dynamics/SEIR.py>` 

.. image:: arbitrary_dynamics/SEIR.png
    :width: 80 %

::


    import EoN
    import networkx as nx
    from collections import defaultdict
    import matplotlib.pyplot as plt
    
    N = 100000
    G = nx.fast_gnp_random_graph(N, 5./(N-1))
    
    H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
    H.add_node('S') #It would still work if this node weren't included.
    H.add_edge('E', 'I', rate = 0.6)   #E->I
    H.add_edge('I', 'R', rate = 0.2)   #I->R
    
    J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
    J.add_edge(('I', 'S'), ('I', 'E'), rate = 0.1) #IS -> EI
    
    IC = defaultdict(lambda: 'S')
    for node in range(200):
        IC[node] = 'I'
    
    return_statuses = ('S', 'E', 'I', 'R')
    
    t, S, E, I, R = EoN.Gillespie_Arbitrary(G, H, J, IC, return_statuses, 
                                            tmax = float('Inf'))    
        
    plt.plot(t, S, label = 'Susceptible') 
    plt.plot(t, E, label = 'Exposed') 
    plt.plot(t, I, label = 'Infected')  
    plt.plot(t, R, label = 'Recovered') 
    plt.legend()
    plt.savefig('SEIR.png')