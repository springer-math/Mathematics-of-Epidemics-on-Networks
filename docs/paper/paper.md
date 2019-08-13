---
title: 'EoN (Epidemics on Networks), software for simulation, analytic approximation, and analysis of epidemics on networks.'
tags:
  - SIR
  - SIS
  - Epidemic
  - Network
  - Stochastic simulation
  - Python
authors:
 - name: Joel C. Miller
   orcid: 0000-0003-4426-0405
   affiliation: "1, 2"
 - name: Tony Ting
   affiliation: "2"
affiliations:
 - name: La Trobe University
   index: 1
 - name: Institute for Disease Modeling
   index: 2
date: 19 July 2019
bibliography: paper.bib
---

# Summary

EoN (EpidemicsOnNetworks) is a pure-python package designed to assist studies of
infectious processes spreading through networks.  It originally rose out of the 
book Mathematics of Epidemics on Networks [], and now consists of over 100 
user-functions.

The main functions provided allow the user to perform and analyze stochastic 
simulations on networkx graphs.  This set of ``simulation'' tools allows the user
to perform:

- Markovian SIS and SIR simulations
- non-Markovian SIS and SIR simulations
- discrete time SIS and SIR simulations where infections last a single time step
- a wide range of Markovian ``simple'' contagions and Markovian ``complex'' contagions.

These algorithms permit transition rates to depend on properties of nodes and
of edges.



perform generation-based SIS and SIR 

EoN provides a set of analytic tools as well as a set of simulation tools.

The "analytic" tools allow us to solve a 
number of ordinary differential equations (ODE) models of SIS 
(Susceptible-Infected-Susceptible) and SIR (Susceptible-Infected-Recovered) 
disease spreading through networks.  These use the Scipy integration tools.
The derivations of the models are described in [].

The "simulation" tools allow us to simulate infection spreading stochastically
across a static network.  The tools for continuous-time simulations come in two
types: event-based and Gillespie.  Both have undergone significant effort to 
make them efficient.  They can typically handle an SIR epidemic spreading on 
hundreds of thousands of nodes within a few seconds on a laptop.  The SIS 
versions are slower because the epidemic may not die out.  

The event-based tools allow for non-Markovian processes, but 
are restricted to SIS and SIR processes.  In earlier versions, these were faster
than the Gillespie tools, and thus they are named `fast_SIR` and `fast_SIS`.

The Gillespie tools make assumptions that the processes are Markovian.  Through
approaches described in [] they are able to peform at similar speed to the
event-based tools.  Although they cannot handle non-Markovian processes, they
are able to handle a much wider variety of spreading processes.  In addition
to ``Gillespie_SIS`` and ``Gillespie_SIR``, there is ``Gillespie_simple_contagion``
and ``Gillespie_complex_contagion`` which allow the user to define new transmission
processes.  Examples are provided in the documentation, including
- SEIR disease (there is an exposed state before becoming infectious)
- SIRS disease (recovered individuals eventually become susceptible again)
- SIRV disease (individuals may get vaccinated) 
- Competing SIR disease
- Cooperative SIR disease

By default the simulations return numpy arrays providing the number of individuals
with each state at each time.  However by setting a flag ``return_full_data=True``,
we can know the exact status of each node at each time, as well as who infected
whom.  There are also methods which use this data to visualize the epidemic at 
a specific time, or to create an animation.  


::

    import networkx as nx
    import EoN
    import matplotlib.pyplot a plt
    from collections import defaultdict

    G = nx.    
    
    #In the below:
    #'SS' means a node susceptible to both diseases
    #'SI' means susceptible to disease 1 and infected with disease 2
    #'RS' means recovered from disease 1 and susceptible to disease 2.
    #etc.
    
    H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
    H.add_node('SS')
    H.add_edge('SI', 'SR', rate = 1)
    H.add_edge('IS', 'RS', rate = 1)
    H.add_edge('II', 'IR', rate = 1)
    H.add_edge('II', 'RI', rate = 1)
    H.add_edge('IR', 'RR', rate = 0.5)
    H.add_edge('RI', 'RR', rate = 0.5)
    
    #In the below the edge (('SI', 'SS'), ('SI', 'SI')) means an
    #'SI' node connected to an 'SS' node can lead to a transition in which
    #the 'SS' node becomes 'SI'.  The rate of this transition is 0.2.
    #
    #Note that `IR` and `RI` nodes are more infectious than other nodes.
    #
    J = nx.DiGraph()    #DiGraph showing transitiona that do require an interaction.
    J.add_edge(('SI', 'SS'), ('SI', 'SI'), rate = 0.2)
    J.add_edge(('SI', 'IS'), ('SI', 'II'), rate = 0.2)
    J.add_edge(('SI', 'RS'), ('SI', 'RI'), rate = 0.2)
    J.add_edge(('II', 'SS'), ('II', 'SI'), rate = 0.2)
    J.add_edge(('II', 'IS'), ('II', 'II'), rate = 0.2)
    J.add_edge(('II', 'RS'), ('II', 'RI'), rate = 0.2)
    J.add_edge(('RI', 'SS'), ('RI', 'SI'), rate = 1)
    J.add_edge(('RI', 'IS'), ('RI', 'II'), rate = 1)
    J.add_edge(('RI', 'RS'), ('RI', 'RI'), rate = 1)
    J.add_edge(('IS', 'SS'), ('IS', 'IS'), rate = 0.2)
    J.add_edge(('IS', 'SI'), ('IS', 'II'), rate = 0.2)
    J.add_edge(('IS', 'SR'), ('IS', 'IR'), rate = 0.2)
    J.add_edge(('II', 'SS'), ('II', 'IS'), rate = 0.2)
    J.add_edge(('II', 'SI'), ('II', 'II'), rate = 0.2)
    J.add_edge(('II', 'SR'), ('II', 'IR'), rate = 0.2)
    J.add_edge(('IR', 'SS'), ('IR', 'IS'), rate = 1)
    J.add_edge(('IR', 'SI'), ('IR', 'II'), rate = 1)
    J.add_edge(('IR', 'SR'), ('IR', 'IR'), rate = 1)
    
    G = nx.grid_2d_graph(100,100) #each node is (u,v) where 0<=u,v<=99
    #we'll initially infect those near the middle
    initial_infections = [(u,v) for (u,v) in G if 45<u<55 and 45<v<55]
    IC = defaultdict(lambda: 'SS')
    for node in initial_infections:
        if u+v<100:
            IC[node] = 'SI'
        else:
            IC[node] = 'IS'
            
    colordict = {'SS': '#009a80','SI':'#ff2000', 'IS': '#5AB3E6'}
    pos = {node:node for node in G}
    tex = False
    sim_kwargs = {'colordict':colordict, 'pos':pos, 'tex':tex}

    return_statuses = ('SS', 'SI', 'SR', 'IS', 'II', 'IR', 'RS', 'RI', 'RR')
    sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax=30, return_full_data=True, sim_kwargs=sim_kwargs)

[Examples here, include cooperative diseases + animation]

The documentation is provided at [XX].  This includes over 30 worked examples.


Dependencies:


# Funding and Support

# References