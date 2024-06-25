EoN module
==========

Introduction
------------
**EoN** (Epidemics on Networks) is a Python package for the simulation of 
epidemics on networks and solving ODE models of disease spread.

The algorithms are based on the book
        
`Mathematics of Epidemics on Networks: from Exact to Approximate Models`_
by Kiss, Miller & Simon  (possibly freely available for `download here`_ depending on your institutional subscription).
        
Please cite the book if using these algorithms.

If you use EoN, please 
- cite the 
`Journal of Open Source Software publication <https://doi.org/10.21105/joss.01731>`_
- `leave a note here`_ 

This helps me **get promoted / get funding**, and it makes me happy when 
people use the software.  A happy developer whose job prospects are improving
because of all the people using the software will have more time to improve 
the software.

For simulations, we assume that contact networks are **NetworkX** 
graphs; see https://networkx.github.io/


**EoN** consists of several sets of algorithms.  

- The first deals with **stochastic simulation of epidemics on networks**.  
  The most significant of these are ``fast_SIS`` and ``fast_SIR`` which 
  usually outperform Gillespie algorithms (also included).  These algorithms 
  are discussed in more detail in the appendix of the book.

- A significant extension of these simulations is a set of tools
  designed to **visualize and animate simulated epidemics**, and
  generally help investigate a given stochastic simulation.
  
- Another set deals with **numerical solution of systems of analytic equations** 
  derived in the book.  For these it is possible to either provide the degree 
  distribution, or simply use a network and let the code determine the degree 
  distribution.

- There are a few additional algorithms which are not described in the
  book, but which we believe will be useful. Most notably, related to 
  visualization and generation of animations.

Distributed under MIT license.  See :download:`license.txt<../license.txt>` 
for full details.


Simulation Toolkit
------------------
This submodule deals with epidemic simulation.  We start with a quick list of 
the functions with links to the individual functions.  A brief description is 
below.


Quick list
^^^^^^^^^^

.. currentmodule:: EoN

.. autosummary::
   :toctree: functions/

   fast_SIR
   fast_nonMarkov_SIR
   fast_SIS
   fast_nonMarkov_SIS
   Gillespie_SIR
   Gillespie_SIS
   Gillespie_Arbitrary
   Gillespie_simple_contagion
   Gillespie_complex_contagion
   basic_discrete_SIR
   basic_discrete_SIS
   discrete_SIR
   percolate_network
   directed_percolate_network
   nonMarkov_directed_percolate_network_with_timing
   nonMarkov_directed_percolate_network
   estimate_SIR_prob_size
   estimate_SIR_prob_size_from_dir_perc
   estimate_directed_SIR_prob_size
   estimate_nonMarkov_SIR_prob_size_with_timing
   estimate_nonMarkov_SIR_prob_size
   get_infected_nodes
   percolation_based_discrete_SIR

Short descriptions
^^^^^^^^^^^^^^^^^^
- Event-based algorithms:

  These algorithms use an efficient approach to simulate epidemics.  ``fast_SIR`` 
  and ``fast_SIS`` assume constant transmission and recovery rates, while
  ``fast_nonMarkov_SIR`` and ``fast_nonMarkov_SIS`` allow the user to specify  
  more detailed rules for transmission.
  
  - **fast_SIR**
  - **fast_nonMarkov_SIR** 
  - **fast_SIS**
  - **fast_nonMarkov_SIS**

- Gillespie Algorithms

  These algorithms simulate epidemics assuming constant transmission and 
  recovery rates.  They are commonly used, but in many cases are slower than
  the event driven methods.  I do not see evidence that they are ever 
  significantly faster.  It is not very practical to get away from the 
  constant rate assumptions so I prefer to avoid them.  However, 
  ``Gillespie_simple_contagion`` allows the user to do SEIR, SIRS, or any of a number
  of other more exotic "simple contagion" scenarios that are not in the event-driven
  code.  ``Gillespie_complex_contagion`` handles complex contagions, in which an
  individual requires multiple partners to have a given state before it changes
  status.  For legacy reasons, ``Gillespie_Arbitrary`` is included, it simply calls
  ``Gillespie_simple_contagion``, and will be removed in future versions.
  
  - **Gillespie_SIR**
  - **Gillespie_SIS**
  - **Gillespie_Arbitrary**
  - **Gillespie_simple_contagion**
  - **Gillespie_complex_contagion**
  

- Discrete-time algorithms

  These algirthms are appropriate for where we separate infection into 
  generations.  We assume infection lasts a single time step.  The ``basic_*`` 
  algorithms assume that transmission occurs with probability p for all edges.
  In contrast ``discrete_SIR`` allows for very general user-specified
  transmission rules.
  
  - **basic_discrete_SIR**
  - **basic_discrete_SIS**
  - **discrete_SIR**


- Percolation-based approaches 
    
  There is a close relation between percolation and SIR disease which is
  described in Chapter 6 of the book.  Many of these algorithms are related
  to demonstrating the equivalence as outlined in the book, and are not really
  the most efficient way to simulate an epidemic.  However, these algorithms
  will be useful for estimating probability and size of epidemics. 
    
  - **percolate_network** (undirected percolation corresponding to fixed 
    transmission probability)
  - **directed_percolate_network** (directed percolation corresponding to 
    constant transmission and recovery rates)
  - **nonMarkov_directed_percolate_network_with_timing** (uses user-generated 
    duration and transmission time distributions)
  - **nonMarkov_directed_percolate_network** (uses user-generated transmission 
    rules)
  - **estimate_SIR_prob_size** (estimates prob/size from an undirected percolated network - only appropriate if constant p)
  - **estimate_SIR_prob_size_from_dir_perc** (estimates epi prob and size from a given percolated network)
  - **estimate_directed_SIR_prob_size** (estimates based on constant transmission and recovery rates)
  - **estimate_nonMarkov_SIR_prob_size_with_timing** (estimates based on user-generated transmission and recovery time distributions)
  - **estimate_nonMarkov_SIR_prob_size** (estimates based on user-generated transmission rules)
  - **get_infected_nodes** (simulates epidemic and returns final infected nodes)
  - **percolation_based_discrete_SIR**

Simulation Investigation toolkit
--------------------------------
We can study simulations in detail through the Simulation_Investigation class.
This includes automated generation of animations.

This is particularly useful if we want to look at time series or at animations
of the network as the disease spreads.

When EoN performs a simulation with ``return_full_data`` set to ``True``, it returns
a Simulation_Investigation object.  At it's core, this has the data about when
each node changed status and what its new status became.  This allows us to 
generate plots of the network at any given instance in time and to produce 
animations.

The basic display produced by a Simulation_Investigation object shows the 
network at a given snapshot in time on the left, and on the right it shows the
time series of S, I, and (if SIR) R.  It has the option to add additional 
curves that might have been calculated by an analytic model, or perhaps
another simulation.

In general, any of the dynamic simulations will produce a Simulation_Investigation
object if we pass it ``return_full_data = True``. 

Some examples appear in section :ref:`visualization`.

Quick list
^^^^^^^^^^

.. currentmodule:: EoN.Simulation_Investigation

.. autosummary::
   :toctree: functions/

   display
   animate
   node_history
   node_status
   get_statuses
   summary
   t
   S
   I
   R
   transmissions
   transmission_tree
   add_timeseries
   update_ts_kwargs
   update_ts_label
   update_ts_color_dict 
   update_ts_tex
   sim_update_kwargs
   sim_update_label
   sim_update_color_dict
   sim_update_tex
   set_pos
   
   
   
   
Short description
^^^^^^^^^^^^^^^^^
- Visualizations

  There are two main commands for visualization.  We can either produce a 
  snapshot at a given time, or produce an animation.  In either case we can
  optionally include plots of S, I, (and R) as functions of time.
  
  - **display** (allows us to plot a graph at a specific time
      point, and to optionally include the calculated time series)

  - **animate** (allows us to plot a graph at many time points
      We can create a visualization such as mp4, or save each
      individual frame of an animation.  The time series are optional.
  
  
- Data about simulation

  Often we'll want to be able to check what happened to specific nodes in the
  network, or we'll want to know what the time history of the outbreak looked
  like


  - **node_history**
  - **node_status** returns the status of a node at a given time
  - **get_statuses** returns the status of a collection of nodes at
     a given time (in a dict).
  - **summary**  returns t, S, I, (and if SIR R) for the population
    (or a subset of the population)
  - **t** 
  - **S** 
  - **I**
  - **R**
  - **transmissions** returns a list of 3-tuples of the form (t, u, v) stating
    that u transmitted to v at time t.
  - **transmission_tree** returns a MultiDiGraph where an edge from u to v with
    attribute time = t means that u transmitted to v at time t.  (For SIR this
    is a tree or a forest).
  

- Details for plotting

  The remaining commands are to do with the specifics of how the plots appear
  
  - **update_ts_kwargs**
  - **update_ts_label**
  - **update_ts_colordict**
  - **update_ts_tex**
  - **sim_update_kwargs**
  - **sim_update_label**
  - **sim_update_colordict**
  - **sim_update_tex**
  - **set_pos**

- Plotting of transmission tree



Analytic Toolkit
----------------
This submodule deals with solution to systems of equations appearing in the book.
The majority of these also have a version that take a graph G.  There are 
additional functions that calculate properties which these need.

Quick list
^^^^^^^^^^
.. currentmodule:: EoN

.. autosummary::
   :toctree: functions

   SIS_individual_based
   SIS_individual_based_pure_IC
   SIS_pair_based
   SIS_pair_based_pure_IC
   SIR_individual_based
   SIR_individual_based_pure_IC
   SIR_pair_based
   SIR_pair_based_pure_IC
   SIS_homogeneous_meanfield
   SIR_homogeneous_meanfield
   SIS_homogeneous_pairwise
   SIS_homogeneous_pairwise_from_graph
   SIR_homogeneous_pairwise
   SIR_homogeneous_pairwise_from_graph
   SIS_heterogeneous_meanfield
   SIS_heterogeneous_meanfield_from_graph
   SIR_heterogeneous_meanfield
   SIR_heterogeneous_meanfield_from_graph
   SIS_heterogeneous_pairwise
   SIS_heterogeneous_pairwise_from_graph
   SIR_heterogeneous_pairwise
   SIR_heterogeneous_pairwise_from_graph
   SIS_compact_pairwise
   SIS_compact_pairwise_from_graph
   SIR_compact_pairwise
   SIR_compact_pairwise_from_graph
   SIS_super_compact_pairwise
   SIS_super_compact_pairwise_from_graph
   SIR_super_compact_pairwise
   SIR_super_compact_pairwise_from_graph
   SIS_effective_degree
   SIS_effective_degree_from_graph
   SIR_effective_degree
   SIR_effective_degree_from_graph
   SIR_compact_effective_degree
   SIR_compact_effective_degree_from_graph
   SIS_compact_effective_degree
   SIS_compact_effective_degree_from_graph
   Epi_Prob_discrete
   Epi_Prob_cts_time
   Epi_Prob_non_Markovian
   Attack_rate_discrete
   Attack_rate_discrete_from_graph
   Attack_rate_cts_time
   Attack_rate_cts_time_from_graph
   Attack_rate_non_Markovian
   Attack_rate_discrete
   EBCM_discrete
   EBCM_discrete_from_graph
   EBCM
   EBCM_uniform_introduction
   EBCM_from_graph
   EBCM_pref_mix
   EBCM_pref_mix_from_graph
   EBCM_pref_mix_discrete
   EBCM_pref_mix_discrete_from_graph
   get_Pk
   get_PGF
   get_PGFPrime
   get_PGFDPrime
   estimate_R0

Short description
^^^^^^^^^^^^^^^^^^^^^^^^^^
These come from the book.  The numbers given below are the equation numbers in the book.

- Chapter 3
  
  This chapter deals with models assuming we know the full network structure.
  
    - System (3.7): SIS model: Closes equations by assuming that knowing the probabilities
      for nodes to have each status is enough to predict impact of their interactions
      (ignores temporal correlation between statuses of neighbors).  
    
       - **SIS_individual_based** 
       - **SIS_individual_based_pure_IC**
       
      The pure_IC version assumes that some nodes begin infected with probability
      1 and the others are susceptible with probability 1.
       
    - System (3.26): assumes that tracking pair correlations is enough.  Many
      more equations than individual-based.
      
       - **SIS_pair_based**
       - **SIS_pair_based_pure_IC**
      
    - System (3.30) SIR equivalent of corresponding SIS model.
    
       - **SIR_individual_based**
       - **SIR_individual_based_pure_IC**
      
    - System (3.39) SIR equivalent of corresponding SIS model.
      
       - **SIR_pair_based**
       - **SIR_pair_based_pure_IC**
    
- Chapter 4
    
  This chapter attempts to approximate the exact dynamics by ignoring
  heterogeneity in degree.
    
    - System (4.8) Assumes dynamics determined by average number of contacts
      and number of nodes of each status.
      
       - **SIS_homogeneous_meanfield**
       
    - System (4.9) As for SIS.
      
       - **SIR_homogeneous_meanfield**
       
    - System (4.10) Assumes dynamics are determined by the average number of 
      contacts, nodes of each status, and pairs of each status.
    
       - **SIS_homogeneous_pairwise**
       - **SIS_homogeneous_pairwise_from_graph** (reads properties from input graph)
      
    - System (4.11) 
      
       - **SIR_homogeneous_pairwise**
       - **SIR_homogeneous_pairwise_from_graph**
    
    
- Chapter 5

  This chapter attempts to approximate the exact dynamics and incorporates
  heterogeneity in degree (at the cost of more complex models)
  
    - System (5.10)
      
       - **SIS_heterogeneous_meanfield**
       - **SIS_heterogeneous_meanfield_from_graph**
       
    - System (5.11)
      
       - **SIR_heterogeneous_meanfield**
       - **SIR_heterogeneous_meanfield_from_graph**
       
    - System (5.13)
      
       - **SIS_heterogeneous_pairwise**
       - **SIS_heterogeneous_pairwise_from_graph**
       
    - System (5.15)
      
       - **SIR_heterogeneous_pairwise**
       - **SIR_heterogeneous_pairwise_from_graph**
       
    - System (5.18)
      
       - **SIS_compact_pairwise**
       - **SIS_compact_pairwise_from_graph**
       
    - System (5.19)
      
       - **SIR_compact_pairwise**
       - **SIR_compact_pairwise_from_graph**
       
    - System (5.20)
      
       - **SIS_super_compact_pairwise**
       - **SIS_super_compact_pairwise_from_graph**
       
    - System (5.22)
      
       - **SIR_super_compact_pairwise**
       - **SIR_super_compact_pairwise_from_graph**
       
    - System (5.36)
      
       - **SIS_effective_degree**
       - **SIS_effective_degree_from_graph**
       
    - System (5.38)
      
       - **SIR_effective_degree**
       - **SIR_effective_degree_from_graph**
       
    - System (5.43)
      
       - **SIR_compact_effective_degree**
       - **SIR_compact_effective_degree_from_graph**
       
    - System (5.44)
      
       - **SIS_compact_effective_degree**
       - **SIS_compact_effective_degree_from_graph**
    
- Chapter 6

  This chapter uses percolation-based techniques to explore epidemic properties.
  
    - System (6.2) Given a degree distribution and uniform transmission probability, find epidemic probability.
      
       - **Epi_Prob_discrete**
       
    - System (6.3) As in 6.2, but assuming constant transmission and recovery rates.
      
       - **Epi_Prob_cts_time**
       
    - System (6.5) As in 6.2, but with user-specified transmission rules
      
       - **Epi_Prob_non_Markovian** 
       
    - System (6.6) Given a degree distribution, initial proportion infected, and transmission probability, find attack rate.
      See also System (6.10).
      
       - **Attack_rate_discrete**
       - **Attack_rate_discrete_from_graph**
       
    - System (6.7) as in 6.6, but assuming constant transmission and recovery rates.
      
       - **Attack_rate_cts_time**
       - **Attack_rate_cts_time_from_graph**
       
    - System (6.8) As in 6.6, but with user-specified transmission rules
      
       - **Attack_rate_non_Markovian**
       
    - System (6.10) See code for System (6.6).
       
    - System (6.11) Perform EBCM calculations for discrete-time.
      
       - **EBCM_discrete**
       - **EBCM_discrete_from_graph**
       
    - System (6.12) Perform EBCM calculations for continuous-time.  
      
       - **EBCM** allows initial status to be degree dependant.
       - **EBCM_uniform_introduction** assumes disease introduced at t=0 uniformly at random
       - **EBCM_from_graph** assumes disease introduced at t=0 uniformly at random in
         network of given degree distribution.

    - exercise 6.21  Deals with the EBCM model assuming preferential mixing.
       - **EBCM_pref_mix**
       - **EBCM_pref_mix_from_graph**
       - **EBCM_pref_mix_discrete**
       - **EBCM_pref_mix_discrete_from_graph**

   
Auxiliary Functions
-------------------
We have a few additional functions which are of value.

Quick List
^^^^^^^^^^

.. currentmodule:: EoN

.. autosummary::
   :toctree: functions
   
   get_time_shift
   subsample
   hierarchy_pos


Short Description
^^^^^^^^^^^^^^^^^

    - **get_time_shift** (allows us to shift plots to eliminate the effect of early-time stochasticity)
    - **subsample** (allows us to take output given at a stochastic
      set of times and get output at given times - particularly useful
      to allow for averaging multiple simulations)
    - **hierarchy_pos** (Provides positions that help visualize the transmission 
      chain from a Simulation_Investigation object)
    
    
.. _Mathematics of epidemics on networks\: from exact to approximate models: http://www.springer.com/us/book/9783319508047
.. _leave a note here: https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/issues/31
.. _download here: https://link.springer.com/book/10.1007%2F978-3-319-50806-1