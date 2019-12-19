---
title: 'EoN (Epidemics on Networks): a fast, flexible Python package for simulation, analytic approximation, and analysis of epidemics on networks'

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
 - name: La Trobe University, Melbourne, Australia
   index: 1
 - name: Institute for Disease Modeling, Seattle, Washington, USA
   index: 2
   
date: 19 July 2019

bibliography: paper.bib


---

# Summary

EoN (EpidemicsOnNetworks) is a pure-python package designed to study
infectious processes spreading in networks.  It rose out of the 
book *Mathematics of Epidemics on Networks* [@kiss:EoN], and now consists of over 
100 user-functions.


EoN provides a set of tools for

- Susceptible-Infected-Susceptible (SIS) and Susceptible-Infected-Recovered 
(SIR) disease
  - Stochastic simulation of disease spread in networkx graphs
    - continuous time Markovian
    - continuous time nonMarkovian
    - discrete time
  - Numerical solution of over 20 differential equation models, including
    - individual and pair-based models
    - pairwise models
    - edge-based compartmental models
- Stochastic simulation of a wide range of Simple and Complex contagions
- Visualization and analysis of stochastic simulations

EoN relies on the networkx package [@hagberg2008exploring].
Its documentation is maintained at 
https://epidemicsonnetworks.readthedocs.io/en/latest/ 
including numerous examples at 
https://epidemicsonnetworks.readthedocs.io/en/latest/Examples.html.  In this 
paper we provide brief descriptions of a few of EoN's main tools.  The 
online documentation gives more detail about how to use them, including examples.

## SIR and SIS disease

### Stochastic simulation
The main stochastic simulation tools allow the user to investigate
standard SIS and SIR dynamics (SEIR/SIRS and other processes are addressed 
within the simple contagion model):

- Markovian SIS and SIR simulations  (``fast_SIS``, ``Gillespie_SIS``, ``fast_SIR``, and ``Gillespie_SIR``).
- non-Markovian SIS and SIR simulations (``fast_nonMarkovian_SIS`` and ``fast_nonMarkovian_SIR``).
- discrete time SIS and SIR simulations where infections last a single time step 
  (``basic_discrete_SIS``, ``basic_discrete_SIR``, and ``discrete_SIR``).

For both Markovian and non-Markovian methods it is possible for the transition 
rates to depend on individual or partnership properties.

The continuous-time stochastic simulations have two different implementations: a 
Gillespie implementation [@gillespie1977exact; @doob1945markoff] and an Event-driven
implementation.  They have similar speed if the 
dynamics are Markovian (depending on the network and disease parameters either
may be faster than the other), but the event-driven implementation can also handle 
non-Markovian dynamics.  In earlier versions, the event-driven simulations were 
consistently faster than the Gillespie simulations, leading to the names
`fast_SIR` and `fast_SIS`.  The Gillespie simulations have been optimized
using ideas from [@holme2014model] and [@cota2017optimized].

The algorithms typically handle an SIR epidemic spreading on 
hundreds of thousands of individuals in well under a minute on a laptop. 
SIS versions are slower because the number of events that happen is typically much
larger.

### Differential Equations Models

EoN also provides tools to numerically solve about 20 differential equations
models for SIS or SIR disease spread in networks using Scipy integration tools.
The various models use different information about the network to predict the 
number infected at different times.  Model derivations and explanations of 
their simplifying assumptions are in [@kiss:EoN].

The models require different information about the network structure which can
be provided as inputs.  However, each model has a version which takes an input
network and measures its properties.


## Simple and Complex Contagions

Other contagious processes in networks have received attention.
Many of these can be classifed as either "simple contagions" or "complex 
contagions".

In a "simple contagion" an individual ``u`` may be induced to change status by
an interaction with its partner ``v``.  This status change occurs with the same
rate regardless of the statuses of other partners of ``u`` (there may be a race
between partners to determine which transmits first).  SIS and SIR 
diseases are special cases of simple contagions.

In a "complex contagion" however, the rate at which ``u`` changes 
from one status to another may depend on the statuses of others in some more
complicated way.  Two infected individuals may cause a susceptible individual
to become infected at some higher rate than would result from them acting 
independently.  This is frequently thought to model social contagions where an
individual may only believe something if multiple partners believe it 
[@centola:cascade].  

Simple and complex contagions are currently implemented only in a 
Gillespie setting, and so they require Markovian assumptions.  Although they 
are reasonably fast, it would typically be feasible to make a bespoke algorithm
that runs significantly faster.

### Simple contagions

EoN provides a function ``Gillespie_simple_contagion`` which allows a user to 
specify the rules governing an arbitrary simple contagion.  The implementation 
requires the user to separate out two distinct ways that 
transitions occur: those that are occur spontaneously from an individual's current state
and those that are induced by a partner.  To help demonstrate, consider an 
"SEIR" epidemic, where individuals begin susceptible, but when they interact 
with infectious partners they may enter an exposed state.  They remain in that 
exposed state for some period of time before transitioning into the infectious 
state independently of the status of any partner.
They remain infectious and eventually transition into the recovered state, again
independently of the status of any partner.  Here the "E" to "I" and "I" to "R"
transitions occur spontaneously given the individual's state, while the "S" to "E" 
transition is induced by a partner.  

### Complex contagions

Complex contagions are implemented through ``Gillespie_complex_contagion`` which 
allows a user to specify the rules governing a relatively arbitrary complex 
contagion.  The one criteria we note is that there is no memory - an individual 
will change from one status to another based on the current statuses of its 
neighbors, and not based on previous interactions with some neighbors who may 
have since changed status.

The Gillespie implementation, requires a user-defined function that 
calculates the rate at which ``u`` will change status given the current system 
state and another function which chooses its new status.  It also
needs a function that determines which nodes have their 
rate change due to ``u``'s transition. 

Once these functions are defined, the Gillespie algorithm is able to perform
the complex contagion simulation.


### Visualization & Analysis

By default simulations return numpy arrays providing the counts of each state.  
However if we set a flag ``return_full_data=True``,
then the simulations return a ``Simulation_Investigation`` object.  This provides
access to complete information about the simulation, including the transmission
chains.

The ``Simulation_Investigation`` 
object can create a snapshot of the network at a given time.
By default the visualization includes the time series (e.g., S, I, and R) 
plotted beside the network snapshot, but there is flexibility about what (or if) other
time series appear.  With appropriate additional packages needed for
matplotlib's animation tools, it can produce animations as well.

## Discussion

EoN provides tools for contagious processes spreading in contact networks, including
SIR and SIS disease and more generally simple and complex contagions.  It also
provides tools for visualizing stochastic simulation output.
Full documentation is available at https://epidemicsonnetworks.readthedocs.io/en/latest/

# Dependencies:

scipy
numpy
networkx
matplotlib

# Related Packages

There are several alternative software packages that allow for simulation of 
epidemics on networks.  Here we briefly review some of these.


### epydemic

The Python package Epydemic simulates SIS and SIR epidemics in 
networks.  It is built on networkx.  It handles both discrete-time or 
continuous-time Markovian simulations for which it uses a 
Gillespie-style algorithm.  It can handle any process which can be simulated using  
``EoN.simple_contagion``.

The documentation is available at https://pyepydemic.readthedocs.io/en/latest/


### Graph-tool

The Python package Graph-tool [@peixoto_graph-tool_2014] serves as a networkx
alternative.  Its underlying processes are written in C++, 
so it is often much faster.

Graph-tool has a number of built-in dynamic models, including the SIS, SIR, 
and SIRS models.  The disease models are currently available only in 
discrete-time versions.

The disease model documentation is available at 
https://graph-tool.skewed.de/static/doc/dynamics.html.

### EpiModel

The R package EpiModel [@jenness:EpiModel] can handle SI, SIS, and SIR 
disease spread.  It is possible to extend EpiModel to other processes.  EpiModel
is built around the StatNet package.  More details are available
at https://www.epimodel.org/



# Funding and Support

The development of EoN has been supported by Global Good and by La Trobe 
University.

# References
