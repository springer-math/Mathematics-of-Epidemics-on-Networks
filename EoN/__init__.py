r'''
EoN (Epidemics on Networks)

EoN is a Python package for the simulation of epidemics on networks 
and ODE models of disease spread.

The algorithms are based on the book
        
`Mathematics of epidemics on networks: from exact to approximate 
models`
by Kiss, Miller & Simon
        http://www.springer.com/book/9783319508047
        
Please cite the book if using these algorithms

For simulations, we assume that input networks are **NetworkX** 
graphs; see https://networkx.github.io/



EoN consists of two sets of algorithms.  

- The first deals with simulation of epidemics on networks.  The most significant of these are `fast_SIS` and `fast_SIR` which significantly outperform Gillespie algorithms (also included).  These algorithms are discussed in more detail in the appendix of the book.


- The second deals with solution of systems of equations derived in the book.  For these it is possible to either provide the degree distribution, or simply use a network and let the code determine the degree distribution.


- There are a few additional algorithms which are not described in the book, but which we believe will be useful. Most notably, the some of the visualization/animation commands.

Distributed under MIT license.  See :download:`license.txt<../license.txt>` for full details.


Auxiliary functions
-------------------
We start with a few useful auxiliary functions

'''

__author__ = "Joel C. Miller, Istvan Z. Kiss, and Peter Simon"
__version__ = "1.0.3"


#__all__ = 

class EoNError(Exception):
    r'''
    this will be the basic error type for EoN
    '''
    pass

def _get_rate_functions_(G, tau, gamma, transmission_weight = None, 
                        recovery_weight=None):
    r'''
    Arguments : 
        G : networkx Graph
            the graph disease spread on

        tau : number
            disease parameter giving edge transmission rate (subject to edge scaling)

        gamma : number (default None)
            disease parameter giving typical recovery rate, 
        
        transmission_weight : string (default None)
            `G.edge[u][v][transmission_weight]` scales up or down the recovery rate.

        recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                `gamma_i = G.node[i][recovery_weight]*gamma`
    Returns : 
        : trans_rate_fxn, rec_rate_fxn
            Two functions such that 
            - `trans_rate_fxn(u,v)` is the transmission rate from u to v and
            - `rec_rate_fxn(u)` is the recovery rate of u.
'''
    if transmission_weight is None:
        trans_rate_fxn = lambda x, y: tau
    else:
        try:
            trans_rate_fxn = lambda x, y: tau*G.edges[x,y][transmission_weight]
        except AttributeError: #apparently you have networkx v1.x not v2.x
            trans_rate_fxn = lambda x, y: tau*G.edge[x][y][transmission_weight]

    if recovery_weight is None:
        rec_rate_fxn = lambda x : gamma
    else:
        rec_rate_fxn = lambda x : gamma*G.node[x][recovery_weight]


    return trans_rate_fxn, rec_rate_fxn




import EoN.auxiliary
from EoN.auxiliary import *
import EoN.simulation
from EoN.simulation import *
import EoN.analytic
from EoN.analytic import *
import EoN.simulation_investigation
from EoN.simulation_investigation import *

'''
These are the systems I still want to include:

(8.1) SIS pairwise contact conserving rewiring
(8.5) SIS eff. deg. contact conserving rewiring
(8.7) SIS pairwise random activation/deletion
(8.13) SIS eff. deg. random activation/deletion
(8.15) SIS pairwise link-status dependent act/del
(8.16) SIS link deactivation-activation on fixed networks.
(8.19) EBCM dynamic network

(9.5) SI^{K}R multistage pairwise for homogeneous
(9.27) SIR pairwise, constant infection duration.
(9.35) SIR homogeneous pairwise, general recovery
(9.36) SIR EBCM non-Markovian trans/recovery

'''

