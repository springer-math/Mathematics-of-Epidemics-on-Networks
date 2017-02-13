r'''
EoN

"Epidemics on Networks"
======

    EoN is a Python package for the simulation of epidemics on networks and ODE models of disease spread.

    The algorithms are based on the book:
            Mathematics of epidemics on networks: from exact to approximate models
            - Kiss, Miller & Simon
            - More information at https://EpidemicsOnNetworks.github.io/EpidemicsOnNetworks/
        Please cite the book if using these algorithms

    For simulations, we assume that input networks are **NetworkX** graphs; see https://networkx.github.io/



This is a preliminary version of the code:
  - The algorithms have not been tested in Python 3 (tested only in 2.7)

  - The figure/equation references are based on the current draft of the book, so numbers may change

  - I believe the simulations all work.

  - At present the ODE models are not fully tested.

  - Additional algorithms may be added beyond those in the book.

The following are not complete:
SIS_pair_based, order of arguments for the two SIR versions and in documentation
Attack_rate_non_Markovian

_SIR_pair_based_initialize_* documentation
SIS_homogeneous_pairwise_from_graph documentation
SIR_homogeneous_pairwise_from_graph documentation
many of the *_from_graph do not have documentation explaining inputs/outputs
many of the descriptions do not state what they return.
-------
Distributed under MIT license.  See license.txt for full details.

'''


import networkx as nx
from collections import defaultdict, Counter
import random
import heapq
import random
import scipy
from scipy import integrate
from scipy.ndimage.interpolation import shift


#######################
#                     #
#   Auxiliary stuff   #
#                     #
#######################
class EoNError(Exception):
    #this will be the basic error type for EoN
    pass

class Event(object): #for fast_SIR and fast_SIS
    r'''
    This class is used in event-driven simulations (fast_SIR and fast_SIS)
    as an event which will be put into a priority queue.  It is sortable
    based on event time.
    '''
    def __init__(self, time, action, node, source = None):
        self.time = time
        self.action = action
        self.node = node
        self.source = source #not needed for fast_SIR or recoveries

    def __lt__(self, other): #used to sort Q
        return self.time < other.time
    

class _ListDict_(object):
    r'''
    The Gillespie algorithm with rejection-sampling will involve a step that
    samples a random element from the set.  This is slow in Python.  So I'm
    introducing a new class based on a stack overflow answer by
    Amber (http://stackoverflow.com/users/148870/amber) for a question by
    tba (http://stackoverflow.com/users/46521/tba) found at
    http://stackoverflow.com/a/15993515/2966723

    Based on some limited tests (in some perhaps-atypical networks),
    the benefit appears to be pretty small.  It may be worth creating this
    data structure in C, but I doubt it.
    '''
    
    def __init__(self):
        self.item_to_position = {}
        self.items = []

    def __len__(self):
        return len(self.items)

    def add(self, item):
        if item in self.item_to_position:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove(self, item):
        position = self.item_to_position.pop(item)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position

    def choose_random(self):
        return random.choice(self.items)



def my_odeint(dfunc, V0, times, args=()):
    r'''For some of the systems odeint will switch to the BDF solver.
    In large enough systems, it then gets stuck trying to estimate the Jacobian.

    This routine has identical inputs to integrate.odeint, but relies on integrate.ode.  It avoids BDF.

    In particular, this seems to be important for SIS heterogeneous pairwise where the number of equations is very large.
    I have found that near equilibrium, this often is interpreted as being a stiff system and it switches to bdf, which
    requires calculating a Jacobian.  In some systems this is impractically large.
    
    See this question: http://stackoverflow.com/q/40317096/2966723,
    with the answer by Phillip: http://stackoverflow.com/users/1881610/phillip
    '''

    r = integrate.ode(lambda t, X: dfunc(X, t, *args))
    r.set_integrator('vode', method='adams')
    r.set_initial_value(V0,times[0])
    V=[V0]
    for time in times[1:]:
        V.append(r.integrate(time))
    V = scipy.array(V)
    return V


def get_Nk_and_IC_as_arrays(G, rho, SIR=True):
    r'''
    Input
    -----
    G : networkx graph
    rho : number between 0 and 1
          fraction of nodes to infect at time 0.
    SIR : boolean
          says whether the system will be SIR or SIS.

    Returns
    -------
    Nk : scipy array
         NUMBER (not proportion) of nodes of each degree.
    Sk0 : scipy array
          NUMBER of susceptible nodes of each degree at t=0, = (1-rho)Nk
    Ik0 : scipy array    
          NUMBER of infected nodes of each degree at t=0,   = rho Nk
    if SIR, also returns
    Rk0 : scipy array
          NUMBER of recovered nodes of each degree at t=0,    = 0 Nk
    '''
    
    degree_count = Counter(G.degree().values())
    maxk = max(degree_count.keys())
    
    Nk = scipy.array([degree_count[k] for k in range(maxk+1)])
    Sk0 = (1-rho)*Nk
    Ik0 = rho*Nk
    Rk0 = 0*Nk
    
    if SIR:
        return Nk, Sk0, Ik0, Rk0
    else:
        return Nk, Sk0, Ik0

def get_NkNl_and_IC_as_arrays(G, rho, withKs = False, SIR=True):
    r'''
    Input
    -----
    G : networkx graph
    rho : number between 0 and 1
          fraction of nodes to infect at time 0.
    withKs : boolean
             flag to say whether we are restricting our attention to just
             those degrees observed in the network or to all degrees.
             If True, then we only consider those degrees that are observed.
             If False, then we treat it as if all degrees from 0 to kmax are observed.
    SIR : boolean
          says whether the system will be SIR or SIS.

    Returns
    -------
    NkNl : 2D scipy array
         NUMBER (not proportion) of edges between each pair of degrees.
    SkSl0 : 2D scipy array
          initial NUMBER of edges between pair of susceptibel nodes of each degree type.
          = (1-rho)^2 NkNl
    SkIl0 : 2D scipy array    
          initial NUMBER of edges from a susceptible to an infected node of the given degrees.
          = rho(1-rho) NkNl

    if not SIR, also returns
    IkIl0 : 2D scipy array
          initial NUMBER of edges between 2 infected nodes.  This is not needed for SIR model.
          = rho^2*NkNl

    if withKs, also returns
    Ks : scipy array
         The observed degrees in the population.
    ''' 
    if withKs:
        Ks = sorted(list(set(G.degree().values())))
        klength = len(Ks)
    else:
        klength = max(G.degree().values())+1
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))
    for u,v in G.edges():
        k = G.degree(u)
        l = G.degree(v)
        NkNl[Ks.index(k)][Ks.index(l)] += 1
        NkNl[Ks.index(l)][Ks.index(k)] += 1
    SkSl0 = (1-rho)*(1-rho)*NkNl
    SkIl0 = (1-rho)*rho*NkNl
    IkIl0 = rho*rho*NkNl
    if withKs:
        if SIR:
            return NkNl, SkSl0, SkIl0, scipy.array(Ks)
        else:
            return NkNl, SkSl0, SkIl0, IkIl0, scipy.array(Ks)
    else:
        if SIR:
            return NkNl, SkSl0, SkIl0
        else:
            return NkNl, SkSl0, SkIl0, IkIl0

def get_Pk(G):
    r'''Used in several places so that we can input a graph and then we can call the methods that depend on the degree distribution

    INPUTS
    ------
    G : networkx Graph

    RETURNS
    -------
    Pk : dict
         Pk[k] is the proportion of nodes with degree k.
    '''

    degree_count = Counter(G.degree().values())
    Pk = {x:degree_count[x]/float(G.order()) for x in degree_count.keys()}
    return Pk

def get_Psi(Pk):
    r'''
    INPUTS
    ------
    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    RETURNS
    psi : function.
          psi(x) = \sum_k Pk[k] x^k
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return lambda x: Pkarray.dot(x**k)

def get_PsiPrime(Pk):
    r'''
    INPUTS
    ------
    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    RETURNS
    psiPrime : function.
          psi'(x) = \sum_k k Pk[k] x^{k-1}
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])

    return lambda x: Pkarray[k].dot(k*x**(k-1))

def get_PsiDPrime(Pk):
    r'''
    INPUTS
    ------
    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    RETURNS
    psiDPrime : function.
          psi''(x) = \sum_k k(k-1)Pk[k] x^{k-2}
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return lambda x: Pkarray[k].dot(k*(k-1)*x**(k-2))



def subsample(report_times, times, status1, status2=None, status3 = None):
    r'''
    If given S, I, and/or R, returns them subsampled at specific times.  
    if more than one is given, does so as a list in order given, but skipping whichever
        was not included (if any not included)
    If only one is given then returns just that.

    if report_times goes longer than times, then this simply assumes the system freezes in the final state.
    There is probably a better way to refactor this.

    INPUTS
    ------
    report_times : iterable (ordered)
                   times at which we want to know state of system
    times : iterable (ordered)
            times at which we have the system state (assumed no change between these times)
    statusX (X one of 1, 2 or 3) : iterable (order corresponds to times)
                          generally S, I, or R
                          number of nodes in given status.
    RETURNS
    -------
    report_statusX (X is 1, 2, or 3) : scipy array
                                       gives statusX subsampled just at report_times.
    '''
    report_status1 = []
    report_status2 = []
    report_status3 = []
    
    next_report_index = 0
    for index, t in enumerate(times):
        while next_report_index<len(report_times) and t>= report_times[next_report_index]:
            report_status1.append(status1[index])
            if status2 is not None:
                report_status2.append(status2[index])
            if status3 is not None:
                report_status3.append(status3[index])
            next_report_index += 1
    while next_report_index<len(report_times):
        report_status1.append(status1[-1])
        if status2 is not None:
            report_status2.append(status2[-1])
        if status3 is not None:
            report_status3.append(status3[-1])
        next_report_index += 1

    return_value = []
    return_value.append(scipy.array(report_status1))
    if status2 is not None:
        return_value.append(scipy.array(report_status2))
    if status3 is not None:
        return_value.append(scipy.array(report_status3))
    if len(return_value)==1:
        return return_value[0]
    else:
        return return_value


def get_time_shift(times, L, threshold):
    '''Identifies the first time at which L crosses a threshold.  Useful for shifting times.
    INPUTS
    ------
    times : list or scipy array (ordered)
            the times we have observations
    L : a list or scipy array
        order of L corresponds to times
    threshold : number
        a threshold value

    RETURNS
    -------
    t : number
        the first time at which L crosses a threshold.
    '''
    for index, t in enumerate(times):
        if L[index]>= threshold:
            break
    return t



##########################
#                        #
#    SIMULATION CODE     #
#                        #
##########################

def _simple_test_transmission_(u, v, p):
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm.

    test_transmission function for basic_discrete_SIR_epidemic.

    This handles the simple case where transmission occurs with probability p.

    INPUTS
    ------
    u : node
        the infected node
    v : node
        the susceptible node
    p : number between 0 and 1
        the transmission probability

    RETURNS
    -------
    True if u will infect v (given opportunity)
    False otherwise
    '''

    return random.random()<p


def discrete_SIR_epidemic(G, test_transmission=_simple_test_transmission_, args=(), initial_infecteds=None, return_node_data = False):
    #tested in test_discrete_SIR_epidemic
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Return details of epidemic curve from a discrete time simulation.
    It assumes that individuals are infected for exactly one unit of time
    and then recover with immunity.  This is defined to handle a user-defined function
         test_transmission(node1,node2,*args)
    which determines whether transmission occurs.  So elaborate rules can be created
    as desired by the user.  By default it uses _simple_test_transmission_, in which case
    args should be entered as (p,)

    INPUTS
    ----------
    G: NetworkX Graph (or some other structure which quacks like a NetworkX Graph)
        The network on which the epidemic will be simulated.
        
    test_transmission: function(u,v,*args)
        (see below for args definition)
        A function that determines whether u transmits to v.
        It returns True if transmission happens and False otherwise.
        The default will return True with probability p, where args=p

        This function can be user-defined.
        It is called like:
            test_transmission(u,v,*args)
        Note that if args is not entered, then args=(), and this call is equivalent to
            test_transmission(u,v)

    args: a list or tuple
        The arguments of test_transmission coming after the nodes.  If simply having transmission with probability p
        it should be entered as args=(p,)   [note the comma is needed to tell Python that this is really a tuple]

    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.

    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    RETURNS
    ------------
    if return_node_data is False:
       the scipy arrays: t, S, I, R
    else:
       the scipy arrays: t, S, I, R and the dicts infection_time and recovery_time

    these arrays give all the times observed and the number in each state at each time.  The dicts give times
        at which each node changed status.
    '''

    if return_node_data:
        infection_time = {}
        recovery_time = {}
    
    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]
    infecteds = initial_infecteds

    N=G.order()
    t = [0]
    S = [N-len(infecteds)]
    I = [len(infecteds)]
    R = [0]
    susceptible = defaultdict(lambda: True)  #equivalent to u.susceptible=True for all nodes.
    for u in infecteds:
        susceptible[u] = False
    while infecteds:
        new_infecteds = []
        for u in infecteds:
            for v in G.neighbors(u):
                if susceptible[v] and test_transmission(u, v, *args):
                    new_infecteds.append(v)
                    susceptible[v] = False
            if return_node_data:
                infection_time[u] = t[-1]
                recovery_time[u] = t[-1]+1
        infecteds = new_infecteds
        R.append(R[-1]+I[-1])
        I.append(len(infecteds))
        S.append(S[-1]-I[-1])
        t.append(t[-1]+1)
    if not return_node_data:
        return scipy.array(t), scipy.array(S), scipy.array(I), scipy.array(R)
    else:
        return scipy.array(t), scipy.array(S), scipy.array(I), scipy.array(R), infection_time, recovery_time



def basic_discrete_SIR_epidemic(G, p, initial_infecteds=None, return_node_data = False):
    #tested in test_basic_discrete_SIR_epidemic   
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm.

    Does a simulation of the simple case of all nodes transmitting
    with probability p independently to each neighbor and then
    recovering.

    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    RETURNS
    -------
    if return_node_data is False:
    the scipy arrays: t, S, I, R
    else:
    the scipy arrays: t, S, I, R and the dicts infection_time and recovery_time

    these scipy arrays give all the times observed and the number in each state at each time.  The dicts give times
        at which each node changed status.
'''

    return discrete_SIR_epidemic(G, _simple_test_transmission_, (p,), initial_infecteds, return_node_data)


def percolate_network(G, p):
    #tested indirectly in test_basic_discrete_SIR_epidemic   

    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Performs bond percolation on the network G with probability p

    INPUTS
    ----------
    G : NetworkX Graph
    p : number between 0 and 1
        the probability of keeping edge

    RETURNS
    -------
    H : NetworkX Graph
        A network with same nodes as G, but with each edge retained independently with probability p.
'''

    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for edge in G.edges_iter():
        if random.random()<p:
            H.add_edge(*edge)
    return H

def _edge_exists_(u, v, H):
    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Tests whether H has an edge from u to v.

    INPUTS
    ----------
    u : node
    v : node
    H : graph

    RETURNS  H.has_edge(u,v)
    -------
    True : if H has the edge
    False : if H does not have the edge
    '''
    return H.has_edge(u,v)

def percolation_based_discrete_SIR_epidemic(G, p, initial_infecteds=None, return_node_data = False):
    #tested in test_basic_discrete_SIR_epidemic   
    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    
    The simple case of all nodes transmitting with probability p independently to each neighbor
    and then recovering, but using a percolation-based approach.  See basic_discrete_SIR_epidemic
    which should produce equivalent outputs.  That algorithm will be faster than this one.  The
    value of this function is that by performing many simulations we can see that the outputs of
    the two are equivalent.  This algorithm leads to a better understanding of the theory.


    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    RETURNS
    ------------
    if return_node_data is False:
    the lists: t, S, I, R
    else:
    the lists: t, S, I, R and the dicts infection_time and recovery_time

    these lists give all the times observed and the number in each state at each time.  The dicts give times
        at which each node changed status.
'''

    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]
    H = percolate_network(G, p)
    return discrete_SIR_epidemic(H, _edge_exists_, [H], initial_infecteds, return_node_data)


def estimate_SIR_prob_size(G, p):
    #tested in test_estimate_SIR_prob_size
    r'''From figure 6.12 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    Provies an estimate of epidemic probability and size assuming a fixed transmission
    probability p.  The estimate is found by performing bond percolation and then finding
    the largest connected component in the remaining network.  This assumes that there is
    a single giant component above threshold.  It will not be an appropriate measure if the
    network is made up of several densely connected components with very weak connections between
    these components.

    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability

    RETURNS
    -------
    P, A : (numbers) estimates of the probability and proportion infected in epidemics
        (the two are equal, but each given for consistency with estimate_directed_SIR_prob_size)
    '''
    H = percolate_network(G, p)
    size = max((len(CC) for CC in nx.connected_components(H)))
    returnval = float(size)/G.order()
    return returnval, returnval


def directed_percolate_network(G, tau, gamma):
    #indirectly tested in test_estimate_SIR_prob_size
    r'''From figure 6.13 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.  This adds node and edge attributes
    which are not at present in the figure.  This option is discussed
    in the text.


    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    tau : number
        transmission rate
    gamma : number
        recovery rate

    RETURNS
    -------
    H : networkx DiGraph  (directed graph)
        a u->v edge exists in H if u would transmit to v if ever infected.
        The edge has a time attribute (time_to_infect) which gives the delay
            from infection of u until transmission occurs.
        each node has a time attribute (duration) which gives the duration of
            u,s infectious period.
    '''
    H = nx.DiGraph()
    for u in G.nodes():
        duration = random.expovariate(gamma)
        H.add_node(u, duration = duration)
        for v in G.neighbors(u):
            time_to_infect = random.expovariate(tau)
            if time_to_infect<duration:
                H.add_edge(u,v,delay=time_to_infect)
    return H
                
def _out_component_(G, source):
    '''rather than following the pseudocode in figure 6.15 of Kiss, Miller & Simon,
       this uses a built-in Networkx command.  

    finds the set of nodes (including source) which are reachable from nodes in source.

    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    source : either a node or an iterable of nodes (set, list, tuple)
        The nodes from which the infections start.  We assume no node
        will ever have a name that is an iterable of other node names.
        It will run, but may not use source user expects.

    RETURNS
    -------
    reachable_nodes : set
        the set of nodes reachable from source (including source).


    Warning: if the graph G has nodes like 1, 2, 3, and (1,2,3), then a
        source of (1,2,3) is potentially ambiguous.  It will interpret
        the source as the single node (1,2,3)
    '''
    if G.has_node(source): 
        source_nodes = {source}
    else:
        source_nodes = set(source)
        
    reachable_nodes = set()

    for node in source_nodes:
        reachable_nodes = reachable_nodes.union(set(nx.descendants(G, node)))

    return reachable_nodes

def _in_component_(G, target):
    r'''creates the _in_component_ by basically reversing _out_component_.

    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    target : a target node (or iterable of target nodes)
        The node whose infection we are interested in.

        In principle target could be an iterable, but in this case we would be finding those possible
        sources whose infection leads to infection of at least one target, not all.

    RETURNS
    -------
    source_nodes : set
        the set of nodes (including target) from which target is reachable

    Warning: if the graph G has nodes like 1, 2, 3, and (1,2,3), then a
        target of (1,2,3) is potentially ambiguous.  It will interpret
        the target as the single node (1,2,3)

    '''
    if G.has_node(target):
        #potential bug if for example G has  nodes like 1, 2, 3, and (1,2,3).  Then a target of (1,2,3) is ambiguous
        target_nodes = {target}
    else:
        target_nodes = set(target)

    source_nodes = set()
    
    for node in target_nodes:
        source_nodes = source_nodes.union(set(nx.ancestors(G, node)))

    return source_nodes


def get_infected_nodes(G, tau, gamma, initial_infecteds=None):
    r'''From figure 6.15 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm

    Finds all eventually infected nodes in a simulation, assuming that the intial infecteds are as given
    and transmission occurs with rate tau and recovery with rate gamma.  Uses a percolation-based
    approach.

    Note that the output of this algorithm is stochastic.
    
    This code has similar run-time whether an epidemic occurs or not.
    There are much faster ways to implement an algorithm giving the same output, for example by actually running an epidemic.
    
    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    tau : number
        transmission rate
    gamma : number
        recovery rate
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.

    RETURNS
    -------
    infected_nodes : set
        the set of nodes infected eventually in a simulation.
    '''    
    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]
    H = directed_percolate_network(G, tau, gamma)
    infected_nodes = _out_component_(G, initial_infecteds)
    return infected_nodes


def estimate_directed_SIR_prob_size(G, tau, gamma):
    #tested in test_estimate_SIR_prob_size
    '''
    Predicts probability and attack rate assuming continuous-time Markovian SIR disease on network G
    
    INPUTS
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    tau : number
        transmission rate
    gamma : number
        recovery rate

    RETURNS
    -------
    PE, AR  :  numbers (between 0 and 1)
        Estimates of epidemic probability and attack rate found by performing directed percolation,
        finding largest strongly connected component and finding its in/out components.
    '''
    
    H = directed_percolate_network(G, tau, gamma)
    return estimate_SIR_prob_size_from_directed_percolation(H)

def estimate_SIR_prob_size_from_directed_percolation(H):
    #indirectly tested in test_estimate_SIR_prob_size
    r'''From figure 6.17 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm

    INPUTS
    ----------
    H: directed graph (assumed to be from directed percolation on previous graph G)

    RETURNS
    -------
    PE, AR  :  numbers
        Estimates of epidemic probability and attack rate found by finding largest strongly
        connected component and finding in/out components.
    '''

    print "note --- I've renamed this from what it was called in the book when it went into production.  I plan to update the name"
    Hscc = max((nx.strongly_connected_components(H)), key = len)
    u = list(Hscc)[0]  #random.choice(Hscc)
    inC = _in_component_(H, u) #includes both H_{IN} and H_{SCC}
    outC = _out_component_(H, u) #includes both H_{OUT} and H_{SCC}
    N=float(H.order())
    PE = len(inC)/N
    AR = len(outC)/N
    return PE, AR
 
def estimate_nonMarkov_SIR_prob_size(G, xi, zeta, transmission):
    '''
    INPUTS
    ----------
    G : NetworkX Graph
        The input graph

    xi : Dict
        xi[u] gives all necessary information to determine what u's infectiousness is
    zeta : Dict
        zeta[v] gives everything needed about v's susceptibility

    transmission : user-defined function
        transmission(xi[u], zeta[v]) determines whether u transmits to v.

    RETURNS
    -------
    PE, AR  :  numbers (between 0 and 1)
        Estimates of epidemic probability and attack rate found by finding largest strongly
        connected component and finding in/out components.

    '''

    H = nonMarkov_directed_percolate_network(G, xi, zeta, transmission)
    return estimate_SIR_prob_size_from_directed_percolation(H)
        
def nonMarkov_directed_percolate_network(G, xi, zeta, transmission):
    r'''From figure 6.18 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    xi and zeta are dictionaries whatever data is needed so that xi[u]
    and zeta[v] are enough to determine the probability of a u-v transmission.

    transmissision is a user-defined function taking xi[u] and zeta[v] and
    returning True if a transmission would occur

    INPUTS
    ----------
    G : NetworkX Graph
        The input graph

    xi : Dict
        xi[u] gives all necessary information to determine what us infectiousness is
    zeta : Dict
        zeta[v] gives everything needed about vs susceptibility

    transmission : user-defined function
        transmission(xi[u], zeta[v]) determines whether u transmits to v.

    RETURNS
    -------
    networkx DiGraph (directed graph) H.  Edge u,v exists in H if it will transmit given the opportunity.
'''
    H = nx.DiGraph()
    for u in G.nodes():
        H.add_node(u)
        for v in G.neighbors(u):
            if transmission(xi[u],zeta[v]):
                H.add_edge(u,v)
    return H

def _find_trans_SIR_(Q, t, tau, source, target, status, pred_inf_time, cutoff_time):
    r'''From figure A.4 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This involves a couple additional dicts because the pseudocode is written as if the nodes
    are a separate class.  I find it easier to have a dict for that since I don't have control
    over the input graph.

    determines if a transmission from source to target will occur and if so puts into Q

    INPUTS
    ----------
    Q : A priority queue of events
    t : current time
    tau : transmission rate
    source : infected node that may transmit
    target : the possibly susceptible node that may receive a transmission
    status : a dict giving the current status of every node
    pred_inf_time : a dict giving a predicted infection time of susceptible nodes (defaults to inf)
    cutoff_time : either tmax or rec_time[source], whichever is first.

    RETURNS
    -------
    Nothing returned

    MODIFIES
    --------
    Q : Adds transmission events to Q
    pred_inf_time : updates predicted infection time of target.
    '''
    if status[target] == 'S':
        delay = random.expovariate(tau)
        inf_time = t + delay
        if inf_time< cutoff_time and inf_time < pred_inf_time[target]:
            event = Event(inf_time, 'transmit', target)
            heapq.heappush(Q, event)
            pred_inf_time[target] = inf_time

def _process_trans_SIR_(G, event, times, S, I, R, Q, status, rec_time, pred_inf_time, tau, gamma, tmax):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    INPUTS
    ----------
    G : NetworkX Graph
    event : event
         has details on node and time
    times : list
        list of times at which events have happened
    S, I, R : lists
        lists of numbers of nodes of each status at each time
    Q : heapq heap
        the queue of events
    status : dict
        dictionary giving status of each node
    rec_time : dict
        dictionary giving recovery time of each node
    pred_inf_time : dict
        dictionary giving predicted infeciton time of nodes 
    tmax : number
        max time allowed
    tau : number
        transmission rate (from node)
    gamma : number
        recovery rate (of node)

    RETURNS
    -------
    nothing returned

    MODIFIES
    --------
    status : updates status of newly infected node
    rec_time : adds recovery time for node
    times : appends time of event
    S : appends new S (reduced by 1 from last)
    I : appends new I (increased by 1)
    R : appends new R (same as last)
    Q : adds recovery and transmission events for newly infected node.
    pred_inf_time : updated for nodes that will receive transmission
                    update happens in _find_trans_SIR_

    Entry requirement:
    -------
    Only enter this if the node is SUSCEPTIBLE and is becoming INFECTED.
    '''
    node = event.node  #The node that is receiving transmission
    time = event.time  #current time
    assert(status[node]=='S')
    status[node] = 'I' 
    times.append(time)
    S.append(S[-1]-1) #one less susceptible
    I.append(I[-1]+1) #one more infected
    R.append(R[-1])   #no change to recovered
    rec_time[node] = time + random.expovariate(gamma)
    if rec_time[node] < tmax:
        newevent = Event(rec_time[node], 'recover', node)
        heapq.heappush(Q,newevent)
    cutoff_time = min(tmax, rec_time[node])
    for v in G.neighbors(node):
        _find_trans_SIR_(Q, time, tau, node, v, status, pred_inf_time, cutoff_time)
    
def _process_recovery_SIR_(event, times, S, I, R, status):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    INPUTS
    ----------
    event : event
         has details on node and time
    times : list
        list of times at which events have happened
    S, I, R : lists
        lists of numbers of nodes of each status at each time
    status : dict
        dictionary giving status of each node


    RETURNS
    ----------
    Nothing

    MODIFIES
    ----------
    status : updates status of newly recovered node
    times : appends time of event
    S : appends new S (same as last)
    I : appends new I (decreased by 1)
    R : appends new R (increased by 1)
    '''
    node = event.node
    time = event.time
    times.append(time)
    S.append(S[-1])   #no change to number susceptible
    I.append(I[-1]-1) #one less infected
    R.append(R[-1]+1) #one more recovered
    status[node] = 'R'
    
    
def fast_SIR(G, tau, gamma, initial_infecteds = None, tmax=float('Inf'), return_node_data = False):
    #tested in test_SIR_dynamics
    r'''From figure A.2 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    fast SIR simulation assuming exponentially distributed infection and recovery times
    

    INPUTS
    ------
    G : NetworkX Graph
       The underlying network
    tau : number
       transmission rate per edge
    gamma : number
       recovery rate per node
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    tmax : number
       maximum time after which simulation will stop.
       default float('Inf') to set to infinity.  Okay for SIR, not for SIS.
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    RETURNS
    -------
    times, S, I, R : each a scipy array
         giving times and number in each status for corresponding time

    OR if return_node_data=True:
    times, S, I, R, infection_time, recovery_time
         first four are scipy arrays as above.  New objects are dicts
         with entries just for those nodes that were infected ever
         infection_time[node] is time of infection
         recovery_time[node] is time of recovery
    
    
    '''
    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]


    Q = []#an empty heap
        
    times, S, I, R= ([0], [G.order()], [0], [0])

    status = defaultdict(lambda: 'S') #node status defaults to 'S'
    rec_time = defaultdict(lambda: -1) #node recovery time defaults to -1
    pred_inf_time = defaultdict(lambda: float('Inf')) #infection time defaults to \infty  --- this could be set to tmax, probably with a slight improvement to performance.
    #Note that if node becomes infected, pred_inf_time is actually inf_time
    #and similarly for rec_time rec_time is correct.  
    for u in initial_infecteds:
        newevent = Event(0, 'transmit', u)
        pred_inf_time[u] = 0
        heapq.heappush(Q,newevent)

    while Q:
        event = heapq.heappop(Q)
        if event.action == 'transmit':
            if status[event.node] == 'S': 
                _process_trans_SIR_(G, event, times, S, I, R, Q, status, rec_time, pred_inf_time,  tau, gamma, tmax)
        else:
            _process_recovery_SIR_(event, times, S, I, R, status)

    #the initial infections were treated as ordinary infection events at time 0.  So
    #each initial infection added an entry at time 0 to lists.  We'd like to get rid
    #these excess events.
    times = times[len(initial_infecteds):]
    S=S[len(initial_infecteds):]
    I=I[len(initial_infecteds):]
    R=R[len(initial_infecteds):]
    if not return_node_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R) # ignoring initial value entered.
    else:
        #strip pred_inf_time and rec_time down to just the values for nodes that became infected
        infection_time = {node:time for (node,time) in pred_inf_time.iteritems() if status[node]!='S'}
        recovery_time = {node:time for (node,time) in rec_time.iteritems() if status[node] !='S'}
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R), infection_time, recovery_time

def fast_nonMarkov_SIR(G, process_trans = _process_trans_SIR_, args = (), initial_infecteds = None, tmax = float('Inf'), return_node_data = False, Q=None):
    r'''A modification of the algorithm in figure A.2 of Kiss, Miller, & Simon to allow for user-defined rules
    governing time of transmission.  Please cite the book if using ethis algorithm.

    This is useful if the transmission rule is non-Markovian in time, or for more elaborate models.  For example if there is a mass action style transmission
    this can be incorporated into the process_trans command defined by user.

    INPUTS
    ------
    G : Networkx Graph
    process_trans : a function that handles a transmission event.
                    Called by process_trans(G, event, times, S, I, R, Q, status, rec_time, pred_inf_time, tmax, *args)
                    must update :   status, rec_time, times, S, I, R,
                    must also update : Q, pred_inf_time.
                    In updating these last two, it calculates the recovery time, and adds the event to Q.  It then calculates
                       predicted times of transmission to neighbors.  If before current earliest prediction, it will add
                       appropriate transmission event to Q and update this prediction.
    args: The final arguments going into process_trans.  If there is some reason to collect data about node that is only calculated when transmission occurs
          it can modify a dict or something similar that is passed as an argument.
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    tmax : number
        stop time
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    
         If Q is defined:
             Then initial_infecteds consists of those nodes infected **PRIOR** to t=0.  Those at t=0 will be handled by being in Q already.
    tmax : (default infinity)
        final time
    return_node_data : boolean (default False)
    Q : If user wants to predefine some events, this can be done.  This can be input as a heap or as a list (it will be heapified and modified)
        User should understand the Event class and use it.  Currently there is no guarantee this is properly supported.
        So much so that right now I'm going to force the user to edit the source code before trying it.  I am reasonably confident it will work.

        When Q is input, initial_infecteds should be the of nodes in I class **PRIOR** to t=0, and the events in Q must have all of their recoveries.  The best way to handle nodes that should be already recovered is to put them in initial_infecteds and give them a recovery event at t=0
    '''

    if process_trans == _process_trans_SIR_:
        try:
            tau, gamma = map(float, args)  #better be able to make args be two floats.
        except:
            raise EoNError("if using default for fast_nonMarkov_SIR, then args should be (tau,gamma)\nConsider just using fast_SIR.")

    status = defaultdict(lambda: 'S') #node status defaults to 'S'
    rec_time = defaultdict(lambda: -1) #node recovery time defaults to -1
    pred_inf_time = defaultdict(lambda: float('Inf')) #infection time defaults to \infty  --- this could be set to tmax, probably with a slight improvement to performance.
    
    if Q is None:
        Q = []
        if initial_infecteds is None:
            initial_infecteds=[random.choice(G.nodes())]
        elif G.has_node(initial_infecteds):
            initial_infecteds=[initial_infecteds]
        for u in initial_infecteds:
            newevent = Event(0, 'transmit', u)
            pred_inf_time[u] = 0
            heapq.heappush(Q,newevent)
        times, S, I, R= ([0], [G.order()], [0], [0])  

    else:
        raise EoNError("inputting Q is not currently tested.\n Email joel.c.miller.research@gmail.com for help.\n I believe this code will work, but you will need to delete this message.")

        if initial_infecteds is None:
            initial_infecteds = []  #assuming that input Q has this taken care of 
        for node in initial_infecteds:
            status[node] = 'I'
            pred_inf_time[node] = -1
        for event in Q:
            if event.action == 'transmit' and event.time<pred_inf_time[event.node]:
                pred_inf_time[event.node] = event.time
            elif event.action == 'recover':
                rec_time[event.node] = event.time
        heapq.heapify(Q)            
        times, S, I, R= ([0], [G.order()-len(initial_infecteds)], [0], [0])  
        
    
    #Note that when finally infected, pred_inf_time is correct
    #and rec_time is correct.  So if return_node_data is true, these are correct

    while Q:
        event = heapq.heappop(Q)
        if event.action == 'transmit':
            if status[event.node] == 'S': 
                process_trans(G, event, times, S, I, R, Q, status, rec_time, pred_inf_time, tmax, *args)
        else:
            _process_recovery_SIR_(event, times, S, I, R, status)

    if not return_node_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R) 
    else:
        #strip pred_inf_time and rec_time down to just the values for nodes that became infected
        infection_time = {node:time for (node,time) in pred_inf_time.iteritems() if status[node]!='S'}
        recovery_time = {node:time for (node,time) in rec_time.iteritems() if status[node] !='S'}
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R), infection_time, recovery_time


def _process_trans_SIS_(G, event, tau, gamma, times, S, I, Q, status, rec_time, tmax):
    r'''From figure A.5 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    INPUTS
    ----------
    G : NetworkX Graph
    event: event
        has details on node and time
    tau : number
        transmission rate (from node)
    gamma : number
        recovery rate (of node)
    times : list
        list of times at which events have happened
    S, I: lists
        lists of numbers of nodes of each status at each time
    Q : heapq heap
        the queue of events
    status : dict
        dictionary giving status of each node
    rec_time : dict
        dictionary giving recovery time of each node
    pred_inf_time : dict
        dictionary giving predicted infeciton time of nodes 
    tmax : max time allowed

    RETURNS
    -------
    nothing returned

    MODIFIES
    --------
    status : updates status of newly infected node
    rec_time : adds recovery time for node
    times : appends time of event
    S : appends new S (reduced by 1 from last)
    I : appends new I (increased by 1)
    Q : adds recovery and transmission events for newly infected node.

    Entry requirement:
    -------
    Only enter this if the node is SUSCEPTIBLE and is becoming INFECTED.
    '''
    node = event.node  #The node that is receiving transmission
    time = event.time  #current time

    assert(status[node]=='S')
    status[node] = 'I'
    I.append(I[-1]+1) #one more infected
    S.append(S[-1]-1) #one less susceptible
    times.append(time)
    rec_time[node] = time + random.expovariate(gamma)
    if rec_time[node] < tmax:
        newevent = Event(rec_time[node], 'recover', node)
        heapq.heappush(Q, newevent)
    for v in G.neighbors(node):
        _find_next_trans_SIS_(Q, time, tau, node, v, status, rec_time, tmax)

def _find_next_trans_SIS_(Q, time, tau, source, target, status, rec_time, tmax):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    determines if a transmission from source to target will occur and if so puts into Q

    INPUTS
    ----------
    Q : A priority queue of events
    t : current time
    tau : transmission rate
    source : infected node that may transmit
    target : the possibly susceptible node that may receive a transmission
    status : a dict giving the current status of every node
    rec_time : a dict giving the recovery time of every node that has been infected.  
    tmax : max time for simulation.

    RETURNS
    -------
    nothing returned

    MODIFIES
    --------
    Q : if a transmission time is potentially valid, add the first event.
        when this transmission occurs later we will consider adding another event.
        note that the event includes the source, so we can later check if same source
        will transmit again.

    Entry requirement:
    -------
    Only enter this if the source node is INFECTED.

    '''

    assert(status[source]=='I')
    if rec_time[target]<rec_time[source]: #if target is susceptible, then rec_time[target]<time
        transmission_time = max(time, rec_time[target])+random.expovariate(tau)
        if transmission_time < rec_time[source] and transmission_time<tmax:
            newEvent = Event(transmission_time, 'transmit', target, source)#{'node': target, 'action': 'transmit', 'source':source}
            heapq.heappush(Q, newEvent)
                #        if inf_time< rec_time[source]:#and inf_time < pred_inf_time[target] and inf_time<tmax:
                #            transmission_time = max(t, rec_time[target])+ran
                #            event = {'node': target, 'action': 'transmit'}
                #            heapq.heappush(Q, (inf_time, event))
                #            pred_inf_time[target] = inf_time

 
def _process_recovery_SIS_(event, times, S, I, status):
    r'''From figure A.4 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    '''
    node = event.node
    time = event.time

    times.append(time)
    S.append(S[-1]+1)   #no change to number susceptible
    I.append(I[-1]-1) #one less infected
    status[node] = 'S'




def fast_SIS(G, tau, gamma, initial_infecteds=None, tmax=100, return_node_data = False):
    r'''From figure A.5 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    INPUTS
    -------------
    G : NetworkX Graph
       The underlying network
    tau : number
       transmission rate per edge
    gamma : number
       recovery rate per node
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    tmax : number
        stop time
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time
    '''

    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]

    times = [0]
    S = [G.order()]
    I = [0]
    Q = []
    status = defaultdict(lambda: 'S') #node status defaults to 'S'
    rec_time = defaultdict(lambda: -1) #node recovery time defaults to -1
    rec_time['initial_condition'] = 0

    infection_times = defaultdict(lambda: []) #defaults to an empty list for each node
    recovery_times = defaultdict(lambda: [])

    for u in initial_infecteds:
        newevent = Event(0, 'transmit', u, source='initial_condition')
        Q.append(newevent)#okay to append rather than heappush since all are at same time
    while Q:
        event = heapq.heappop(Q)
        assert(event.time<= tmax)
        node = event.node
        time = event.time
        if event.action == 'transmit':
            source = event.source
            if status[node] == 'S':
                _process_trans_SIS_(G, event, tau, gamma, times, S, I, Q, status, rec_time, tmax)
                infection_times[node].append(time)
            if source != 'initial_condition':
                _find_next_trans_SIS_(Q, time, tau, source, node, status, rec_time, tmax)
        else:
            _process_recovery_SIS_(event, times, S, I, status)
            recovery_times[node].append(time)
    if not return_node_data:
        return scipy.array(times), scipy.array(S), scipy.array(I)
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), infection_times, recovery_times

def _Gillespie_Initialize_(G, initial_infecteds, infection_times, return_node_data, SIR = True):
    '''Initializes the network'''
    times = [0]
    S = [G.order()-len(initial_infecteds)]
    I = [len(initial_infecteds)]
    R = [0]
    status = defaultdict(lambda:'S') #by default all are susceptible
    infected = list(initial_infecteds)
    infected_neighbor_count = defaultdict(lambda:0)#
    risk_group = defaultdict(lambda:_ListDict_()) 
    for node in initial_infecteds:
        status[node]='I'
    for node in initial_infecteds:
        for neighbor in G.neighbors(node):
            if status[neighbor]=='S':
                infected_neighbor_count[neighbor] += 1
                if infected_neighbor_count[neighbor]>1:
                    risk_group[infected_neighbor_count[neighbor]-1].remove(neighbor)
                risk_group[infected_neighbor_count[neighbor]].add(neighbor)
    if return_node_data:
        for node in initial_infecteds:
            infection_times[node].append(0)
    if SIR:
        return times, S, I, R, status, infected, infected_neighbor_count, risk_group
    else:
        return times, S, I, status, infected, infected_neighbor_count, risk_group
    
def _Gillespie_Infect_(G, S, I, R, times, infected, current_time, infected_neighbor_count, risk_group, status, infection_times, return_node_data, SIR=True):
    '''
    Chooses the node to infect.

    First chooses which risk_group the node is in.  Then choose the node.

    An alternative which may be cleaner (not sure about faster?) is to create a data type which tracks maximum
    risk.  Then choose a random node and do a rejection sampling step.  If rejected, then select again.  May need a very careful rejection sampling to account for repeated selection.  We could probably define this as a method of the class.
    '''
    r = random.random()*sum(n*len(risk_group[n]) for n in risk_group.keys())
    for n in risk_group.keys():
        r-= n*len(risk_group[n])
        if r<0:
            break
    #we've got n now

    recipient = risk_group[n].choose_random()
                                                    #OLD VERSION choose random element from dict
                                                    #based on http://stackoverflow.com/a/24949742/2966723
                                                    #question by http://stackoverflow.com/users/2237265/jamyn
                                                    #answer by http://stackoverflow.com/users/642757/scott-ritchie
                                                    #
                                                    #CURRENT VERSION
                                                    #http://stackoverflow.com/a/15993515/2966723
    assert(status[recipient]=='S')
    risk_group[n].remove(recipient)
    infected.append(recipient)
    infection_times[recipient].append(current_time)
    status[recipient]='I'
    S.append(S[-1]-1)
    I.append(I[-1]+1)
    times.append(current_time)
    if SIR:
        R.append(R[-1])

    for neighbor in G.neighbors(recipient):
        if status[neighbor]=='S':
            if infected_neighbor_count[neighbor]>0:
                risk_group[infected_neighbor_count[neighbor]].remove(neighbor)
            infected_neighbor_count[neighbor]+=1
            risk_group[infected_neighbor_count[neighbor]].add(neighbor)
    
            

def _Gillespie_Recover_SIR_(G, S, I, R, times, infected, current_time, status, infected_neighbor_count, risk_group, recovery_times, return_node_data):
    r''' Changes: S, I, R, infected, times'''
    assert(I[-1]==len(infected))
    index = random.randint(0,I[-1]-1)
    infected[index], infected[-1] = infected[-1], infected[index] #http://stackoverflow.com/a/14088129/2966723
    recovering_node = infected.pop()

    I.append(I[-1]-1)
    status[recovering_node]='R'
    S.append(S[-1])
    R.append(R[-1]+1)
    times.append(current_time)
    for neighbor in G.neighbors(recovering_node):
        if status[neighbor] == 'S': #neighbor susceptible, its risk just got smaller
            risk_group[infected_neighbor_count[neighbor]].remove(neighbor)
            infected_neighbor_count[neighbor] -= 1
            if infected_neighbor_count[neighbor]>0:
                risk_group[infected_neighbor_count[neighbor]].add(neighbor)            
    if return_node_data:
        recovery_times[recovering_node].append(current_time)

def _Gillespie_Recover_SIS_(G, S, I, times, infected, current_time, status, infected_neighbor_count, risk_group, recovery_times, return_node_data):
    r''' x'''
    assert(I[-1]==len(infected))
    index = random.randint(0,I[-1]-1)
    infected[index], infected[-1] = infected[-1], infected[index] #http://stackoverflow.com/a/14088129/2966723
    recovering_node = infected.pop()

    I.append(I[-1]-1)
    status[recovering_node]='S'
    S.append(S[-1]+1)
    times.append(current_time)
    infected_neighbor_count[recovering_node] = 0
    for neighbor in G.neighbors(recovering_node):  #there is probably a good way to count the number of infected neighbors
        if status[neighbor] == 'I':
            infected_neighbor_count[recovering_node] += 1
        else: #neighbor susceptible, its risk just got smaller
            risk_group[infected_neighbor_count[neighbor]].remove(neighbor)
            infected_neighbor_count[neighbor] -= 1
            if infected_neighbor_count[neighbor]>0:
                risk_group[infected_neighbor_count[neighbor]].add(neighbor)
    if infected_neighbor_count[recovering_node]>0:
        risk_group[infected_neighbor_count[recovering_node]].add(recovering_node)
    if return_node_data:
        recovery_times[recovering_node].append(current_time)

def Gillespie_SIR(G, tau, gamma, initial_infecteds=None, tmax=float('Inf'), return_node_data = False):
    #tested in test_SIR_dynamics
    r'''
    Performs a Gillespie-based SIR simulation.

    Assumes that the network is unweighted.  Thus the risks are
    quantized: equal to tau times the number of infected neighbors of
    a node.
    
    This would not be as good if the edges were weighted, but we could put the at_risk nodes into bands.
    At present the network is treated as unweighted, so the possible risk rates are discretized.
    
    The event-driven simulation is almost certainly faster in all cases, and the benefit would increase
    if the network were weighted.  I think the coding would also be easier.

    INPUTS
    ------
    G : NetworkX Graph
       The underlying network
    tau : number
       transmission rate per edge
    gamma : number
       recovery rate per node
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    tmax : number
        stop time
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time

    '''

    infection_times = defaultdict(lambda: []) #defaults to an empty list for each node
    recovery_times = defaultdict(lambda: [])

    tau = float(tau)  #just to avoid integer division problems.
    gamma = float(gamma)
    
    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]

    times, S, I, R, status, infected, infected_neighbor_count, risk_group = _Gillespie_Initialize_(G, initial_infecteds, infection_times, return_node_data)

    total_trans_rate = tau*sum(n*len(risk_group[n]) for n in risk_group.keys())
    total_rec_rate = gamma*len(infected)
    total_rate = total_rec_rate + total_trans_rate
    next_time = times[-1] + random.expovariate(total_rate)
    while next_time<tmax and infected:
        r = random.random()*total_rate
        if r<total_rec_rate:
            #a recovery occurs
            _Gillespie_Recover_SIR_(G, S, I, R, times, infected, next_time, status, infected_neighbor_count, risk_group, recovery_times, return_node_data)
            total_rec_rate = gamma*I[-1]
        else:
            #an infection occurs
            _Gillespie_Infect_(G, S, I, R, times, infected, next_time, infected_neighbor_count, risk_group, status, infection_times, return_node_data, SIR=True)
        total_trans_rate = tau*sum(n*len(risk_group[n]) for n in risk_group.keys())
        total_rate = total_rec_rate + total_trans_rate
        if total_rate >0:  
            next_time += random.expovariate(total_rate)

    if not return_node_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R)
    else:
        #need to change data type of infection_times and recovery_times
        infection_time = {node: L[0] for node, L in infection_times.items()}
        recovery_time = {node: L[0] for node, L in recovery_times.items()}
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R), infection_time, recovery_time



def Gillespie_SIS(G, tau, gamma, initial_infecteds=None, tmax=100, return_node_data = False):
    r'''
    This could be made more efficient if we divide the at_risk nodes into groups based on how at risk they are. 
    This would not be as good if the edges were weighted, but we could put the at_risk nodes into bands.

    Warning: self-edges will cause this to die.  You can remove self-edges by G.remove_edges_from(G.selfloop_edges())

    INPUTS
    ------
    G : NetworkX Graph
       The underlying network
    tau : number
       transmission rate per edge
    gamma : number
       recovery rate per node
    initial_infecteds: node or iterable of nodes
       if a single node, then this node is initially infected
       if an iterable, then whole set is initially infected
       if None, then a randomly chosen node is initially infected.
    tmax : number
        stop time
    return_node_data: boolean
        Tells whether the infection and recovery times of each individual node
        should be returned.  It is returned in the form of two dicts, infection_time and recovery_time
        infection_time[node] is the time of infection and recovery_time[node] is the recovery time
    '''

    infection_times = defaultdict(lambda: []) #defaults to an empty list for each node
    recovery_times = defaultdict(lambda: [])

    tau = float(tau)  #just to avoid integer division problems.
    gamma = float(gamma)
    
    if initial_infecteds is None:
        initial_infecteds=[random.choice(G.nodes())]
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]

    times, S, I, status, infected, infected_neighbor_count, risk_group = _Gillespie_Initialize_(G, initial_infecteds,  infection_times,  return_node_data, SIR=False)
    #note that at this point times, S, and I must all be lists since we will be appending to them

    total_trans_rate = tau*sum(n*len(risk_group[n]) for n in risk_group.keys())
    total_rec_rate = gamma*len(infected)
    total_rate = total_rec_rate + total_trans_rate
    next_time = times[-1] + random.expovariate(total_rate)
    while next_time<tmax and infected:
        r = random.random()*total_rate
        if r<total_rec_rate:
            #a recovery occurs
            _Gillespie_Recover_SIS_(G, S, I, times, infected, next_time, status, infected_neighbor_count, risk_group, recovery_times, return_node_data)
            total_rec_rate = gamma*I[-1]
        else:
            #an infection occurs
            _Gillespie_Infect_(G, S, I, [], times, infected, next_time, infected_neighbor_count, risk_group, status, infection_times, return_node_data, SIR=False)
            #updates variables as needed and calculates new max_trans_rate
        total_trans_rate = tau*sum(n*len(risk_group[n]) for n in risk_group.keys())
        total_rate = total_rec_rate + total_trans_rate
        if total_rate>0:
            next_time += random.expovariate(total_rate)
        else:  #occurs if everyone recovered
            next_time = float('Inf')
        #        print next_time, I[-1]

    if not return_node_data:
        return scipy.array(times), scipy.array(S), scipy.array(I)
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), infection_times, recovery_times




##################
#                #
#    ODE CODE    #
#                #
##################
'''All system numbers below are based on the current draft of the book.  They are subject to change'''

'''Code will be ordered so that when SIS and SIR versions both exist, the corresponding functions are adjacent to each other.'''



########      INDIVIDUAL BASED code -
########  given node and who its neighbors are, we track probability of
########  having given status based on probabilities of neighbors.  Assumes
########  independence.
    
def _dSIS_individual_based_(Y, t, G, nodelist, trans_rate, rec_rate):
    N = len(nodelist)
    dY = scipy.zeros(N)
    for index, (node, Yi) in enumerate(zip(nodelist,Y)):
        #This would probably be faster if it were done as a
        #matrix multiplication, but then I'd need to have
        #the matrix Gij explicitly included.  Perhaps that
        #would be better.  Networkx has something to create
        #numpy sparse matrices.  Perhaps that works?
        #No plan to do premature optimization.  Let's get it
        #working and then see if it's slow.
        dY[index] = sum(trans_rate(node,nbr)*(1-Y[node])*Y[nbr] for nbr in G.neighbors(node)) - rec_rate(node)*Yi
    return dY

def _dSIR_individual_based_(V, t, G, nodelist, index_of_node, trans_rate, rec_rate):
    '''    <\dot{X}_i> = - tau sum_j g_{ij} <Xi><Yj>
    <\dot{Y}_i> = tau sum_j g_{ij} <Xi><Yj> - gamma_i <Y_i>
    Z_i = 1-X_i-Y_i
    '''
    N = len(nodelist)
    X = V[:N]
    Y = V[N:]
    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    for index, (node, Xi, Yi) in enumerate(zip(nodelist,X, Y)):
        #This would probably be faster if it were done as a
        #matrix multiplication, but then I'd need to have
        #the matrix Gij explicitly included.  Perhaps that
        #would be better.  Networkx has something to create
        #numpy sparse matrices.  Perhaps that works?
        #No plan to do premature optimization.  Let's get it
        #working and then see if it's slow.
        
        dX[index] = -Xi*sum(trans_rate(node,nbr)*Y[index_of_node[nbr]] for nbr in G.neighbors(node))
        dY[index] =  -dX[index] - rec_rate(node)*Yi
    dV = scipy.concatenate((dX,dY), axis=0)
    return scipy.array(dV)

def SIS_individual_based(G, nodelist, Y0, tau, gamma=None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    #tested in test_SIS_individual_based
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology

    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    INPUTS
    -------
    G : Networkx graph
    
    Y0 : scipy array
         the array of initial infection probabilities

    nodelist : list
         list of nodes in G in the same order as in Y0

    tau : number
          transmission rate of disease

    gamma : number      (default None)
            global recovery rate  (incompatible with recovery_label!=None)

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    edge_label : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}

    recovery_label : string       (default None)
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error

    return_full_data       (default False)
            If True, returns times, Ss, Is
            if False, returns times, S, I

    RETURNS
    ---------

    '''
    if gamma is None:
        if recovery_label is None:
            EoNError("need gamma or recovery_label defined.")
        else:
            rec_rate = lambda x : G.node[x][recovery_label]
    else:
        if recovery_label is not None:
            EoNError("only one of gamma and recovery_label can be defined.")
        rec_rate = lambda x : gamma
    if edge_label is not None:
        trans_rate = lambda x, y: tau*G.edge[x][y][edge_label]
    else:
        trans_rate = lambda x, y: tau

    times = scipy.linspace(tmin, tmax, tcount)
    Y = integrate.odeint(_dSIS_individual_based_, Y0, times, args = (G, nodelist, trans_rate, rec_rate))
    Is = Y.T
    Ss = scipy.ones(len(Is))[:,None]-Is 
    
    if return_full_data:
        return times, Ss, Is
    else:
        return times, sum(Ss), sum(Is)

def SIR_individual_based(G, nodelist, X0, Y0, tau, gamma=None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    '''Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:

    INPUTS
    -------
    G : Networkx graph

    X0 : scipy array
         the array of initial susceptibility probabilities
    Y0 : scipy array
         the array of initial infection probabilities

    nodelist : list
         list of nodes in G in the same order as in X0 and Y0

    tau : number
          transmission rate of disease

    gamma : number      (default None)
            global recovery rate  (incompatible with recovery_label!=None)

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    edge_label : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}

    recovery_label : string       (default None)
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error

    return_full_data       (default False)
            If True, returns times, S, I, R, Ss, Is, Rs
            if False, returns times, S, I, R

    RETURNS
    ---------

    '''

    if gamma is None:
        if recovery_label is None:
            EoNError("need gamma or recovery_label defined.")
        else:
            rec_rate = lambda x : G.node[x][recovery_label]
    else:
        if recovery_label is not None:
            EoNError("only one of gamma and recovery_label can be defined.")
        rec_rate = lambda x : gamma
    if edge_label is not None:
        trans_rate = lambda x, y: tau*G.edge[x][y][edge_label]
    else:
        trans_rate = lambda x, y: tau

    index_of_node = {}
    for i, node in enumerate(nodelist):
        index_of_node[node] = i

    N = len(X0)
    times = scipy.linspace(tmin, tmax, tcount)
    V0 = scipy.concatenate((X0,Y0), axis=0)
    V = integrate.odeint(_dSIR_individual_based_, V0, times, args = (G, nodelist, index_of_node, trans_rate, rec_rate))
    Ss = V.T[:N]
    S = Ss.sum(axis=0)
    Is = V.T[N:]
    I = Is.sum(axis=0)
    Rs = scipy.ones(N)[:, None] - Ss - Is
    R = Rs.sum(axis=0)
    if return_full_data:
        return times, S, I, R, Ss, Is, Rs
    else:
        return times, S, I, R


def SIS_individual_based_pure_IC(G, index_nodes, tau, gamma=None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIS_individual_based is that this one assumes a "pure initial condition", that is, we know exactly what the statuses of the nodes are at the initial time.  
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    INPUTS
    -------
    G : Networkx graph
    index_nodes : list or set
      the set of nodes initially infected
    nodelist : list
         list of nodes in G in the same order as in Y0
    tau : number
          transmission rate of disease
    gamma : number      (default None)
            global recovery rate  (incompatible with recovery_label!=None)

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    edge_label : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}

    recovery_label : string       (default None)
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error

    return_full_data : boolean      (default False)


    RETURNS:
    --------
            if return_full_data is True,
                returns times, Ss, Is, nodelist
                    the reason it also returns nodelist is that the order of nodelist
                    coming from G.nodes() isn't guaranteed to be the same in future calls
            if return_full_data is False,
                returns times, S, I
    '''
    nodelist = G.nodes()
    #make Y0[u] be 1 if infected 0 if not
    Y0 = scipy.array([1 if u in index_nodes else 0 for u in nodelist])
    if return_full_data:
        times, Ss, Is = SIS_individual_based(G, nodelist, Y0, tau, gamma, tmin, tmax, tcount, edge_label, recovery_label, return_full_data)
        return times, Ss, Is, nodelist
    else:
        SIS_individual_based(G, nodelist, Y0, tau, gamma, tmin, tmax, tcount, edge_label, recovery_label, return_full_data)



def SIR_individual_based_pure_IC(G, index_nodes, initial_susceptible, tau, gamma=None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    '''Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIR_individual_based is that this one assumes a "pure initial condition", that is, we know exactly what the statuses of the nodes are at the initial time.  
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    INPUTS
    -------
    G : Networkx graph
    index_nodes : list or set
      the set of nodes initially infected
    initial_susceptible : list or set
      initially susceptible nodes
    nodelist : list
         list of nodes in G in the same order as in Y0
    tau : number
          transmission rate of disease
    gamma : number      (default None)
            global recovery rate  (incompatible with recovery_label!=None)
    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    edge_label : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}

    recovery_label : string       (default None)
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error

    return_full_data : boolean      (default False)

    RETURNS:
            if return_full_data is True,
                returns times, S, I, R, Ss, Is, Rs, nodelist
                    the reason it also returns nodelist is that the order of nodelist
                    coming from G.nodes() isn't guaranteed to be the same in future calls
            if return_full_data is False,
                returns times, S, I, R
    
    '''
    nodelist = G.nodes()
    N = len(nodelist)
    #make Y0[u] be 1 if infected 0 if not
    Y0 = scipy.array([1 if u in index_nodes else 0 for u in nodelist])
    X0 = scipy.array([1 if u in initially_susceptible else 0 for u in nodelist])
    
    if return_full_data:
        times, S, I, R, Ss, Is, Rs = SIR_individual_based(G, nodelist, X0, Y0, tau, gamma, tmin, tmax, tcount, edge_label, recovery_label, True)
        return times, S, I, R, Ss, Is, Rs, nodelist
    else:
        return SIR_individual_based(G, nodelist, Y0, tau, gamma, tmin, tmax, tcount, edge_label, recovery_label, False)











########   PAIR BASED

def _dSIS_pair_based_(V, t, G, nodelist, index_of_node, trans_rate, rec_rate):
    '''
    <\dot{Y}_i> = tau \sum_j g_{ij} <XiYj>  -  gamma_i <Yi>
    <\dot{XY}_ij> = tau sum_{k \neq i} g_{jk} <XiXj><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XiYj>/<Xi>
                   - tau g_{ij}<XiYj> - gamma_j <XiYj>  + ***gamma_i <YiYj>***
    <\dot{XX}_ij> = - tau sum_{k\neq i} g_{jk} <XiXj><XjYk>/<Xj>
                    - tau sum_{k neq j} g_{ik} <YkXi><XiXj>/<Xi>
                    **** + \gamma_i <YiXj> + gamma_j <XiYj>****

    <Xi>=1-<Yi>
    <YiYj> = 1 - <XiXj> - <XiYj> - <XjYi>
    <YiXj> = <XjYi>

    (Starred terms differ from SIR)
    
    The equations as coded involve all pairs rather than just the
    pairs that are in edges.  Those that are not part of an edge are
    set to zero and their derivatives are zero.  So the code could run
    faster if we took these out of the calculation.  That is a
    potential future improvement.

    I think for most cases this is a small contribution.
    All derivatives are initialized to 0, and then the loop only makes changes
    for those terms where an edge exists.
    '''
    N=G.order()
    Y = V[0:N] #infecteds
    X = 1-Y    #susceptibles
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) #there are places where we divide by X[i] which may = 0.
                                                         #In those cases the numerator is (very) 0, so it's easier
                                                         #to set this up as mult by inverse with a dummy value when
                                                         #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[N: N+N**2]
    XX = V[N+N**2:]

    print X.shape, Y.shape, XY.shape, XX.shape, N
    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T #not really needed, but helps keep consistent with equations as written.

    YY = 1 - XY-XX-YX

    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, but I think this sequence of for loops is easier to read, or at least understand.
    #I expect this isn't the bottleneck.  Will avoid (premature) optimization for now.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dY[i] += trans_rate(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate(u,v)+rec_rate(v))*XY[i,j] + rec_rate(u)*YY[i,j]
            dXX[i,j] +=  rec_rate(u)*YX[i,j] + rec_rate(v)*XY[i,j]
            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(u):
                if w == v: #skip these
                    continue
                #so w != v. 
                k= index_of_node[w]

                dXY[i,j] += trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  -  trans_rate(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j] -  trans_rate(u,w) * YX[k,i] * XX[i,j]*Xinv[i]


                #following gives an implementation accounting for triangles.  Have removed that.
                #if G.has_edge(w,v): #it's a triangle
                #    dXY[i,j] += trans_rate(v,w) * XX[i,j] * XY[j,k] * XY[i,k]*(Xinv[i]*Xinv[j]*Yinv[k])  -  trans_rate(u,w) * YX[k,i] * XY[i,j] * YY[j,k]*(Xinv[i]*Yinv[j]*Yinv[k])
                #    dXX[i,j] += -trans_rate(v,w) * XX[i,j] * XY[j,k] * XY[i,k]*(Xinv[i]*Xinv[j]*Yinv[k]) -  trans_rate(u,w) * YX[k,i] * XX[i,j] * YX[k,j]*(Xinv[i]*Xinv[j]*Yinv[k])
                #else: #it's not a triangle
                #    dXY[i,j] += trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  -  trans_rate(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                #    dXX[i,j] += -trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j] -  trans_rate(u,w) * YX[k,i] * XX[i,j]*Xinv[i]  

    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dY[:,None], dXY, dXX), axis=0).T[0]
    print t, sum(dV)
    return dV

def _dSIR_pair_based_(V, t, G, nodelist, index_of_node, trans_rate, rec_rate):
    '''
    <\dot{X}_i> = -tau sum_j g_{ij} <XiYj>
    <\dot{Y}_i> = tau \sum_j g_{ij} <XiYj>  -  gamma_i <Y_i>
    <\dot{XY}_ij> = tau sum_{k \neq i} g_{jk} <XiXj><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XiYj>/<Xi>
                   - tau g_{ij}<XiYj> - gamma_j <XiYj> 
    <\dot{XX}_ij> = -tau sum_{k\neq j} gik <YkXi><XiXj>/<Xi>
                    -tau sum_{k neq i} gjk <XiXj><XjYk>/<Xj>
    <>
    The equations as coded involve all pairs rather than just the
    pairs that are in edges.  Those that are not part of an edge are
    set to zero and their derivatives are zero.  So the code could run
    faster if we took these out of the calculation. I think for most
    cases this is a small contribution.  Before I forced the initial
    conditions for these nonedges to be 0, they caused quite a bit of
    numerical headaches.
    '''
    N=G.order()
    X = V[0:N] #susceptibles
    Y = V[N:2*N] #infecteds
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) #there are places where we divide by X[i] which may = 0.
                                                         #In those cases the numerator is (very) 0, so it's easier
                                                         #to set this up as mult by inverse with a dummy value when
                                                         #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[2*N: 2*N+N**2]
    XX = V[2*N+N**2:]

    #print X.shape, Y.shape, XY.shape, XX.shape, N
    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T #not really needed, but helps keep consistent with equations as written.

    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, but I think this sequence of for loops is easier to read, or at least understand.
    #I expect it to run quickly regardless.  Will avoid (premature) optimization for now.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dX[i] += -trans_rate(u,v)*XY[i,j]
            dY[i] += trans_rate(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate(u,v)+rec_rate(v))*XY[i,j] 

            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(v):
                if w == u: #skip these
                    continue
                #so w != u.  
                k= index_of_node[w]
                #i corresponds to u, j to v and k to w.
                dXY[i,j] += trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  
                dXX[i,j] += -trans_rate(v,w) * XX[i,j] * XY[j,k]*Xinv[j] 
            for w in G.neighbors(u):
                if w == v:
                    continue #skip these
                k = index_of_node[w]
                dXY[i,j] += -  trans_rate(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -  trans_rate(u,w) * YX[k,i] * XX[i,j]*Xinv[i]

    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dX[:, None], dY[:,None], dXY, dXX), axis=0).T[0]
    print t, sum(dV), sum(dX)+sum(dY), sum(dXY), sum(dXX)
    return dV


def SIS_pair_based(G, nodelist, Y0, tau, gamma=None, XY0=None, XX0 = None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    r'''
    Encodes System (3.26) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIS disease model spreading on a given graph.  It captures the dependence with pairs, but not triples.

    It does not include corrections for triangles (or any other cycles).  The corrections for triangles are provided in the text, but not implemented here.

    There are some inefficiencies in the implementation:
        we track all pairs, rather than just those pairs in edges, but this is unlikely to significantly affect the calculation time.  This makes it much easier to vectorize things.
        We track pairs in both directions: e.g., XX[1,2] and XX[2,1].


    INPUT:
    ------
    G : Networkx graph
    nodelist : list
         list of nodes in G in the same order as in Y0
    Y0 : scipy array
         the array of initial infection probabilities for each node in order as in nodelist
    tau : number
          transmission rate of disease
    gamma : number (default None)
            global recovery rate  (incompatible with recovery_label!=None)
    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is infected.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    edge_label : string
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}
    recovery_label : string
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error.
    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    RETURNS
    -------
    if return_full_data is True:
        returns times, S, I, Xs, Ys, XY, XX
    if False:
        returns times, S, I

'''
    if gamma is None:
        if recovery_label is None:
            EoNError("need gamma or recovery_label defined.")
        else:
            rec_rate = lambda x : G.node[x][recovery_label]
    else:
        if recovery_label is not None:
            EoNError("only one of gamma and recovery_label can be defined.")
        rec_rate = lambda x : gamma
    if edge_label is not None:
        trans_rate = lambda x, y: tau*G.edge[x][y][edge_label]
    else:
        trans_rate = lambda x, y: tau

    times = scipy.linspace(tmin,tmax,tcount)


    N = len(Y0)
    X0=1-Y0

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    else:
        if XY0.shape != (N,N):
            EoNError("incompatible lengths for XY0 and Y0")

    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    else:
        if XX0.shape != (N,N):
            EoNError("incompatible lengths for XX0 and Y0")
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs that aren't
    XX0 = XX0*A  #in an edge, but this led to the error with odeint.  Multiplying by
                 #A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((Y0[:,None], XY0, XX0), axis=0).T[0]
    print 'V0shape', V0.shape
    index_of_node = {node:i for i, node in enumerate(nodelist)}

    V = integrate.odeint(_dSIS_pair_based_, V0, times, args = (G, nodelist, index_of_node, trans_rate, rec_rate))
    Ys = V.T[0:N]
    I = Ys.sum(axis=0)
    Xs = scipy.ones(N)[:,None]-Ys
    S = Xs.sum(axis=0)
    if return_full_data:
        XY = V.T[N: N+N**2]
        XX = V.T[N+N**2:]
        XY.shape = (N,N,tcount)
        XX.shape = (N,N,tcount)
        return times, S, I, Xs, Ys, XY, XX
    else:
        return times, S, I






def SIR_pair_based(G, nodelist, Y0, tau, gamma=None, X0 = None, XY0=None, XX0 = None, tmin = 0, tmax = 100, tcount = 1001, edge_label=None, recovery_label=None, return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIR disease model spreading on a given graph.  It captures the dependence with pairs, but not triples.
    It will be exact for a tree.

    There are NO CORRECTIONS for the existence of TRIANGLES or any other CYCLES.
    Some corrections for triangles are provided in the text, but not implemented here.
    
    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology
    

    INPUT:
    ------
    G : Networkx graph
    nodelist : list
         list of nodes in G in the same order as in Y0
    Y0 : scipy array
         the array of initial infection probabilities for each node in order as in nodelist
    tau : number
          transmission rate of disease
    gamma : number (default None)
            global recovery rate  (incompatible with recovery_label!=None)
    X0 : scipy array (default None)
            probability a random node is initially susceptible.
            the probability of initially recovered will be 1-X0-Y0.  By default we assume no
            initial recoveries, so X0=1-Y0 will be assumed in this case.
    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is infected.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    edge_label : string
            the label for a weight given to the edges.
            G.edge[i][j][edge_label] = g_{ij}
    recovery_label : string
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error.
    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    RETURNS
            if return_full_data is True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R
    '''


    if gamma is None:
        if recovery_label is None:
            EoNError("need gamma or recovery_label defined.")
        else:
            rec_rate = lambda x : G.node[x][recovery_label]
    else:
        if recovery_label is not None:
            EoNError("only one of gamma and recovery_label can be defined.")
        rec_rate = lambda x : gamma
    if edge_label is not None:
        trans_rate = lambda x, y: tau*G.edge[x][y][edge_label]
    else:
        trans_rate = lambda x, y: tau

    times = scipy.linspace(tmin,tmax,tcount)
    if X0 is None:
        X0 = 1-Y0
    N = len(Y0)

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    else:
        if XY0.shape != (N,N):
            EoNError("incompatible lengths for XY0 and Y0")
    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    else:
        if XX0.shape != (N,N):
            EoNError("incompatible lengths for XX0 and Y0")
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs that aren't
    XX0 = XX0*A  #in an edge, but this led to the error with odeint.  Multiplying by
                 #A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((X0[:,None], Y0[:,None], XY0, XX0), axis=0).T[0]
    #print V0.shape
    index_of_node = {node:i for i, node in enumerate(nodelist)}
    #    index_of_node = {}
    #    for i, node in enumerate(nodelist):
    #        index_of_node[node] = i


    V = integrate.odeint(_dSIR_pair_based_, V0, times, args = (G, nodelist, index_of_node, trans_rate, rec_rate))#, mxstep=10)#(times[1]-times[0])/1000)
    Xs = V.T[0:N]
    S = Xs.sum(axis=0)
    Ys = V.T[N:2*N]
    I = Ys.sum(axis=0)
    Zs = scipy.ones(N)[:,None]-Xs-Ys
    R = Zs.sum(axis=0)
    if return_full_data:
        XY = V.T[2*N: 2*N+N**2]
        XX = V.T[2*N+N**2:]
        #print len(XX)
        XY.shape = (N,N,tcount)
        XX.shape = (N,N,tcount)
        return times, S, I, R, Xs, Ys, Zs, XY, XX
    else:
        return times, S, I, R


def _dSIR_pair_based2_(V, t, G, nodelist, index_of_node, edgelist, index_of_edge, trans_rate, rec_rate):
    #X, Y, XY, YX, XX
    N=len(nodelist)
    E=len(edgelist)
    
    X = V[0:N]
    Y = V[N:2*N]
    def Xinv(u):
        '''there are places where we divide by X[i] which may = 0.  In those cases the numerator is (very) 0, so it's easier
        to set this up as mult by inverse with a dummy value when it is 1/0.
        '''
        index = index_of_node[u]
        if X[index]==0:
            return 0
        else:
            return 1/X[index]
    def Yinv(u):
        index = index_of_node[u]
        if Y[index]==0:
            return 0
        else:
            return 1/Y[index]
    #Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
    #Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[2*N: 2*N+E]
    YX = V[2*N+E:2*N+2*E]
    XX = V[2*N+2*E:]

    #print X.shape, Y.shape, XY.shape, XX.shape, N

    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    dXY = scipy.zeros(E)
    dYX = scipy.zeros(E)
    dXX = scipy.zeros(E)

    def my_XY(u,v): #either XY[i,j] or YX[i,j], depending on which is defined.  Might be cleaner with Try, Except
        if index_of_edge.has_key((u,v)):
            index = index_of_edge[(u,v)]
            return XY[index]
        else:
            index = index_of_edge[(v,u)]
            return YX[index]
    def my_YX(u,v):
        return my_XY(v,u)
    def my_XX(u,v):
        if index_of_edge.has_key((u,v)):
            index = index_of_edge[(u,v)]
        else:
            index = index_of_edge[(v,u)]
        return XX[index]
    
    #I could make the below more efficient, but I think this sequence of for loops is easier to read, or at least understand.
    #I expect it to run quickly regardless.  Will avoid (premature) optimization for now.
    for i, u in enumerate(nodelist):
        dY[i] += -rec_rate(u)*Y[i] 
    for edgeindex, (u,v) in enumerate(edgelist):
        i = index_of_node[u]
        j = index_of_node[v]

        dX[i] += -trans_rate(u,v)*XY[edgeindex] 
        dY[i] += trans_rate(u,v)*XY[edgeindex]  

        dX[j] += -trans_rate(u,v)*YX[edgeindex]
        dY[j] += trans_rate(u,v)*YX[edgeindex]

        dXY[edgeindex] += - (trans_rate(u,v)+rec_rate(v))*XY[edgeindex]
        dYX[edgeindex] += - (trans_rate(u,v)+rec_rate(v))*YX[edgeindex]

        for w in G.neighbors(u):
            if w==v:  #skip this case
                continue
            k=index_of_node[w]
            dXY[edgeindex] += -trans_rate(w,u)*my_YX(w,u) * my_XY(u,v)*Xinv(u)
            dYX[edgeindex] += trans_rate(w,u)*my_XX(v,u)*my_XY(u,w)*Xinv(u)
            dXX[edgeindex] += -trans_rate(w,u) * my_XX(u,v) * my_XY(u,w)*Xinv(u) 
        for w in G.neighbors(v):
            if w==u:
                continue
            #w transmits to v.
            dXY[edgeindex] += trans_rate(w,v)*my_XX(u,v) *my_XY(v,w)*Xinv(v)
            dYX[edgeindex] += -trans_rate(w,v)*my_XY(v,w)*my_YX(u,v)*Xinv(v)
            dXX[edgeindex] += -trans_rate(w,v)*my_XX(u,v)*my_XY(v,w)*Xinv(v)

    dV = scipy.concatenate((dX, dY, dXY, dYX, dXX), axis=0)
    print t, sum(dV)
    return dV


def _SIR_pair_based_initialize_node_data(G, rho, nodelist, X0, Y0):
    #inputs must define either rho or Y0.  In first case nodelist is optional.
    if (rho and Y0 is not None) or (not rho and Y0 is None):
        raise EoNError("need rho or Y0 defined for initial condition, but not both")
    if Y0 is not None and nodelist is None:
        raise EoNError("order in Y0 is ambiguous if nodelist is not given.")

    if not nodelist: #then rho is defined but not Y0
        nodelist = list(G.nodes())
    if rho: #Y0 not defined
        Y0 = rho*scipy.ones(len(nodelist))
        X0 = 1-Y0 #assume X0=0
    else:  #Y0 is defined
        if not X0:
            X0=1-Y0  #assume Z0=0
        #otherwise X0 is given and Z0=1-X0-Y0
    return nodelist, X0, Y0

def _SIR_pair_based_initialize_edge_data(G, edgelist, nodelist, XY0, YX0, XX0, X0, Y0, index_of_node):
    if (not XY0 is None or YX0 is None or XX0 is None) and (XY0 is not None or YX0 !=None  or XX0 is not None):  #at least one defined and one not defined
        raise EoNError("must define all of XY0, YX0, and XX0 or none of them")
    if not edgelist:
        if XY0:
            raise EoNError("order in XY0, YX0, and XX0 is ambiguous if edgelist is not given.")
        else:
            edgelist = list(G.edges())

    if XY0:
        #test that XY0 <= X0Y0, same for  YX0 and XX0
        for index,(u,v) in enumerate(edgelist):
           i_u = index_of_node[u]
           i_v = index_of_node[v]
           if XY0[index] >X0[i_u]*Y0[i_v] or YX0[index]>Y0[i_u]*X0[I_v] or XX0[index]>X0[i_u]*X0[i_v]:
               raise EoNError("edge probabilities inconsistent with node probabilities")
    else:
        XY0 = scipy.array([X0[index_of_node[u]]*Y0[index_of_node[v]] for u,v in edgelist])
        YX0 = scipy.array([Y0[index_of_node[u]]*X0[index_of_node[v]] for u,v in edgelist])
        XX0 = scipy.array([X0[index_of_node[u]]*X0[index_of_node[v]] for u,v in edgelist])
    return edgelist, XY0, YX0, XX0

def get_rate_functions(G, tau, gamma=None, recovery_label=None, edge_rate_label = None):
    r'''
    INPUT:
    -----
    G : networkx Graph
        the graph disease spread on
    tau : number
        disease parameter giving edge transmission rate (subject to edge scaling)
    gamma : number (default None)
        disease parameter giving typical recovery rate, If given then recovery_label must be None
    recovery_label : key (default None)
        G.node[node][recovery_label] is recovery rate
    edge_rate_label : key (default None)
        G.edge[u][v][edge_rate_label] scales up or down the recovery rate.
'''
    if (gamma is None and recovery_label is None) or (gamma is not None and recovery_label is not None):
        raise EoNError("need exactly one of gamma and recovery_label defined.")
    
    if gamma:
        rec_rate = lambda x : gamma
    else:
        rec_rate = lambda x : G.node[x][recovery_label]
        
    if edge_rate_label:
        trans_rate = lambda x, y: tau*G.edge[x][y][edge_rate_label]
    else:
        trans_rate = lambda x, y: tau
    return rec_rate, trans_rate
    
def SIR_pair_based2(G, tau, gamma=None, rho = None, nodelist=None, X0=None, Y0=None, edgelist = None, XY0=None, YX0 = None, XX0 = None, tmin = 0, tmax = 100, tcount = 1001, edge_rate_label=None, recovery_label=None, return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIR disease model spreading on a given graph.  It captures the dependence with pairs, but not triples.
    It will be exact for a tree.

    There are NO CORRECTIONS for the existence of TRIANGLES or any other CYCLES.
    Some corrections for triangles are provided in the text, but not implemented here.

    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology

    <\dot{Y}_i> = tau \sum_j g_{ij} <XY>>  -  gamma_i <Y_i>
    <\dot{XY}> = tau sum_{k \neq i} g_{jk} <XX><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XY>/<Xi>
                   - tau g_{ij}<XY> - gamma_j <XIYj> 
    <\dot{XX}> = 
    <>
    The equations as coded involve all pairs rather than just the pairs that are in edges.  Those that are not part of an edge are set to zero and
    their derivatives are zero.  So the code could run faster, but I think for most cases this is a small contribution.  Before I forced the initial
    conditions for these nonedges to be 0, they caused quite a bit of numerical headaches.

    -------
    G : Networkx graph
    nodelist : list
         list of nodes in G in the same order as in Y0
    Y0 : scipy array
         the array of initial infection probabilities for each node in order as in nodelist
    tau : number
          transmission rate of disease
    gamma : number (default None)
            global recovery rate  (incompatible with recovery_label!=None)
    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is infected.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced randomly
            according to Y0.
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    edge_rate_label : string
            the label for a weight given to the edges.
            G.edge[i][j][edge_rate_label] = g_{ij}
    recovery_label : string
            a label for a weight given to the nodes for their recovery rates
            G.node[i][recovery_label] = gamma_i
            We cannot define both gamma and recovery_label.  This will raise an error.
    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R
    '''

    #note: we do not test whether the nodelist and edgelist are in fact lists of nodes or edges
    nodelist, X0, Y0 = _SIR_pair_based_initialize_node_data(G, rho, nodelist, X0, Y0)

    index_of_node = {node:i for i, node in enumerate(nodelist)}

    edgelist, XY0, YX0, XX0 = _SIR_pair_based_initialize_edge_data(G, edgelist, nodelist, XY0, YX0, XX0, X0, Y0, index_of_node)

    index_of_edge = {edge:i for i, edge in enumerate(edgelist)}
    N = len(nodelist)
    E = len(edgelist)
    #now we define functions which give the transmission rate of edges and recovery rate of
    #nodes.  
    rec_rate, trans_rate = get_rate_functions(G, tau, gamma, recovery_label, edge_rate_label)

    times = scipy.linspace(tmin,tmax,tcount)

    
    V0 = scipy.concatenate((X0, Y0, XY0, YX0, XX0),axis=0)
    V = integrate.odeint(_dSIR_pair_based2_, V0, times, args = (G, nodelist, index_of_node, edgelist, index_of_edge, trans_rate, rec_rate))#, mxstep=10)#(times[1]-times[0])/1000)

    #dX, dY, dXY, dYX, dXX
    Xs = V.T[0:N]
    S = Xs.sum(axis=0)
    Ys = V.T[N:2*N]
    I = Ys.sum(axis=0)
    Zs = scipy.ones(N)[:,None]-Xs-Ys
    R = Zs.sum(axis=0)
    if return_full_data:
        XY = V.T[2*N: 2*N+E]
        YX = V.T[2*N+E:2*N+2*E]
        XX = V.T[2*N+2*E:]
        YY = 1 - XY - YX-XX
        return times, S, I, R, Xs, Ys, Zs, XY, YX, XX, YY,edgelist, nodelist
    else:
        return times, S, I, R




######    HOMOGENEOUS MEANFIELD

def _dSIS_homogeneous_meanfield_(X, t, n_over_N, tau, gamma):
    S, I = X
    dSdt = gamma * I - tau*n_over_N*S*I
    dIdt = -dSdt
    dX = scipy.array([dSdt,dIdt])
    return dX

def _dSIR_homogeneous_meanfield_(X, t, n_over_N, tau, gamma):
    S, I = X
    dSdt = -tau*n_over_N*S*I
    dIdt = tau*n_over_N*S*I - gamma*I
    dX = scipy.array([dSdt, dIdt])
    return dX

           
def SIS_homogeneous_meanfield(S0, I0, n, tau, gamma, tmin=0, tmax=100, tcount=1001):
    '''Encodes System (4.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{S}] = \gamma [I] - tau n[S][I]/N
    [\dot{I}] = \tau n[S][I]/N - \gamma [I]


    INPUTS
    --------
    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    n : integer
        degree of all nodes.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports

    '''

    N=S0+I0
    X0=scipy.array([S0,I0])
    times=scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_meanfield_, X0,times, args=(float(n)/N, tau, gamma))
    S, I= X.T
    return times, S, I

def SIR_homogeneous_meanfield(S0, I0, R0, n, tau, gamma, tmin=0, tmax=100, tcount=1001):
    '''Encodes System (4.9) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{S}] = - tau n[S][I]/N
    [\dot{I}] = \tau n[S][I]/N - \gamma [I]
    [\dot{R}] = \gamma [I]


    INPUTS
    --------
    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    R0 : number
         initial number recovered
    n : integer
        degree of each node
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.
    
    
    '''

    N=S0+I0+R0
    X0= scipy.array([S0,I0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_meanfield_, X0, times, args=(float(n)/N, tau, gamma))
    S, I= X.T
    R = N-S-I
    return times, S, I, R


####   HOMOGENEOUS PAIRWISE

def _dSIS_homogeneous_pairwise_(X, t, N, n, tau, gamma):
    S, SI, SS = X 
    I = N-S
    II = N*n - SS - 2*SI
    nm1_over_n = (n-1.)/n
    dSdt = gamma * I - tau*SI
    #dIdt = -dSdt
    dSIdt = gamma*(II-SI) + tau*nm1_over_n*SI*(SS-SI)/S - tau*SI
    dSSdt = 2*gamma*SI - 2*tau*nm1_over_n*SI*SS/S
    #dIIdt = -2*gamma*II + 2*tau*nm1_over_n*SI**2/S + 2*tau*SI
    dX = scipy.array([dSdt,dSIdt,dSSdt])
    return dX 


def _dSIR_homogeneous_pairwise_(X, t, n, tau, gamma):
    S,I,SI,SS = X
    nm1_over_n = (n-1.)/n
    dSdt = -tau*SI
    dIdt = tau*SI - gamma*I
    dSIdt = -gamma*SI + tau*nm1_over_n * SI*(SS-SI)/S - tau*SI
    dSSdt = -2*tau*nm1_over_n*SI*SS/S
    dX =  scipy.array([dSdt,dIdt,dSIdt,dSSdt])
    return dX

def SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    r'''Encodes System (4.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{S}] = gamma [I] - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I] = -[\dot{S}]
    [\dot{SI}] = \gamma([II]-[SI])+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] - \tau [SI]
    [\dot{SS}] = 2\gamma[SI] - 2\tau ((n-1)/n) [SI][SS]/[S]
    [\dot{II}] = -2\gamma[II] + 2\tau((n-1)/n) [SI]^2/[S] + 2\tau[SI]

    conserved quantities: [S]+[I]
                          [SS]+2[SI]+[II]

    n([S]+[I]) should equal [SS]+2[SI]+[II].  So II will be calculated based on this.
    
    INPUTS
    --------
    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    SI0 : number
          initial number of SI edges
    SS0 : number
          initial number of SS edges
    n : integer
          (common) degree of nodes.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.
    
    '''
    N = S0+I0
                #twoM = SS0+2*SI0+II0 #M is number of edges 
		#    n = float(twoM)/N   #allowing for float is better than forcing an integer.  
                                 #it's reasonable that this approximate model might be
                                 #used for such a case.  Especially since the equations
                                 #should permit for use of proportions rather than integers.
    X0 = scipy.array([S0, SI0, SS0])
    
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_pairwise_, X0, times, args=(N, n, tau, gamma))
    S, SI, SS= X.T
    I = N-S
    
    if return_full_data:
	II = N*n - SS-2*SI
        return times, S, I, SI, SS, II
    else:
        return times, S, I
    
    
def SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes System (4.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{S}] = - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I]
    [\dot{R}] = \gamma [I]    ;    [R] = N-[S]-[I]
    [\dot{SI}] = -\gamma [SI]+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] - \tau [SI]
    [\dot{SS}] = - 2\tau ((n-1)/n) [SI][SS]/[S]

    conserved quantities: [S]+[I]+[R]  also [SS]+[II]+[RR] + 2([SI] + [SR] + [IR])

    INPUTS
    ---------
    S0 : Initial number suusceptible
    I0 : Initial number infected
    R0 : Initial number recovered
    SI0 : Initial number of SI edges
    SS0 : Initial number of SS edges
    n : Degree of nodes
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I, R or all calculated data.
                       if True, then returns times, S, I, R, SI, SS
    '''
    N = S0+I0+R0
    #    if n*N < SS0+2*SI0:
    #	    print("Warning: for these inputs, require a negative number of edges that include R nodes, blindly running ahead")
    X0 = scipy.array([S0, I0, SI0, SS0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_pairwise_, X0, times, args=(n, tau, gamma))
    S, I, SI, SS = X.T
    R = N-S-I
    if return_full_data:
        return times, S, I, R, SI, SS
    else:
        return times, S, I, R


def SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()
    S0 = (1-rho)*N
    I0 = rho*N
    SI0 = (1-rho)*N*n*rho
    SS0 = (1-rho)*N*n*(1-rho)
    return SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin, tmax, tcount, return_full_data)

def SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()
    S0 = (1-rho)*N
    I0 = rho*N
    SI0 = (1-rho)*N*n*rho
    SS0 = (1-rho)*N*n*(1-rho)
    R0=0
    return SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin, tmax, tcount, return_full_data)




#######     HETEROGENEOUS MEAN FIELD

def _dSIS_heterogeneous_meanfield_(X, t, kcount, tau, gamma):
    ks = scipy.arange(kcount)
    S = scipy.array(X[:kcount])
    I = scipy.array(X[kcount:])
    pi_I = ks.dot(I)/ ks.dot(I+S)

    Sdot = gamma*I - tau * ks*S*pi_I
    Idot = tau*ks*S*pi_I - gamma*I
    dX = scipy.concatenate((Sdot,Idot), axis=0)
    return dX
def _dSIR_heterogeneous_meanfield_(X, t, S0, Nk, tau, gamma):
    theta = X[0]
    Rk = X[1:]
    ks = scipy.arange(len(Rk))
    Sk = S0*(theta**ks)
    Ik = Nk - Sk - Rk
    pi_I = ks.dot(Ik)/ks.dot(Nk)
    dRkdt = gamma*Ik
    dThetadt = - tau *pi_I * theta

    dX = scipy.concatenate(([dThetadt],dRkdt), axis=0)
    return dX




def SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes System (5.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    Sk0 is an array (or a list). It is not a dict.  Sk0[k] is the *number* of
    nodes that are susceptible and have degree k (even if some degrees 
    missing).  A dict like this can be converted into an array by
    Sk0 = scipy.array([Sk0dict.get(k,0) for k in xrange(max(Sk0dict.keys())+1)])

    Ik0 is similar to Sk0.


    [\dot{S}_k] = \gamma [I_k] - \tau k [S_k] \pi_I
    [\dot{I}_k] = -(above)
    \pi_I = \sum_k k [I_k] / \sum_k  k [N_k]



    INPUTS
    ---------
    Sk0 : scipy array
          number susceptible for each k
    Ik0 : scipy array
          number infected for each k
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I or all calculated data.
                       if True, returns t, S, I, Sk, Ik

    '''
    if len(Sk0) != len(Ik0):
        raise EoNError('length of Sk0 not equal to length of Ik0')

    kcount = len(Sk0)
    
    X0 = scipy.concatenate((Sk0,Ik0), axis=0)
    Nk = Sk0+Ik0
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIS_heterogeneous_meanfield_, X0, times, args=(kcount,tau,gamma))
    Sk = scipy.array(X.T[:kcount])
    Ik = scipy.array(X.T[kcount:])
    S = Sk.sum(axis=0)
    I = Ik.sum(axis=0)
    if not return_full_data:
	return times, S, I
    else:
	return times, S, I, Sk, Ik

    
def SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes System (5.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    Ik0 and Rk0 are similar to Sk0.

    [S_k] = [S_k](0) theta^k
    [I_k] = [N_k] - [S_k] - [R_k]
    [\dot{R}_k] = \gamma [I_k]
    pi_I = \sum_k k[I_k]


    INPUTS
    ------------
    Sk0 : array
          Sk0[k] is the number of
    nodes that are susceptible and have degree k (even if some degrees 
    missing).

    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.

    RETURNS
    -------
    if return_full_data is False
        times, S, I, R (all scipy arrays)
    if True,
        times, Sk, Ik, Rk (the Xk are scipy 2D arrays)
    
    '''
    if len(Sk0) != len(Ik0) or len(Sk0) != len(Rk0):
        raise EoNError('length of Sk0, Ik0, and Rk0 must be the same')

    theta0=1
    Nk = Sk0+Ik0 +Rk0
    X0 = scipy.concatenate(([theta0],Rk0), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_heterogeneous_meanfield_, X0, times, args = (Sk0, Nk, tau, gamma))
    
    theta = X.T[0]
    Rk = X.T[1:]
    ks = scipy.arange(len(Rk))
    L=(theta[None,:]**ks[:,None])
    Sk=Sk0[:,None]*L
    Ik = Nk[:,None]-Sk-Rk
    if not return_full_data:
	return times, sum(Sk), sum(Ik), sum(Rk)
    else:
	return times, Sk, Ik, Rk
    
def SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Takes a graph and an initial proportion infected rho.  Calculates Sk0 and Ik0 and calls the heterogeneous meanfield model'''
    if rho is None:
        rho = 1./G.order()
    #Pk = get_Pk(G)
    #maxk = max(Pk.keys())
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)#G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin, tmax, tcount, return_full_data)

def SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Takes a graph and an initial proportion infected rho.  Calculates Sk0 and Ik0 and calls the heterogeneous meanfield model

    RETURNS
    -------
    if return_full_data is False
        times, S, I, R (all scipy arrays)
    if True,
        times, Sk, Ik, Rk (the Xk are scipy 2D arrays)
    '''
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)#G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False)



#######      HETEROGENEOUS PAIRWISE
def _dSIS_heterogeneous_pairwise_(X, t, Nk, NkNl, tau, gamma, Ks):
    '''
    Gives the derivatives

    INPUTS
    ------
    X : current values of variables
    t : current time
    Nk : The number of nodes of each degree; Nk[i] is the number of nodes with degree Ks[i]
    NkNl : number of edges of various types.  NkNl[i,j] corresponds to Ks[i] and Ks[j].
    tau : transmission rate
    gamma : recovery rate
    Ks : a scipy array --- gives the observed degrees in increasing order.
    '''
    kcount = len(Ks)
    Sk = X[:kcount]
    SkSl = X[kcount: kcount + kcount**2]
    SkIl = X[kcount+kcount**2: kcount + 2*kcount**2]

    SkSl.shape=(kcount,kcount)
    SkIl.shape=(kcount,kcount)

    Ik = Nk - Sk
    SkI = SkIl.sum(1)#sum(SkIl.T)
    IkSl = SkIl.T
    IkIl = NkNl - SkSl - SkIl - IkSl

    Ls = Ks #easier to keep as l for readability of eqns

    SlI = SkI #easier in this form for readability
    ISk = SkI

    tmpSk = 1.*Sk #the 1.* is to make it a copy, not the same object
    tmpSk[tmpSk==0] = 1

    SkSlI = SkSl * ((Ls-1))*SlI/ (Ls*tmpSk)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (Ks*tmpSk)).T
    ISkSl = SkSlI.T 
    IkSlI = ISkIl.T 
    
    dSk = gamma*Ik - tau *SkI
    dSkSl = gamma*(SkIl + IkSl) - tau*(SkSlI + ISkSl)
    dSkIl = gamma*(IkIl-SkIl) + tau*(SkSlI - ISkIl - SkIl)

    dSkIl.shape=(kcount**2,1)
    dSkSl.shape=(kcount**2,1)
    dX = scipy.concatenate((dSk[:,None], dSkSl, dSkIl),axis=0).T[0]
    return dX

def _dSIR_heterogeneous_pairwise_(X, t, tau, gamma, Nk, Ks):
    kcount = len(Ks)
    Sk = X[:kcount]
    Ik = X[kcount:2*kcount]
    SkSl = X[2*kcount:2*kcount+kcount**2]
    SkIl = X[2*kcount+kcount**2:2*kcount+2*kcount**2]

    SkSl.shape=(kcount,kcount)
    SkIl.shape=(kcount,kcount)

    SkI = SkIl.sum(1)
    IkSl = SkIl.T

    Ls = Ks #easier to keep as l for readability of eqns
    SlI = SkI
    ISk = SkI

    tmpSk = 1.*Sk #the 1.* is to make it a copy, not the same object
    tmpSk[tmpSk==0] = 1
    tmpSl = tmpSk

    SkSlI = SkSl * ((Ls-1))*SlI/ (Ls*tmpSl)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (Ks*tmpSk)).T
    ISkSl = SkSlI.T 
    IkSlI = ISkIl.T 

    dSk = -tau*SkI
    dIk = tau*SkI - gamma*Ik
    dSkIl = -gamma*SkIl + tau*(SkSlI - ISkIl - SkIl)
    dSkSl = -tau*(SkSlI + ISkSl)

    dSkIl.shape=(kcount**2,1)
    dSkSl.shape = (kcount**2,1)

    dX = scipy.concatenate((dSk[:,None], dIk[:,None], dSkSl, dSkIl),axis=0).T[0]
    return dX



def SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data = False, Ks = None):
    '''Encodes System (5.13) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{Sk}] = gamma [Ik] - tau [SkI]
    [\dot{Ik}] = tau [SkI] - gamma [Ik]   = - [\dot{Sk}]
    [\dot{SkIl}] = gamma([IkIl] - [SkIl]) + tau([SkSlI] - [ISkIl] - [SkIl])
    [\dot{SkSl}] = gamma([SkIl]+[IkSl]) - tau ([SkSlI] + [ISkSl])
    [\dot{IkIl}] = -2 gamma[IkIl] + tau ([SkIl]+[IkSl]+[ISkIl]+[IkSlI])
    [AlSkI] = ((k-1)/k) [AlSk][SkI]/[Sk]
    [ISkAl] = ((k-1)/k) [ISk][SkAl]/[Sk]

    So: [SkSlI] = ((l-1)/l) [SkSl][SlI]/[Sl] 
        [ISkIl] = ((k-1)/k) [ISk] [SkIl]/[Sk]
        [ISkSl] = ((k-1)/k) [ISk] [SkSl]/[Sk]
        [IkSlI] = ((l-1)/l) [IkSl][SlI]/[Sl]

    conserved quantities : [Sk]+[Ik],  [SkSl] + [SkIl] + [IkIl] + [IkSl]

    identities: [IkSl] = [SlIk], so where IkSl needed, use SlIk.T

    INPUTS
    ------------
    Sk0 : array.  Sk0[k] is the number of
                 nodes that are susceptible and have degree k.  If one is empty, it becomes 0.
                (if Ks is defined, the definition changes slightly, see below)

    Ik0 : array
         similar to Sk0, but for infected.
                (if Ks is defined, the definition changes slightly, see below)

    SkSl0, SkIl0, and IkIl0 : 2D arrays
         SkIl0[k][l] is [S_kI_l] - see below for constraints these should satisfy
         related to Sk0 and Ik0.  The code does not enforce these constraints.
                (if Ks is defined, the definition changes slightly, see below)

    tau : number
          transmission rate

    gamma : number
            recovery rate

    tmin : number (default 0)
           minimum report time

    tmax : number (default 100)
           maximum report time 

    tcount : integer (default 1001)
             number of reports

    return_full_data : boolean (default False)
                       If True, return times, Sk, Ik, SkIl, SkSl, IkIl
                       If False, return times, S, I

    Ks : scipy array. (default None)
         (helps prevent memory errors) if some degrees are not
         observed, then the corresponding entries of these arrays are
         zero.  This can lead to memory errors in the case of a
         network with many missing degrees.  So Ks is an (assumed)
         ordered vector stating which Ks are actually observed.  Then
         the Sk0[i] is the number of nodes that are susceptible and
         have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.

    In principle, there are constraints relating Sk with SkSl and SkIl and similarly relating Ik with IkIl and SkIl.T.  No attempt is made to enforce these.  It is assumed the user will ensure acceptible inputs.

    We could also remove Sk0 and Ik0 as inputs and infer them from the others, but for consistency with elsewhere, this is not done here.
    '''

    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
    times = scipy.linspace(tmin,tmax,tcount)
    Nk = Sk0+Ik0
    kcount = len(Nk)
    NkNl = SkSl0 + SkIl0 + IkIl0 + SkIl0.T
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)
    X0 = scipy.concatenate((Sk0[:,None], SkSl0, SkIl0), axis=0).T[0]

    X = my_odeint(_dSIS_heterogeneous_pairwise_, X0, times, args = (Nk, NkNl, tau, gamma, Ks))

    kcount = len(Nk)
    Sk = X.T[:kcount]
    S = Sk.sum(axis=0)
    Ik = Nk[:,None] - Sk
    I = Ik.sum(axis=0)
    if return_full_data:
	SkSl = X.T[kcount:kcount+kcaount**2]
	SkIl = X.T[kcount+kcount**2:]
        SkSl.shape = (kcount,kcount,tcount)
        SkIl.shape = (kcount,kcount,tcount)
	IkIl = NkNl - SkSl - SkIl - SkIl.T
	return times, S, I, Sk, Ik, SkIl, SkSl, IkIl
    else:
	return times, S, I

def SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False, Ks = None):
    '''Encodes System (5.15) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [\dot{S}_k] = -tau [S_k I]
    [\dot{I}_k] = tau [S_k I] - gamma [I_k]
    [\dot{R}_k] = gamma [I_k]  (but using Rk=Nk-Sk-Ik for this equation)
    [\dot{S_kI_l}] = -gamma[S_k I_l] + tau([S_k S_l I] - [I S_k I_l] - [S_k I_l])
    [\dot{S_kS_l}] = -tau([S_k S_l I] + [I S_k S_l])

    [A_l S_k I] = ((k-1)/k) [A_l S_k] [S_k I]/ [S_k]
    [I S_k A_l] = ((k-1)/k) [I S_k] [S_k A_l]/ [S_k]

    INPUTS
    ----------
    Sk0 : scipy array, Sk0[k] is number of degree k susceptible at time 0.
    Ik0 : scipy array
    Rk0 : scipy array
    SkIl0 : scipy 2D array
    SkSl0 : scipy 2D array

    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       If True, return times, Sk, Ik, Rk, SkIl, SkSl
                       If False, return times, S, I, R


    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.
    Ks : scipy array. (default None)
         (helps prevent memory errors) if some degrees are not
         observed, then the corresponding entries of these arrays are
         zero.  This can lead to memory errors in the case of a
         network with many missing degrees.  So Ks is an (assumed)
         ordered vector stating which Ks are actually observed.  Then
         the Sk0[i] is the number of nodes that are susceptible and
         have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.        

'''

    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
#    print "Ks is ", Ks
    times = scipy.linspace(tmin,tmax,tcount)

    Nk = Sk0+Ik0+Rk0
    kcount = len(Ks)
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)

    X0 = scipy.concatenate((Sk0[:,None], Ik0[:,None], SkSl0, SkIl0), axis=0).T[0]
    X = integrate.odeint(_dSIR_heterogeneous_pairwise_, X0, times, args = (tau, gamma, Nk, Ks))

    Sk = X.T[:kcount]
    S = Sk.sum(axis=0)
    Ik = X.T[kcount:2*kcount]
    I = Ik.sum(axis=0)
    SkIl = X.T[2*kcount:2*kcount+kcount**2]
    SkSl = X.T[2*kcount+kcount**2: 2*kcount+2*kcount**2]

    Rk = Nk[:,None] - Sk - Ik
    R = Rk.sum(axis=0)
    
    if return_full_data:
        SkIl.shape = (kcount,kcount,tcount)
        SkSl.shape = (kcount,kcount,tcount)
        return times, S, I, R, Sk, Ik, Rk, SkIl, SkSl
    else:
        return times, S, I, R


def SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data = False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)#G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])

    NkNl, SkSl0, SkIl0, IkIl0, Ks = get_NkNl_and_IC_as_arrays(G, rho, withKs = True, SIR=False)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    return SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, gamma, tmin, tmax, tcount, return_full_data, Ks = scipy.array(Ks))

    
def SIR_heterogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin=0, tmax=100, tcount = 1001, return_full_data=False):
    if rho==None:
        rho = 1./G.order()    
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    NkNl, SkSl0, SkIl0, Ks = get_NkNl_and_IC_as_arrays(G, rho, withKs = True, SIR = True)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    Rk0 = scipy.array([Rk0[k] for k in Ks])
    return SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma, tmin = tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data, Ks = Ks)




#######    COMPACT PAIRWISE


def _dSIS_compact_pairwise_(X, t, Nk, twoM, tau, gamma):
    SI, SS= X[-2:]
    Sk = X[:-2]

    Ik = Nk - Sk
    II = twoM - SS - 2*SI

    ks = scipy.arange(len(Sk))
    SX = float(SI + SS)
    Q = (1./SX**2) * (ks*(ks-1)).dot(Sk)

    dSk = gamma *Ik - tau * ks *Sk *SI/SX
    dSI = gamma*(II-SI) + tau*(SS-SI)*(SI)*Q - tau*SI
    dSS = 2*gamma*SI - 2*tau*SS*SI*Q

    dX = scipy.concatenate((dSk, [dSI, dSS]), axis=0)
    return dX

def _dSIR_compact_pairwise_(X, t, N, tau, gamma):
    SS, SI, R = X[-3:]
    Sk = X[:-3]
    ks = scipy.arange(len(Sk))
    SX = ks.dot(Sk)
    Q = (ks*(ks-1)).dot(Sk)/SX**2
    I = N - sum(Sk) - R
    
    dSk = -tau * ks * Sk *SI/SX
    dSS = -2*tau*SS*SI*Q
    dSI = -gamma*SI + tau*(SS-SI)*SI*Q-tau*SI
    dR = gamma*I

    dX = scipy.concatenate((dSk, [dSS,dSI,dR]),axis=0)
    return dX 


def SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.18) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [dot{S}_k] = gamma[I_k] - tau k [S_k] [SI]/[SX]
    [dot{I}_k] = tau * k *[S_k] [SI]/SX - gamma [I_k] = -[dot{S}_k]
    [dot{SI}] = gamma([II]-[SI]) + tau([SS]-[SI])[SI]Q - tau[SI]
    [dot{SS}] = 2 gamma[SI] - 2 tau [SS] [SI] Q
    [dot{II}] = 2 tau[SI] = 2 gamma[II] + 2 tau[SI]^2Q
    [SX] = sum_k k [S_k]
    Q = (1/[SX]^2) sum_k (k-1)k[S_k]

    Input
    Sk0 : scipy array
          number susceptible for each k
    Ik0 : scipy array
          number infected for each k
    SI0 : number
          number of SI edges
    SS0 : number
          number of SS edges
    II0 : number
          number of II edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       if True, return times, S, I, Sk, Ik, SI, SS, II
		       if False,  return times, S, I




    conserved quantities:  [Sk]+[Ik]     ;     SS + II + 2SI
    '''
    Nk = Sk0+Ik0
    twoM = SS0+II0+2*SI0
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.concatenate((Sk0, [SI0,SS0]), axis=0)
    X = integrate.odeint(_dSIS_compact_pairwise_, X0, times, args = (Nk, twoM, tau, gamma))

    SI, SS = X.T[-2:]
    Sk = X.T[:-2]
    S = Sk.sum(axis=0)

    Ik = Nk[:,None] - Sk
    I = Ik.sum(axis=0)

    II = twoM - SS - 2*SI

    if return_full_data:
	return times, S, I, Sk, Ik, SI, SS, II
    else:
	return times, S, I

def SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.19) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [dot{S}_k] = -tau k [S_k] [SI]/[SX]
    [dot{SS}] = -2 tau [SS] [SI] Q
    [dot{SI}] = -gamma [SI] + tau([SS]-[SI])[SI]Q - tau [SI]
    [dot{R} = gamma [I]
    [SX] = sum_k k[S_k]
    Q = (1/[SX]^2) sum_k (k-1) k [S_k]
    [S] = sum [S_k]
    I = N-[S]-R

    INPUTS
    ------
    Sk0 : scipy array
          initial number of suscetibles of each degree k
    I0 : number
         initial number infected
    R0 : number
         initial number recovered
    SS0 : number
          initial number of SS edges
    SI0 : number
          initial number of SI edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.

    RETURNS
    -------
    '''
    
    times = scipy.linspace(tmin,tmax,tcount)
    N = I0+R0+sum(Sk0)
    X0 = scipy.concatenate((Sk0, [SS0, SI0, R0]), axis=0)
    X = integrate.odeint(_dSIR_compact_pairwise_, X0, times, args = (N, tau, gamma))
    SI, SS, R = X.T[-3:]
    Sk = X.T[:-3]
    S = Sk.sum(axis=0)
    I = N - R - S
    if return_full_data:
	return times, Sk, I, R, SS, SI
    else:
	return times, S, I, R


def SIS_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    N = G.order()
    maxk = max(Pk.keys())
    Nk = G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    Sk0 = (1-rho)*Nk
    Ik0 = rho*Nk
    SI0 = sum(Nk[k]*k*(1-rho)*rho for k in range(maxk+1))
    SS0 = sum(Nk[k]*k*(1-rho)*(1-rho) for k in range(maxk+1))
    II0 = sum(Nk[k]*k*rho*rho for k in range(maxk+1))
    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, tmax, tcount, return_full_data)

    
def SIR_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]
    
    I0 = sum(Ik0)
    R0 = 0

    SX0 = scipy.dot(Sk0,ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0
    
    return SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)




#######SUPER COMPACT PAIRWISE

def _dSIS_super_compact_pairwise_(X, t, tau, gamma, N, k_ave, ksquare_ave, kcube_ave):
    '''    
    [\dot{I}] = tau [SI] - gamma [I]
    [\dot{SS}] = 2 gamma [SI] - 2 tau [SI] [SS] Q
    [\dot{SI}] = gamma ([II]-[SI]) + tau [SI] ([SS]-[SI])Q - tau [SI]
    [\dot{II}] = -2 gamma [II] + 2 tau [SI]^2 Q + 2 tau [SI]
    Q = ((<K^2>(<K^2>-n_S<K>) + <K^3>(n_S-<K>))/(n_S(<K^2>-<K>^2)) - 1)/n_S[S]
    n_S = ([SI] + [SS])/(N-[I])
    '''

    I, SS, SI, II = X
    S = N-I
    
    n_S = (SS+SI)/(S)
    
    Q = ((ksquare_ave*(ksquare_ave-n_S*k_ave)+kcube_ave*(n_S-k_ave))/(n_S*(ksquare_ave-k_ave**2))-1)/(S*n_S)
    dIdt = tau*SI - gamma*I
    dSSdt = 2*gamma*SI - 2*tau*SI*SS*Q
    dSIdt = gamma*(II - SI) + tau*SI*(SS-SI)*Q - tau*SI
    dIIdt = -2*gamma*II + 2*tau*SI**2*Q + 2*tau*SI

    dX = scipy.array([dIdt, dSSdt, dSIdt, dIIdt])
    return dX

def _dSIR_super_compact_pairwise_(X, t, tau, gamma, psihat, psihatPrime, psihatDPrime, N):
    theta = X[0]
    SS = X[1]
    SI = X[2]
    R = X[3]
    S = N * psihat(theta)
    I = N - S - R

    Q = psihatDPrime(theta)/(N*psihatPrime(theta)**2)
    dThetadt = - tau*SI/(N*psihatPrime(theta))
    dSSdt = -2*tau*SS*SI*Q
    dSIdt = - gamma*SI + tau*(SS-SI)*SI*Q - tau*SI
    dRdt = gamma*I

    dX = scipy.array([dThetadt, dSSdt, dSIdt, dRdt])
    return dX


def SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, k_ave, ksquare_ave, kcube_ave, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.20) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    INPUTS
    -------

    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
             tells whether to just return times, S, I, R or all calculated data.

    RETURNS
    -------
    
    '''
    X0 = [I0, SS0, SI0, II0]
    N = S0+I0
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_super_compact_pairwise_, X0, times, args = (tau, gamma, N, k_ave, ksquare_ave, kcube_ave))
    I, SS, SI, II = X.T
    S = N-I
    if return_full_data:
        return times, S, I, SS, SI, II
    else:
        return times, S, I



def SIR_super_compact_pairwise(SS0, SI0, R0, N, tau, gamma, psihat, psihatPrime, psihatDPrime, tmin = 0, tmax = 100, tcount = 1001, return_full_data = False):
    '''Encodes system (5.22) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{theta} = -tau [SI]/N*psihat(theta)
    [dot{SS}] = -2 tau [SS] [SI] Q
    [dot{SI}] = -gamma[SI] + tau([SS]-[SI])[SI]Q - tau*[SI]
    [dot{R}] = gamma*[I]
    [S] = N psihat(theta)
    [I] = N-[S]-[R]
    Q = psihat_xx(theta)/N(psihat_x(theta))^2


    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
            tells whether to just return times, S, I, R or all calculated data.
'''
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.array([1., SS0, SI0, R0])
    X = integrate.odeint(_dSIR_super_compact_pairwise_, X0, times, args = (tau, gamma, psihat, psihatPrime, psihatDPrime, N))
    theta, SS, SI, R = X.T
    S = N*psihat(theta)
    I = N-S-R
    if return_full_data:
	return times, S, I, R, SS, SI
    else:
	return times, S, I, R

def SIS_super_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]

    S0 = sum(Sk0)
    I0 = sum(Ik0)
    R0 = 0

    SX0 = scipy.dot(Sk0,ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0
    II0 = scipy.dot(Nk,ks)-SX0
    Pk = get_Pk(G)
    Pks = scipy.array([Pk.get(k,0) for k in ks])
    k_ave = scipy.dot(Pks, ks)
    ksquare_ave = scipy.dot(Pks, ks*ks)
    kcube_ave = scipy.dot(Pks, ks*ks*ks)
    
    return SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, k_ave, ksquare_ave, kcube_ave, tmin=tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)

def SIR_super_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]

    N = G.order()
    R0 = 0

    SX0 = Sk0.dot(ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0

    Pk = get_Pk(G)
    def psihat(x): #probably faster if vectorized, but need to be careful with broadcasting...
        return (1-rho)*sum(Pk[k]*x**k for k in Pk)
    def psihatPrime(x):
        return (1-rho)*sum(k*Pk[k]*x**(k-1) for k in Pk)
    def psihatDPrime(x):
        return (1-rho)*sum(k*(k-1)*Pk[k]*x**(k-2) for k in Pk)

    return  SIR_super_compact_pairwise(SS0, SI0, R0, N, tau, gamma, psihat, psihatPrime, psihatDPrime, tmin = tmin, tmax = tmax, tcount = tcount, return_full_data = return_full_data)






########     EFFECTIVE DEGREE


def _dSIS_effective_degree_(X, t, original_shape, tau, gamma):
    '''
    \dot{S}_{s,i} = - tau i S_{s,i} + gamma*I_{s,i}
                    + gamma((i+1)S_{s-1,i+1}-iS_{s,i})
                    + tau[ISS]((s+1)S_{s+1,i-1} - sS_{s,i})/[SS]

    \dot{I}_{s,i} = tau i S_{s,i} - gamma I_{s,i}
                    + gamma((i+1)I_{s-1,i+1} - iI_{s,i})
                    + tau([ISI]/[SI] + 1)((s+1)I_{s+1,i-1} - sI_{s,i})
    S = sum S_{s,i}
    I = sum I_{s,i}
    '''
    #note that dSIR_effective_degree has some commands for vectorizing tentatively coded (and commented out)
    ksq = original_shape[0]*original_shape[1]
    Ssi = X[:ksq]
    Isi = X[ksq:]
    Ssi.shape = original_shape
    Isi.shape = original_shape

    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    ISI = sum([sum([i*(i-1)*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    SI = sum([sum([i*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])

    g1 = scipy.zeros(original_shape)
    g2 = scipy.zeros(original_shape)
    t1 = scipy.zeros(original_shape)
    t2 = scipy.zeros(original_shape)
    
    dSsi = scipy.zeros(original_shape)
    dIsi = scipy.zeros(original_shape)
    for s in xrange(original_shape[0]):
        for i in xrange(original_shape[1]):
            if s==0 or i+1 == original_shape[1]:
                Ssm1ip1 = 0
                Ism1ip1 = 0
            else:
                Ssm1ip1 = Ssi[s-1,i+1]
                Ism1ip1 = Isi[s-1,i+1]
            if i==0 or s+1 == original_shape[0]:
                Ssp1im1 = 0
                Isp1im1 = 0
            else:
                Ssp1im1 = Ssi[s+1,i-1]
                Isp1im1 = Isi[s+1,i-1]
            
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*Isi[s,i] + gamma*((i+1)*Ssm1ip1 - i*Ssi[s,i]) + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
            dIsi[s,i] =  tau*i*Ssi[s,i] - gamma*Isi[s,i] + gamma*((i+1)*Ism1ip1 - i*Isi[s,i]) + tau*(ISI/SI + 1)*((s+1)*Isp1im1 - s*Isi[s,i])# 

    dSsi.shape = (original_shape[0]*original_shape[1])
    dIsi.shape = (original_shape[0]*original_shape[1])

    dX = scipy.concatenate((dSsi, dIsi), axis=0)
    return dX


def _dSIR_effective_degree_(X, t, N, original_shape, tau, gamma):
    R = X[-1]
    Ssi = X[:-1]
    Ssi.shape=original_shape
    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    
    #commenting out commands for vectorizing this.  I should do this eventually, but not now.  Apply to SIS version as well.
    #ultimately I think right option is to pad with zeros and then cut out appropriate section.
    #ivec = scipy.array(range(original_shape[1]))
    #ivec.shape = (1,original_shape[1])
    #svec = scipy.array(range(original_shape[0]))
    #svec.shape = (original_shape[0],1)
    #
    #Ssip1 = Ssi[:,1:]
    #scipy.pad(Ssip1, pad_width=((0,0),(0,1)), mode = 'constant', constant_values=0)
    #Ssp1im1 = Ssi[1:,:-1]
    #scipy.pad(Ssp1im1, pad_width=((0,1),(1,0)), mode = 'constant', constant_values=0)
    #dSsi = - tau* ivec*Ssi + gamma*((ivec+1)*Ssip1 - i*Ssi) + tau *ISS*((svec+1)*Sp1im1 - svec*Ssi)/SS

    dSsi = scipy.zeros(original_shape)
    for s in xrange(original_shape[0]):
        for i in xrange(original_shape[1]):
            if i+1 == original_shape[1]:
                Ssip1 = 0
            else:
                Ssip1 = Ssi[s,i+1]
            if s+1 == original_shape[0] or i == 0:
                Ssp1im1=0
            else:
                Ssp1im1 = Ssi[s+1,i-1]            
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*((i+1)*Ssip1 - i*Ssi[s,i]) + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
    S = Ssi.sum() 
    I = N-S-R
    dR = gamma*I

    dSsi.shape = (original_shape[0]*original_shape[1])
    dX = scipy.concatenate((dSsi, [dR]), axis=0)
    return dX



def SIS_effective_degree(Ssi0, Isi0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.36) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    INPUTS
    ---------
    Ssi0 and Isi0 : (square) numpy 2D arrays of same shape.
                      Entries are initial number susceptible or infected with
		      given initial number of susceptible/infected neighbors.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       if True, return times, S, I, Ssi, Isi
                       if False, return times, S, I
                       
                       tells whether to just return times, S, I, R or all calculated data.
    '''
    times = scipy.linspace(tmin,tmax,tcount) 
    original_shape = Ssi0.shape
    ksq = original_shape[0]*original_shape[1]
    Ssi0.shape = (1,ksq)
    Isi0.shape = (1,ksq)
    
    X0= scipy.concatenate((Ssi0[0], Isi0[0]), axis=0)
    X = integrate.odeint(_dSIS_effective_degree_, X0, times, args = (original_shape, tau, gamma))
    Ssi = X.T[0:ksq]
    Isi = X.T[ksq:]
    S = Ssi.sum(axis=0)
    I = Isi.sum(axis=0)
    if return_full_data:
        Ssi.shape = (original_shape[0],original_shape[1],tcount)
        Isi.shape = (original_shape[0],original_shape[1],tcount)
        return times, S, I, Ssi, Isi
    else:
        return times, S, I

def SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.38) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_{s,i} = - tau i S_{s,i}  + gamma((i+1)S_{s,i+1} - i S_{s,i})
                    + tau [ISS]((s+1)S_{s+1,i-1} - sS_{s,i})/[SS]
    \dot{R} = gamma I
    S = \sum_{s,i} S_{s,i}
    I = N-S-R

    INPUTS
    ---------
    S_si0 : (square) numpy 2-D array
            S_{s,i} at time 0
    I0 : number
         number of infected individuals at time 0
    R0 : number
         number of recovered individuals at time 0
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.
    '''
    times = scipy.linspace(tmin,tmax, tcount)
    N = S_si0.sum()+I0+R0
    original_shape = S_si0.shape
    S_si0.shape = (original_shape[0]*original_shape[1]) #note this makes it array([[values...]])
    R0=scipy.array([R0])
    R0.shape=(1)
    X0 = scipy.concatenate((S_si0, R0), axis=0)
    X = integrate.odeint(_dSIR_effective_degree_, X0, times, args = (N, original_shape, tau, gamma))

    R = X.T[-1]
    S_si = X.T[:-1]
    S = S_si.sum(axis=0)
    I = N - R - S
    if return_full_data:
        S_si.shape = (original_shape[0], original_shape[1], tcount)
	return  times, S, I, R, S_si
    else:
        return times, S, I, R


def SIS_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    S_si0 = scipy.zeros((len(Sk0),len(Sk0)))
    I_si0 = scipy.zeros((len(Sk0),len(Sk0)))

    for s in range(len(Sk0)):
        for i in range(len(Sk0)-s):
            S_si0[s,i] = Sk0[s+i]*scipy.special.binom(s+i,i)*(rho**i)*(1-rho)**s
            I_si0[s,i] = Ik0[s+i]*scipy.special.binom(s+i,i)*(rho**i)*(1-rho)**s

    return SIS_effective_degree(S_si0, I_si0, tau, gamma, tmin = tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)


def SIR_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    S_si0= scipy.zeros((len(Sk0),len(Sk0)))
    for s in range(len(Sk0)):
        for i in range(len(Sk0)-s):
            S_si0[s,i] = Sk0[s+i]*scipy.special.binom(s+i,i)*(rho**i)*(1-rho)**s
    I0 = sum(Ik0)
    R0 = sum(Rk0)
    return SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)




#######     COMPACT EFFECTIVE DEGREE


def SIS_compact_effective_degree(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    '''This model is identical to the SIS compact pairwise model, so it simply calls SIS_compact_pairwise()'''

    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, tmax, tcount, return_full_data)

def SIS_compact_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    return SIS_compact_pairwise_from_graph(G, tau, gamma, rho = rho, tmin = tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)


def _dSIR_compact_effective_degree_(X, t, N, tau, gamma):
    Skappa = X[:-2]
    R, SI = X[-2:]
    I = N- R- Skappa.sum()
    kappas = scipy.arange(len(Skappa))
    effectiveI = float(SI) /Skappa.dot(kappas)
    dSkappa = effectiveI*(-(tau+gamma)*kappas*Skappa + gamma*shift(kappas*Skappa,-1))
    dSI = -(tau+gamma)*SI + tau*(effectiveI-2*effectiveI**2)*sum(kappas*(kappas-1)*Skappa)

    dR = gamma*I
    dX = scipy.concatenate((dSkappa, [dR, dSI]), axis=0) 
    return dX
    
def SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, tmin=0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.43) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_kappa = <I> [-(tau+gamma) kappa S_kappa + gamma(kappa+1)S_{kappa+1}
    [\dot{SI}] = -(tau+gamma)[SI] + tau(<I> - 2 <I>^2) sum_{kappa} kappa(kappa-1) S_kappa
    \dot{R} = gamma I
    <I> = [SI]/sum_kappa kappa S_kappa
    S = sum_kappa S_kappa
    I = N - S - R

    INPUTS
    ---------
    Skappa0 : scipy array
              from S_0(0) up to S_kappamax(0) of number susceptible with each effective degree
              Skappa = number of nodes that are susceptible and have kappa non-recovered neighbors
    I0 : number
         number of infected individuals at time 0
    R0 : number
         initial number recovered
    SI0 : number
          initial number of SI edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or all calculated data.
    
    RETURNS (if return_full_data==False):
    --------
    times : scipy.array of times
    S : scipy.array of number susceptible
    I : scipy.array of number infected
    R : scipy.array of number recovered

    if return_full_data==True
    --------------------------
    times : as before
    Skappa : array of arrays, each subarray gives particular S_kappa
    I : number infected
    R : number recovered
    SI : number of SI edges
    '''

    N = Skappa0.sum()+I0+R0
    X0= scipy.concatenate((Skappa0,[R0,SI0]), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_compact_effective_degree_, X0, times, args = (N, tau, gamma))
    Skappa = X.T[:-2]
    S = Skappa.sum(axis=0)
    R, SI = X.T[-2:]
    I = N - S -R
    if return_full_data:
        return times, S, I, R, Skappa, SI
    else:
	return times, S, I, R

def SIR_compact_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    Skappa0 = Sk0
    I0 = sum(Nk)*rho
    R0 = 0
    SI0 = sum([k*Skappa0[k]*rho for k in range(len(Sk0))])
    return SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, tmin=tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)






#######################################
#    EBCM and other related results   #
#######################################

def Epi_Prob_discrete(Pk, p, number_its = 100):
    '''Encodes System (6.2) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    INPUTS
    ---------
    Pk : scipy array Pk[k] is probability a node has degree k.

    p : transmission probability

    number_its : number of iterations before assumed converged.
                 default value is 100

    RETURNS
    ----------
    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        return (ks*Pk).dot(x**(ks-1))
    
    alpha = 1-p
    k_ave = psiPrime(1)
    for counter in range(number_its):
        alpha = 1-p +p *psiPrime(alpha)/k_ave
    return 1- psi(alpha)

def Epi_Prob_cts_time(Pk, tau, gamma, umin=0, umax = 10, ucount = 1001, number_its = 100):
    '''Encodes System (6.3) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The equations are rescaled by setting $u=\gamma T$.  Then it becomes

    P = 1- \int_0^\infty \psi(\alpha(u/\gamma)) e^{-u} du
    alpha_d(u/\gamma) = 1- p(u/\gamma)
                        + p(u/\gamma)\int_0^\infty (\psiPrime(\alpha(\hat{u}/\gamma))/<K>)e^{-u}du

    where p(u/\gamma) = 1 - e^{-\tau u/\gamma}

    Define \hat{p}(u) = p(u/\gamma), and \hat{\alpha}(u) = \alpha(u/\gamma)
    and then drop hats to get

    P = 1-\int_0^\infty \psi(\alpha(u)) e^{-u} du
    \alpha(u) = 1-p(u) + p(u) \int_0^\infty (\psiPrime(\alpha(u))/<K>)e^{-u} du

    with initial guess \alpha_1(u) = e^{-\tau u/\gamma} and p(u) = 1-e^{-\tau u/\gamma}
    
    INPUTS
    ---------
    Pk : scipy array Pk[k] is probability a node has degree k.

    tau : transmission rate

    gamma : recovery rate

    umin : minimal value of \gamma T used in calculation
    umax : maximum value of \gamma T used in calculation
    ucount : number of points taken for integral.
           So this integrates from umin to umax using simple Riemann sum.

    number_its : number of iterations before assumed converged.
                 default value is 100

    RETURNS
    ----------
    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        ks*Pk
        x**(ks-1)
        return (ks*Pk).dot(x**(ks-1))
    us = scipy.linspace(umin, umax, ucount) 
    alpha = scipy.e**(-tau*us/gamma)  #initial guess for alpha(u)
    p = 1- scipy.e**(-tau*us/gamma)    #initial guess for p(u)
    exp_neg_u = scipy.e**(-us)        #e^{-u}
    for counter in xrange(number_its):
        alpha = 1 - p + p* (psiPrime(alpha)/kave).dot(exp_neg_u)/(ucount-1.)
    return 1 - psi(alpha).dot(exp_neg_u)/(ucount - 1)


def Epi_Prob_non_Markovian(Pk, tau, gamma, Pxidxi, po, umin=0, umax = 10, ucount = 1001, number_its = 100):
    '''Encodes system (6.5) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    INPUTS
    ---------
    Pk : scipy array Pk[k] is probability a node has degree k.

    tau : transmission rate

    gamma : recovery rate

    Pxidxi : a dict.  Returns P(xi)dxi for user-selected xi.  The algorithm will replace the integral with
           \sum_{xi \in Pxidxi.keys()} \psi(\alpha(xi)) Pxidxi(xi)

    po : a function.
         returns p_o(xi), the probability a node will transmit to a random neighbor given xi.
    umin : minimal value of \gamma T used in calculation
    umax : maximum value of \gamma T used in calculation
    ucount : number of points taken for integral.
           So this integrates from umin to umax using simple Riemann sum.

    number_its : number of iterations before assumed converged.
                 default value is 100

    RETURNS
    ----------
    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        return (ks*Pk).dot(x**(ks-1))
    xis = Pxidxi.keys()
    alpha = {xi: 1-po(xi) for xi in xis}
    for counter in xrange(number_its):
        newalpha = {}
        for xi in xis:
            newalpha[xi] = 1 - po(xi)  + po(xi)*sum(psiPrime(alpha[xihat])*Pxidxi(xihat) for xihat in xis)/kave
        alpha = newalpha
    return 1 - sum(psi(alpha[xi])*Pxidxi[xi] for xi in xis)

def Attack_rate_discrete(Pk, p, number_its=100, rho = None, Sk0=None, phiS0=None, phiR0=0):
    '''Encodes systems (6.6) and (6.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    To use system (6.6), leave rho and Sk0 as None.

    INPUTS
    ------
    Pk : dict
        Pk[k] is the probability a randomly selected node has degree k.
    tau : number
        per-edge transmission rate.
    gamma : number
        per-node recovery rate
    number_its : The solution is found iteratively, so this determines the number of iterations.
    rho : Number (default None)
    Sk0 : dict (default None)
          only one of rho and Sk0 can be defined.  The other (or both) should remain None.
          rho gives the fraction of nodes randomly infected.
          Sk0 is a dict such that Sk0[k] is the probability that a degree k node is susceptible at start.
    phiS0 : number (default None)
          Should only be used if Sk0 is not None.  If it is None, then assumes that initial introduction
          is randomly introduced
    phiR0 : number (default 0)
          As with phiS0, only used if Sk0 is not None.
    OUTPUT
    ------
    A : number
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoNError("at most one of rho and Sk0 can be defined")
    if Sk0 is None:
        if rho is None or rho == 0:
            return Epi_Prob_discrete(Pk, p, number_its)
        else:
            Sk0 = {k: (1-rho) for k in Pk.keys()}
    def psihat(x):
        return sum(Pk[k]*Sk0[k]*x**k for k in Pk.keys())
    def psihatPrime(x):
        return sum(k*Pk[k]*Sk0[k]*x**(k-1) for k in Pk.keys())

    if phiS0 == None:
        phiS0 = psihatPrime(1)/sum(k*Pk[k] for k in Pk.keys())
    if phiR0 == None:
        phiR0 = 0

    theta  = 1
    for counter in xrange(number_its):
        theta = 1-p + p*(phiR0 +  phiS0*psihatPrime(theta)/psihatPrime(1))
    return 1 - psihat(theta)

def Attack_rate_discrete_from_graph(G, p, number_its = 100, rho = None, Sk0 = None):
    Pk = get_Pks(G)
    return Attack_rate_discrete_from_graph(Pk, p, number_its = number_its, rho = rho, Sk0 = Sk0)

def Attack_rate_cts_time(Pk, tau, gamma, number_its =100, rho = None, Sk0 = None, phiS0=None, phiR0=0):
    #tested in test_SIR_final_sizes
    '''Encodes system (6.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    
    This system predicts the fraction of nodes infected if an epidemic occurs in a
    Configuration Model network assuming a continuous-time Markovian SIR disease.  This gives the
    limit of the attack rate of epidemics as the initial fraction infected approaches 0.

    If we look for the limit of a nonzero initial fraction infected, we introduce rho or Sk0

    INPUTS
    ------
    Pk : dict
        Pk[k] is the probability a randomly selected node has degree k.
    tau : number
        per-edge transmission rate.
    gamma : number
        per-node recovery rate
    number_its : The solution is found iteratively, so this determines the number of iterations.
    rho : Number (default None)
    Sk0 : dict (default None)
          only one of rho and Sk0 can be defined.  The other (or both) should remain None.
          rho gives the fraction of nodes randomly infected.
          Sk0 is a dict such that Sk0[k] is the probability that a degree k node is susceptible at start.

    OUTPUT
    ------
    A : number
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoNError("at most one of rho and Sk0 can be defined")
    if Sk0 is None:
        if rho is None:
            rho = 0
        Sk0 = {k: (1-rho) for k in Pk.keys()}
    def psihat(x):
        return sum(Pk[k]*Sk0[k]*x**k for k in Pk.keys())
    def psihatPrime(x):
        return sum(k*Pk[k]*Sk0[k]*x**(k-1) for k in Pk.keys())

    if phiS0 == None:
        phiS0 = psihatPrime(1)/sum(k*Pk[k] for k in Pk.keys())
    if phiR0 == None:
        phiR0 = 0

    kave = sum(Pk[k]*k for k in Pk.keys())
    omega = gamma/(gamma+tau)
    for counter in range(number_its):
        print omega
        omega = gamma/(gamma+tau) + tau*phiS0*psihatPrime(omega)/(psihatPrime(1)*(gamma+tau)) + tau*phiR0/(gamma+tau)
    return 1 - psihat(omega)

def Attack_rate_cts_time_from_graph(G,  tau, gamma, number_its =100, rho=None, Sk0 = None):
    '''
    Given a graph, predicts the attack rate for Configuration Model networks with the given degree distribution.
    First calculates the degree distribution and then calls Attack_rate_cts_time
    '''
    Pk = get_Pk(G)
    return Attack_rate_cts_time(Pk, tau, gamma, number_its = number_its, rho=rho, Sk0=Sk0)
    
def Attack_rate_non_Markovian():
    '''Encodes system (6.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    '''
    pass

def EBCM_discrete(N, psihat, psihatPrime, p, phiS0, phiR0=0, R0=0, tmax = 100, return_full_data = False):
    #tested in test_basic_discrete_SIR_epidemic
    '''Encodes system (6.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    theta(t) = (1-p) + p(phi_R(0) + phi_S(0) psihatPrime(theta(t-1))/psihatPrime(1))
    R(t) = R(t-1) + I(t-1)
    S(t) = N psihat(theta(t))
    I(t) = N-S-R

    INPUTS
    ------
    N : number
        size of population
    psihat : function
        psihat(x) = \sum_k S(k,0) x^k
    psihatPrime : function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    p : number
        per edge transmission probability
    phiS0 : number
        initial proportion of edges (of susceptible nodes) connecting to susceptible nodes
    phiR0 : number
        initial proportion of edges (of susceptible nodes) connecting to recovered nodes
    R0 : number
        number of recovered nodes at time 0
    tmax : number
        maximum time
    return_full_data : boolean
        if False, return t, S, I, R and if True return t, S, I, R, and theta

    OUTPUTS
    -------
    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    times = [0]
    theta = [1]
    R = [R0]
    S = [N*psihat(1)]
    I = [N-S[-1]-R[-1]]

    for time in range(1,tmax+1):
        times.append(time)
        newtheta = (1-p) + p *(phiR0 + phiS0*psihatPrime(theta[-1])/psihatPrime(1))
        newR = R[-1]+I[-1]
        newS = N*psihat(newtheta)
        newI = N-newR-newS
        theta.append(newtheta)
        R.append(newR)
        S.append(newS)
        I.append(newI)
    if not return_full_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R)
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), scipy.array(R), scipy.array(theta)

def EBCM_discrete_uniform_introduction(N, psi, psiPrime, p, rho, tmax=100, return_full_data=False):
    #tested in test_basic_discrete_SIR_epidemic
    '''Handles the case that the disease is introduced uniformly as opposed to depending on degree.

    INPUTS
    ------
    N : Number
        number of nodes
    psi : function
        psi(x) = \sum P(k) x^k
    psiPrime : function
        psiPrime(x)=d psi(x)/dx = \sum kP(k) x^{k-1}
    p : number
        per edge transmission probability
    rho : number
        initial proportion of infected nodes
    tmax : number
        maximum time
    return_full-data : boolean
        if False, return t, S, I, R and if True return t, S, I, R, and theta

    OUTPUTS
    -------
    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM_discrete(N, psihat, psihatPrime, p, 1-rho, tmax=tmax, return_full_data=return_full_data)


def EBCM_discrete_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    #tested in test_basic_discrete_SIR_epidemic
    '''
    Takes a given graph, finds the degree distribution (from which it gets psi),
    assumes a constant proportion of the
    population is infected at time 0, and then uses the discrete EBCM model.
    
    INPUTS
    ------
    G : Networkx Graph
    p : number
        per edge transmission probability
    rho : number
        initial proportion of infected nodes
    tmax : number
        maximum time
    return_full-data : boolean
        if False, return t, S, I, R and if True return t, S, I, R, and theta

    OUTPUTS
    -------
    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
 '''
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    def psi(x):
        return sum(Pk[k]*x**k for k in Pk)
    def psiPrime(x):
        return sum(k*Pk[k]*x**(k-1) for k in Pk)
    return EBCM_discrete_uniform_introduction(G.order(), psi, psiPrime, p, rho, tmax=tmax, return_full_data=False)


def _dEBCM_(X, t, N, tau, gamma, psihat, psihatPrime, phiS0, phiR0):
    theta = X[0]
    R = X[1]
    
    dtheta = -tau*theta + tau*phiS0*psihatPrime(theta)/psihatPrime(1) + gamma*(1-theta) + tau*phiR0

    S = N*psihat(theta)
    I = N-S-R
    dR = gamma*I
    return scipy.array([dtheta, dR])
    
def EBCM(N, psihat, psihatPrime, tau, gamma, phiS0, phiR0=0, R0=0, tmin=0, tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (6.12) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    note : R0 is R(0), not the reproductive number


    INPUTS
    ------
    N : number
        size of population
    psihat : function
        psihat(x) = \sum_k S(k,0) x^k
    psihatPrime : function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    tau : number
        per edge transmission rate
    gamma : number
        per node recovery rate
    phiS0 : number
        initial proportion of edges (of susceptible nodes) connecting to susceptible nodes
    phiR0 : number
        initial proportion of edges (of susceptible nodes) connecting to recovered nodes
    R0 : number
        number of recovered nodes at time 0
    tmin : number
        start time
    tmax : number
        stop time
    tcount : integer
        number of distinct times to calculate
    return_full_data : boolean
        if False, return t, S, I, R and if True return t, S, I, R, and theta

    OUTPUTS
    -------
    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    times = scipy.linspace(tmin, tmax, tcount)
    X0 = scipy.array([1, R0])
    X = integrate.odeint(_dEBCM_, X0, times, args = (N, tau, gamma, psihat, psihatPrime, phiS0, phiR0))
    theta = X.T[0]
    R = X.T[1]
    S = N*psihat(theta)
    I = N-S-R
    if not return_full_data:
        return times, S, I, R
    else:
        return times, S, I, R, theta

def EBCM_uniform_introduction(N, psi, psiPrime, tau, gamma, rho, tmin=0, tmax=100, tcount=1001, return_full_data=False):
    r'''Handles the case that the disease is introduced uniformly as opposed to depending on degree.

    INPUTS
    ------
    N : number
        size of population
    psi : function
        psihat(x) = \sum_k S(k,0) x^k
    psiPrime : function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    tau : number
        per edge transmission rate
    gamma : number
        per node recovery rate
    rho : number
        initial proportion infected
    tmin : number
        start time
    tmax : number
        stop time
    tcount : integer
        number of distinct times to calculate
    return_full_data : boolean
        if False, return t, S, I, R and if True return t, S, I, R, and theta

    OUTPUTS
    -------
    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
'''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM(N, psihat, psihatPrime, tau, gamma, 1-rho, tmin=tmin, tmax=tmax, tcount=tcount, return_full_data=return_full_data)

def EBCM_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, tcount=1001, return_full_data=False):
    #tested in test_SIR_dynamics
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    def psi(x):
        return sum(Pk[k]*x**k for k in Pk)
    def psiPrime(x):
        return sum(k*Pk[k]*x**(k-1) for k in Pk)
    return EBCM_uniform_introduction(G.order(), psi, psiPrime, tau, gamma, rho, tmax=tmax, return_full_data=return_full_data)


'''
These are the systems I want to include based on their numbering in the book:

coded (3.7) SIS individual based
(3.30) SIR individual based
NOT coded (3.26) SIS pair based
(3.39) SIR pair based

chapter 4?

(5.13) SIS heterogeneous pairwise
(5.15) SIR heterogeneous pairwise
(5.18) SIS compact pairwise
(5.19) SIR compact pairwise
(5.20) SIS super compact pairwise
(5.22) SIR super compact pairwise
(5.36) SIS effective degree
(5.38) SIR effective degree
(5.43) SIR compact effective degree
(5.44) SIS compact effective degree = SIS compact pairwise

(6.2) Epidemic probability discrete time
(6.3) Epidemic probability continuous time
(6.5) Epidemic probability non-Markovian
(6.6) Epidemic size discrete time
(6.7) Epidemic size continuous time
(6.8) Epidemic size non-Markovian
(6.10) Epidemic size discrete time (large IC)
(6.11) Discrete-time EBCM model
(6.12) Continuous time EBCM model

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

add models that take in graph, measure degree distribution and run EBCM
similarly for EBCM with neighbor degrees (see barabasi_SIR.py)

consider explicitly defining toast graph etc.

SIR_pair_based --- which of 2 versions to keep, and comments need to explain it a bit better.
'''
