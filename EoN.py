import networkx as nx
from collections import defaultdict
import random
import heapq
import random

r'''
EoN

"Epidemics on Networks"
======

    EoN is a Python package for the simulation of epidemics on networks.

    We assume that networks are created using the NetworkX package.

    The algorithms are based on the book:
            Mathematics of network epidemics: from exact to approximate models
            - Kiss, Miller & Simon
        Please cite the book if using these algorithms


This is a preliminary version:
  - The algorithms have not been tested in Python 3 (tested only in 2.7)
  - The figure references are based on the current draft of the book, which may still be edited, so figure numbers may change.
  - Additional algorithms may be added

'''

def discrete_SIR_epidemic(G, test_transmission, parameters=[], initial_infecteds=None):
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Return details of epidemic curve from a discrete time simulation.
    It assumes that individuals are infected for exactly one unit of time
    and then recover with immunity.  This is defined to handle a user-defined function
         test_transmission(node1,node2,*parameters)
    which determines whether transmission occurs.  So more elaborate rules can be created
    as desired by the user.

    Parameters:
    ----------
    G: NetworkX Graph (or some other structure which quacks like a NetworkX Graph)
        The network on which the epidemic will be simulated.
        
    test_transmission: function(u,v,*parameters)
        (see below for parameters definition)
        A function that determines whether u transmits to v.
        It returns True if transmission happens and False otherwise.
        This function can be user-defined.
          If parameters is None, then it is called from this function as
            test_transmission(u,v)
          otherwise it is called as
            test_transmission(u,v,*parameters)

    parameters: a list or tuple
        The parameters of test_transmission coming after the nodes.

    initial_infecteds: list (could be any iterable)
        The initially infected nodes.  



    Returns
    ------------
    the lists: t, S, I, R
    these lists give all the times observed and the number in each state at each time.
    '''
    
    if initial_infecteds is None:
        infecteds = [random.choice(G.nodes())]
    else:
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
                if susceptible[v] and test_transmission(u, v, *parameters):
                    new_infecteds.append(v)
                    susceptible[v] = False
        infecteds = new_infecteds
        R.append(R[-1]+I[-1])
        I.append(len(infecteds))
        S.append(S[-1]-I[-1])
        t.append(t[-1]+1)
    return t, S, I, R


def simple_test_transmission(u, v, p):
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm.

    test_transmission function for basic_discrete_SIR_epidemic.

    This handles the simple case where transmission occurs with probability p.
    '''
    #p = parameterlist[0]
    return random.random()<p

def basic_discrete_SIR_epidemic(G, p, initial_infecteds=None):
    r'''From figure 6.8 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm.

    The simple case of all nodes transmitting with probability p independently to each neighbor and then recovering.

    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability
    initial_infecteds : list (could be any iterable)
        The initially infected nodes. 

    Returns
    -------
    the lists: t, S, I, R
    these lists give all the times observed and the number in each state at each time.
'''

    return discrete_SIR_epidemic(G, simple_test_transmission, [p], initial_infecteds)


def percolate_network(G, p):
    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Performs bond percolation on the network G.

    Parameters
    ----------
    G : NetworkX Graph
    p : number
        probability of keeping edge

    Returns
    -------
    H : NetworkX Graph
        A network with same nodes as G, but with each edge retained independently with probability p.
'''

    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for edge in G.edges_iter():
        if random.random()<p:
            H.add_edge(edge)
    return H

def edge_exists(u, v, H):
    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.

    Tests wether H has an edge from u to v.
'''
    return H.has_edge(u,v)

def percolation_based_discrete_SIR_epidemic(G, p, initial_infecteds):
    r'''From figure 6.10 of Kiss, Miller, & Simon.  Please cite the book
    if using this algorithm.


    
    The simple case of all nodes transmitting with probability p independently to each neighbor
    and then recovering, but using a percolation-based approach.  See basic_discrete_SIR_epidemic
    which should produce equivalent outputs.  That algorithm will be faster than this one.  The
    value of this function is that by performing many simulations we can see that the outputs of
    the two are equivalent.  This algorithm leads to a better understanding of the theory.


    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability
    initial_infecteds : list (could be any iterable)
        The initially infected nodes. 

    Returns
    -------
    the lists: t, S, I, R
    these lists give all the times observed and the number in each state at each time.
'''

    H = percolate_network(G, p)
    t, I = discrete_SIR_epidemic(H, edge_exists, [H], initial_infecteds)


def estimate_SIR_prob_size(G, p):
    r'''From figure 6.12 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    Provies an estimate of epidemic probability and size assuming a fixed transmission
    probability p.  The estimate is found by performing bond percolation and then finding
    the largest connected component in the remaining network.  This assumes that there is
    a single giant component above threshold.  It will not be an appropriate measure if the
    network is made up of several densely connected components with very weak connections between
    these components.

    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    p : number
        transmission probability

    Returns
    -------
    P, A : estimates of the probability and proportion infected in epidemics
        (the two are equal, but each given for consistency with estimate_directed_prob_size)
    '''
    H = percolate_network(G, p)
    size = max((len(CC) for CC in nx.connected_components(H)))
    returnval = float(size)/G.order()
    return returnval, returnval


def directed_percolate_network(G, tau, gamma):
    r'''From figure 6.13 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.  This adds node and edge attributes
    which are not at present in the figure.  This option is discussed
    in the text.


    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    tau : number
        transmission rate
    gamma : number
        recovery rate

    Returns
    -------
    H : networkx DiGraph  (directed graph)
        a u->v edge exists in H if u would transmit to v if ever infected.
        The edge has a time attribute (time_to_infect) which gives the delay
            from infection of u until transmission occurs.
        each node has a time attribute (duration) which gives the duration of
            u's infectious period.
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
                
def out_component(G, source):
    '''rather than following the pseudocode in figure 6.15 of Kiss, Miller & Simon, this uses a built in networkx command.  I plan to improve this algorithm.

    finds the set of nodes (including source) which are reachable from nodes in source.

    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    source : either a node or an iterable of nodes (set, list, tuple)
        The nodes from which the infections start.

    Returns
    -------
    reachable_nodes : set
        the set of nodes reachable from source (including source).
    '''
    try:
        #testing whether this is an iterable
        iterator = iter(source)
    except TypeError:
        #It's not an iterable.  It "must" be a node.
        if G.has_node(source):
            source_nodes = set([source])
    else:
        #it's an iterable.  
        source_nodes = set(source)
    reachable_nodes = set([])
    for node in source_nodes:
        reachable_nodes = reachable_nodes.union(set(nx.descendants(G, node)))
    return reachable_nodes

def in_component(G, target):
    r'''creates the in_component by basically reversing out_component.

    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    target : a target node
        The node whose infection we are interested in.

        In principle target could be an iterable, but in this case we would be finding those possible
        sources whose infection leads to infection of at least one target, not all.

    Returns
    -------
    source_nodes : set
        the set of nodes (including target) from which target is reachable
    '''
    try:
        #testing whether this is an iterable
        iterator = iter(target)
    except TypeError:
        #It's not an iterable.  It "must" be a node.
        if G.has_node(target):
            target_nodes = set([target])
    else:
        #it's an iterable.  
        target_nodes = set(target)
    source_nodes = set([])
    for node in target_nodes:
        source_nodes = source_nodes.union(set(nx.ancestors(G, node)))
    return source_nodes


def get_infected_nodes(G, initial_infecteds, tau, gamma):
    r'''From figure 6.15 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm

    Finds all eventually infected nodes in a simulation, assuming that the intial infecteds are as given
    and transmission occurs with rate tau and recovery with rate gamma.  Uses a percolation-based
    approach.

    This code has similar run-time whether an epidemic occurs or not.  
    
    Parameters
    ----------
    G : NetworkX Graph
        The network the disease will transmit through.
    initial_infecteds : either a node or an iterable of nodes (set, list, tuple)
        The nodes from which the infections start.
    tau : number
        transmission rate
    gamma : number
        recovery rate

    Returns
    -------
    infected_nodes : set
        the set of nodes infected eventually in a simulation.
    '''    
    H = directed_percolate_network(G, tau, gamma)
    infected_nodes = out_component(G, initial_infecteds)
    return infected_nodes


def estimate_directed_prob_size(H):
    r'''From figure 6.17 of Kiss, Miller, & Simon.  Please cite the book if using this algorithm

    Parameters
    ----------
    H: directed graph (assumed to be from directed percolation on previous graph G)

    Returns:
    PE, AR  :  numbers
        Estimates of epidemic probability and attack rate found by finding largest strongly
        connected component and finding in/out components.
    '''
    Hscc = max((nx.strongly_connected_components(H)), key = len)
    u = random.choice(Hscc)
    inC = in_component(H, u) #includes both H_{IN} and H_{SCC}
    outC = out_component(H, u) #includes both H_{OUT} and H_{SCC}
    N=float(H.order())
    PE = len(inC)/N
    AR = len(outC)/N
    return PE, AR
    
        
def nonMarkov_directed_percolate_network(G, xi, zeta, transmission):
    r'''From figure 6.18 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    xi and zeta are dictionaries where xi[u] and zeta[v] are enough to determine
    a u-v transmission.

    transmissision is a user-defined function taking xi[u] and zeta[v] and
    returning True if a transmission would occur

    Parameters
    ----------
    G : NetworkX Graph
        The input graph

    xi : Dict
        xi[u] gives all necessary information to determine what u's infectiousness is
    zeta : Dict
        zeta[v] gives everything needed about v's susceptibility

    transmission : user-defined function
        transmission(xi[u], zeta[v]) determines whether u transmits to v.
'''
    H = nx.DiGraph()
    for u in G.nodes():
        H.add_node(u)
        for v in G.neighbors(u):
            if transmission(xi[u],zeta[v]):
                H.add_edge(u,v)
    return H

#def process_recovery():

    #class Event(object):
    #    def __init__(self, node, time, action):
    #        self.node = node
    #        self.time = time
    #        self.action = action

def process_trans_SIR(G, node, time, tau, gamma, times, S, I, R, Q, status, rec_time, pred_inf_time, tmax):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    Parameters
    ----------
    G : NetworkX Graph

    node : node
        The node that is receiving transmission
    time : number
        current time
    tau : number
        transmission rate (from node)
    gamma : number
        recovery rate (of node)
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
    tmax : max time allowed

    Returns
    -------
    nothing returned

    Modifies
    --------
    times
    S
    I
    R
    Q
    '''
    
    status[node] = 'I'
    times.append(time)
    S.append(S[-1]-1) #one less susceptible
    I.append(I[-1]+1) #one more infected
    R.append(R[-1])   #no change to recovered
    rec_time[node] = time + random.expovariate(gamma)
    if rec_time[node] < tmax:
        newevent = {'node':node, 'action':'recover'}
        heapq.heappush(Q,(rec_time[node],newevent))
    for v in G.neighbors(node):
        find_trans_SIR(Q,time, tau, node, v, status, rec_time, pred_inf_time, tmax)

def find_trans_SIR(Q, t, tau, source, target, status, rec_time, pred_inf_time, tmax):
    r'''From figure 6.18 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    determines if a transmission from source to target sill occur and if so puts into Q

    Parameters
    Q
    t
    tau
    source
    target
    status
    rec_time
    pred_inf_time
    tmax
    '''
    if status[target] == 'S':
        delay = random.expovariate(tau)
        inf_time = t + delay
        if inf_time< rec_time[source] and inf_time < pred_inf_time[target] and inf_time<tmax:
            event = {'node': target, 'action': 'transmit'}
            heapq.heappush(Q, (inf_time, event))
            pred_inf_time[target] = inf_time

    
    
    
def process_recovery_SIR(node, t, times, S, I, R, status):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    '''
    times.append(t)
    S.append(S[-1])   #no change to number susceptible
    I.append(I[-1]-1) #one less infected
    R.append(R[-1]+1) #one more recovered
    status[node] = 'R'
    
    
def fast_SIR(G, tau, gamma, source=[1], tmax=1000):
    r'''From figure A.2 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

fast SIR simulation assuming exponentially distributed infection and recovery times based on figure A.1 of the appendix of Kiss, Miller & Simon.  Please cite the book if you use this algorithm.'''

    '''
    inputs:
    G: a graph
    tau: transmission rate per edge
    gamma: recovery rate per node
    source: list containing the initial infected nodes
    tmax: maximum time after which simulation will stop.
          Use float('Inf') to set to infinity
    '''
    
    times, S, I, R= ([0], [G.order()], [0], [0])
    Q = []#an empty heap

    status = defaultdict(lambda: 'S') #node status defaults to 'S'
    rec_time = defaultdict(lambda: -1) #node recovery time defaults to -1
    pred_inf_time = defaultdict(lambda: float('Inf')) #infection time defaults to \infty  --- this could be set to tmax, probably with a slight improvement to performance.
    for u in source:
        newevent = {'node':u, 'action':'transmit'}
        pred_inf_time[u] = 0
        Q.append((0,newevent))#okay to append rather than heappush since all are at same time

    while Q:
        time, event = heapq.heappop(Q)
        if event['action'] == 'transmit':
            if status[event['node']] == 'S':
                process_trans_SIR(G, event['node'], time, tau, gamma, times, S, I, R, Q, status, rec_time, pred_inf_time)
        else:
            process_recovery_SIR(event['node'], time, times, S, I, R, status)
    return times, S, I, R
            

    


def fast_SIS(G, tau, gamma, source=[1], tmax=1000):
    r'''From figure A.4 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    '''
    ''' fast SIS simulation assuming exponentially distributed infection and recovery times based on figure A.3 of the appendix of Kiss, Miller & Simon.  Please cite the book if you use this algorithm.

    
    inputs:
    G: a graph
    tau: transmission rate per edge
    gamma: recovery rate per node
    source: list containing the initial infected nodes
    tmax: maximum time after which simulation will stop.
    '''
    
    pass #I need to write up this algorithm cleanly still.
