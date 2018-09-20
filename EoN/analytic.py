# -*- coding: utf-8 -*-
from scipy import integrate
from scipy.ndimage.interpolation import shift
import scipy
import networkx as nx
import EoN
from collections import defaultdict, Counter

#######################
#                     #
#   Auxiliary stuff   #
#                     #
#######################ƒ



def _my_odeint_(dfunc, V0, times, args=()):
    r'''For some of the systems odeint will switch to the BDF solver.
    In large enough systems, it then gets stuck trying to estimate the 
    Jacobian.  This will even cause segmentation faults.

    This routine has identical inputs to integrate.odeint, but relies on 
    integrate.ode.  It avoids BDF.

    In particular, this seems to be important for SIS heterogeneous 
    pairwise where the number of equations is very large.  I have found 
    that near equilibrium, this often is interpreted as being a stiff 
    system and it switches to bdf, which requires calculating a 
    Jacobian.  In some systems this is impractically large.
    
    See this question: 
        http://stackoverflow.com/q/40317096/2966723,
    with the answer by 
        Phillip: http://stackoverflow.com/users/1881610/phillip
        
    INPUT and OUTPUT are as integrate.odeint
    
    :Arguments: 
    Identical to integrate.odeint
        
    :Returns: 
    Identical to integrate.odeint
    '''

    r = integrate.ode(lambda t, X: dfunc(X, t, *args))
    r.set_integrator('vode', method='adams')
    r.set_initial_value(V0,times[0])
    V=[V0]
    for time in times[1:]:
        V.append(r.integrate(time))
    V = scipy.array(V)
    return V

def _initialize_node_status_(G, initial_infecteds, initial_recovereds = None):
    if initial_recovereds is None:
        initial_recovereds = []
    intersection = set(initial_infecteds).intersection(set(initial_recovereds))
    if  intersection:
        raise EoN.EoNError("{} are in both initial_infecteds and initial_recovereds".format(intersection))
        
    status = defaultdict(lambda : 'S')
    for node in initial_infecteds:
        if not G.has_node(node):
            raise EoN.EoNError("{} not in G".format(node))
        status[node] = 'I'
    for node in initial_recovereds:
        if not G.has_node(node):
            raise EoN.EoNError("{} not in G".format(node))
        status[node] = 'R'
    return status

def _count_edge_types_(G, initial_infecteds, initial_recovereds = None, SIR=True):
    status = _initialize_node_status_(G, initial_infecteds, initial_recovereds)
    SS0 = 0
    SI0 = 0
    II0 = 0
    for u,v in G.edges():
        if status[u] == 'S':
            if status[v] == 'S':
                SS0 += 2
            elif status[v] == 'I':
                SI0 += 1
        elif status[u]=='I':
            if status[v] == 'S':
                SI0 += 1
            elif status[v] == 'I':
                II0 += 2
    if SIR:
        return SS0, SI0
    else:
        return SS0, SI0, II0
        
#################################
#                               #
#   Degree Distribution stuff   #
#                               #
#################################


def _get_Nk_and_IC_as_arrays_(G, initial_infecteds = None, initial_recovereds = None, rho=None, SIR=True):
    r'''
    Given the graph and initial proportion infected this finds the 
    initial conditions and number of nodes of each degree as needed
    by many of the differential equations models.
    
    :Arguments: 

    **G** networkx Graph
        The contact network
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.    
               
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
            
    **rho** float between 0 and 1  (default None)
            the fraction to be randomly infected at time 0
            If None, then rho=1/N is used where N = G.order()

    **SIR** (boolean)
            says whether the system will be SIR or SIS.

    :Returns: 

    Nk (scipy array)
       NUMBER (not proportion) of nodes of each degree.

    Sk0 (scipy array)
        NUMBER of susceptible nodes of each degree at t=0, 
        = (1-rho)Nk

    Ik0 (scipy array)
        NUMBER of infected nodes of each degree at t=0,   
        = rho Nk

    if SIR, also returns
    Rk0 (scipy array)
        NUMBER of recovered nodes of each degree at t=0,    
        = 0 Nk
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    if SIR is False and initial_recovereds is not None:
        raise EoN.EoNError("cannot define initial_recovereds for SIS")
        
    Nk = Counter(dict(G.degree()).values())
    maxk = max(Nk.keys())
    Nk = scipy.array([Nk[k] for k in range(maxk+1)])
    
    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, initial_recovereds)
        Sk0 = 0*Nk
        Ik0 = 0*Nk
        Rk0 = 0*Nk
        for node in G.nodes():
            k = G.degree(node)
            if status[node] == 'S':
                Sk0[k] += 1
            elif status[node] == 'I':
                Ik0[k] += 1
            else:# status[node] == 'R'
                Rk0[k] += 1  
    else:
        if rho is None:
            rho = 1./G.order()
        Sk0 = (1-rho)*Nk
        Ik0 = rho*Nk
        Rk0 = 0*Nk
    
    if SIR:
        return Nk, Sk0, Ik0, Rk0
    else:
        return Nk, Sk0, Ik0

def _get_NkNl_and_IC_as_arrays_(G, initial_infecteds=None, initial_recovereds = None, 
                                rho=None, withKs = False, SIR=True):
    r'''
    In some of the differential equations models, we need to know how
    many edges exist between nodes of given degrees in a graph.  
    
    This finds that and the initial conditions for numbers of the
    various edges assuming a fraction rho is initially infected.
    
    :Arguments: 
    **G** networkx Graph
    **initial infecteds**  node or iterable of nodes   (default None)
            if a single node, then this node is initially infected
            if an iterable, then whole set is initially infected
            if None, then choose randomly based on rho.  If rho is also
            None, a random single node is chosen.
            If both initial_infecteds and rho are assigned, then there
            is an error.       
    **initial_recovereds**  iterable of nodes (default None)
            this whole collection is made recovered.
            Currently there is no test for consistency with initial_infecteds.
            Understood that everyone who isn't infected or recovered initially
            is initially susceptible.
    **rho**  float between 0 and 1   (default None)
            the fraction to be randomly infected at time 0
            If None, then rho=1/N is used where N = G.order()
        withKs (boolean)
            flag to say whether we are restricting our attention to 
            just those degrees observed in the network or to all 
            degrees.
            If True, 
             then we only consider those degrees that are observed.
            If False, 
             then we treat it as if all degrees from 0 to kmax are 
             observed.
        SIR (boolean)
            says whether the system will be SIR or SIS.

    :Returns: 
        
        NkNl (2D scipy array)
            NUMBER (not proportion) of edges between each pair of degrees.
        SkSl0 (2D scipy array)
            initial NUMBER of edges between pair of susceptibel nodes of 
            each degree type.
            = (1-rho)^2 NkNl
        SkIl0 (2D scipy array)
            initial NUMBER of edges from a susceptible to an infected node 
            of the given degrees.
            = rho(1-rho) NkNl

        if not SIR, also returns

        IkIl0 (2D scipy array)
            initial NUMBER of edges between 2 infected nodes.  This is not 
            needed for SIR model.
            = rho^2*NkNl

        if withKs, also returns
            Ks (scipy array)
                The observed degrees in the population.
        ''' 
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")

    if withKs:
        Ks = sorted(list(set(dict(G.degree()).values())))
        klength = len(Ks)
    else:
        klength = max(dict(G.degree()).values())+1
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))

    for u,v in G.edges():
        k = G.degree(u)
        l = G.degree(v)
        NkNl[Ks.index(k)][Ks.index(l)] += 1
        NkNl[Ks.index(l)][Ks.index(k)] += 1

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, initial_recovereds)
        SkSl0 = 0*NkNl
        SkIl0 = 0*NkNl
        IkIl0 = 0*NkNl
        for u,v in G.edges():
            k = G.degree(u)
            l = G.degree(v)
            if status[u] == 'S':
                if status[v] == 'S':
                    SkSl0[Ks.index(k)][Ks.index(l)] += 1
                    SkSl0[Ks.index(l)][Ks.index(k)] += 1
                elif status[v] == 'I':
                    SkIl0[Ks.index(k)][Ks.index(l)] +=1
            elif status[u] == 'I':
                if status[v] == 'S':
                    SkIl0[Ks.index(l)][Ks.index(k)] += 1
                elif status[v] == 'I':
                    IkIl0[Ks.index(k)][Ks.index(l)] += 1
                    IkIl0[Ks.index(l)][Ks.index(k)] += 1
    else: 
        if rho is None:
            rho = 1./G.order()      
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
    r'''
    Used in several places so that we can input a graph and then we 
    can call the methods that depend on the degree distribution

    :Arguments: 

    **G** networkx Graph
    
    :Returns: 

    **Pk** dict
        Pk[k] is the proportion of nodes with degree k.
    '''

    Nk = Counter(dict(G.degree()).values())
    Pk = {x:Nk[x]/float(G.order()) for x in Nk.keys()}
    return Pk

def get_PGF(Pk):
    r'''
    Given a degree distribution (as a dict), returns the probability 
    generating function
    
    :Arguments:

    **Pk** dict
        Pk[k] is the proportion of nodes with degree k.

    :Returns: 
        
    **psi** function
            :math:`\psi(x) = \sum_k P_k[k] x^k`
    '''
    maxk = max(Pk.keys())
    ks = scipy.linspace(0,maxk, maxk+1)
    Pkarray = scipy.array([Pk.get(k,0) for k in ks])
    return lambda x: Pkarray.dot(x**ks)

def get_PGFPrime(Pk):
    r'''
    Given a degree distribution (as a dict) returns the function
    :math:`\psi'(x)`
    
    :Arguments: 

    **Pk** dict
        Pk[k] is the proportion of nodes with degree k.

    :Returns: 

    **psiPrime** (function)
        :math:`\psi'(x) = \sum_k k P_k[k] x^{k-1}`
    '''
    maxk = max(Pk.keys())
    ks = scipy.linspace(0,maxk, maxk+1)
    Pkarray = scipy.array([Pk.get(k,0) for k in ks])

    return lambda x: Pkarray.dot(ks*x**(ks-1))

def get_PGFDPrime(Pk):
    r'''
    Given a degree distribution (as a dict) returns the function 
    :math:`\psi''(x)`
    
    :Arguments: 

    **Pk** dict
        Pk[k] is the proportion of nodes with degree k.

    :Returns: 
        
    **psiDPrime** function
        :math:`\psi''(x) = \sum_k k(k-1)P_k[k] x^{k-2}`
    '''
    maxk = max(Pk.keys())
    ks = scipy.linspace(0,maxk, maxk+1)
    Pkarray = scipy.array([Pk.get(k,0) for k in ks])
    return lambda x: Pkarray.dot(ks*(ks-1)*x**(ks-2))


def get_Pnk(G):

    r'''    
    :Arguments: 

    **G** networkx Graph

    :Returns: 

    **Pnk** dict
        Pnk[k1][k2] is the proportion of neighbors of degree k1 nodes 
        that have degree k2.
    '''
    Pnk = {k1:defaultdict(int)  for k1 in dict(G.degree()).values()}
    Nk = Counter(dict(G.degree()).values())

    for node in G.nodes():
        k1 = G.degree(node)
        nbr_degrees = [G.degree(nbr) for nbr in G.neighbors(node)]
        for k2 in nbr_degrees:
            Pnk[k1][k2] += 1./(k1*Nk[k1])
    return Pnk
    
    
def estimate_R0(G, tau = None, gamma = None, transmissibility = None):
    r'''
    provides the estimate of the reproductive number R_0 = T <K^2-K>/<K>
    
    This handles the Markovian continuous time case with T = tau/(tau+gamma)
    and it handles other cases by the user inputting an average transmissibility

    For <K^2-K>/<K> it measures the network.
    
    
    :Arguments: 
    
    **G** networkx Graph

    **tau** positive float (default None)
        transmission rate
        
        Either both tau and gamma must be given or
        transmissibility is given.

    **gamma** positive float (default None)
        recovery rate 

        Either both tau and gamma must be given or
        transmissibility is given.

    **transmissibility** positive float (default None)
        average transmission probability 
        
        Either both tau and gamma must be given or
        transmissibility is given.
 
        
    :Returns:
    **R_0** float
        Reproductive number :math:`\mathcal{R}_0= T \langle K^2-K\rangle/\langle K\rangle`
    '''
    
    if transmissibility is None:
        if tau is None or beta is None:
            raise EoN.EoNError("not enough information give to estimate transmission probability")
        else:
            transmissibility = tau/(tau+gamma)
    Pk = get_Pk(G)
    psiDPrime = get_PGFDPrime(Pk)
    psiPrime = get_PGFPrime(Pk)
    return transmissibility * psiDPrime(1.)/psiPrime(1.)
            
##################
#                #
#    ODE CODE    #
#                #
##################


########      INDIVIDUAL BASED code -
########  given node and who its neighbors are, we track probability of
########  having given status based on probabilities of neighbors.  Assumes
########  independence.
    
def _dSIS_individual_based_(Y, t, G, nodelist, trans_rate_fxn, rec_rate_fxn):
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
        dY[index] = sum(trans_rate_fxn(node,nbr)*(1-Y[node])*Y[nbr] 
                            for nbr in G.neighbors(node)) - rec_rate_fxn(node)*Yi
    return dY

def _dSIR_individual_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, 
                            rec_rate_fxn):
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
        
        dX[index] = -Xi*sum(trans_rate_fxn(node,nbr)*Y[index_of_node[nbr]] 
                                for nbr in G.neighbors(node))
        dY[index] =  -dX[index] - rec_rate_fxn(node)*Yi
    dV = scipy.concatenate((dX,dY), axis=0)
    return scipy.array(dV)

def SIS_individual_based(G, tau, gamma, rho = None, Y0=None, nodelist = None, tmin = 0, 
                            tmax = 100, tcount = 1001, transmission_weight=None, 
                            recovery_weight=None, return_full_data = False):
    #tested in test_SIS_individual_based
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention 
        strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology

    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :Arguments: 

    **G** networkx Graph
        
    **tau** positive float
        transmission rate of disease

    **gamma** number 
        global recovery rate 
            
    **rho**  float between 0 and 1 (default None)
        initial uniformly random probability of being infected.
        Cannot define both rho and Y0.

    **Y0** scipy array (default None)
        the array of initial infection probabilities.
        If Y0 is defined, nodelist must also be defined.
        If Y0 is defined, rho cannot be defined.

    **nodelist** list (default None)
        list of nodes in G in the same order as in Y0.  If rho is
        defined and nodelist is not, then nodelist= G.nodes()
        Only affects returned values if return_full_data=True.

    **tmin** number       (default 0)
        minimum report time

    **tmax** number       (default 100)
        maximum report time

    **tcount** integer       (default 1001)
        number of reports

    **transmission_weighta** string       (default None)
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight** string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates so
        gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**       (default False)
        If True, returns times, Ss, Is
        if False, returns times, S, I

    :Returns: 
        
    if return_full_data is True:
        returns **times, Ss, Is**
        where times is a scipy array of times, Ss is a 2D scipy array
        Ss[i,j] gives probability nodelist[i] is susceptible at time
        times[j].
        Similarly for Is.
    if return_full data is False:
        returns **times, S, I**
        all are scipy arrays.  gives times, and expected number 
        susceptible and expected number infected.
             
    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN as EoN
        import scipy
        
        G = nx.configuration_model([3,10]*1000)
        N = G.order()
        rho = 1./N
        t, S, I = EoN.SIS_individual_based(G, 0.3, 1, rho=rho, tmax = 20)
    '''
    if nodelist is None:
        if Y0 is not None:
            raise EoN.EoNError("cannot define Y0 without defining nodelist")
        nodelist = G.nodes()

    if rho is None:
        if Y0 is None:
            raise EoN.EoNError("must define one of rho and Y0")
    else:
        if Y0 is not None:
            raise EoN.EoNError("cannot define both rho and Y0")
        else:
            Y0 = rho*scipy.ones(len(nodelist))
    
                            
    trans_rate_fxn, rec_rate_fxn = EoN._get_rate_functions_(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin, tmax, tcount)
    Y = integrate.odeint(_dSIS_individual_based_, Y0, times,  
                                    args =(G, nodelist, trans_rate_fxn, rec_rate_fxn))
    Is = Y.T
    Ss = scipy.ones(len(Is))[:,None]-Is 
    
    if return_full_data:
        return times, Ss, Is
    else:
        return times, sum(Ss), sum(Is)


def SIR_individual_based(G, tau, gamma, rho = None, Y0 = None, X0= None, 
                            nodelist = None, tmin = 0, 
                            tmax = 100, tcount = 1001, transmission_weight=None, 
                            recovery_weight=None, return_full_data = False):
    '''
    Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:

    :Arguments: 

    **G** networkx Graph

    **tau** positive float
        transmission rate of disease
    
    **gamma** number      (default None)
        global recovery rate  
        
    **rho**  number between 0 and 1 (default None)
        probability random node is infected.  Cannot be defined along with X0 and Y0
        At least one of rho and Y0 must be defined.
    **Y0**  scipy array (default None)
        the array of initial infection probabilities.  If not defined, set
        to be rho uniformly.
    **X0**  scipy array (default None)
        the array of initial susceptibility probabilities.  If not defined
        set to be 1-Y0.  Cannot define X0 without Y0.

    **nodelist**  list  (default None)
        list of nodes in G in the same order as in X0 and Y0.
        Only relevant to returned data if return_full_data=True
    
    
    **tmin**  number       (default 0)
        minimum report time

    **tmax**  number       (default 100)
        maximum report time

    **tcount**  integer       (default 1001)
        number of reports

    **transmission_weight**  string       (default None)
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**    (default False)
        If True, returns times, S, I, R, Ss, Is, Rs
        if False, returns times, S, I, R

    :Returns: 
    
    if return_full_data is True:
        returns **times, Ss, Is, Rs**
        where times is a scipy array of times, Ss is a 2D scipy array
        Ss[i,j] gives probability nodelist[i] is susceptible at time
        times[j].
        Similarly for Is ans Rs
    if return_full data is False:
        returns **times, S, I, R**
        all are scipy arrays.  gives times, and expected number 
        susceptible, expected number infected, and expected number
        recovered
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        import matplotlib.pyplot as plt
        
        G = nx.configuration_model([3,10]*10000)
        tau = 0.3
        gamma = 1
        N = G.order()
        rho = 1./N
    
        t, S, I, R = EoN.SIR_individual_based(G, tau, gamma, rho=rho, tmax = 20)
        plt.plot(t,I)
    '''

    if rho is None and Y0 is None:
        raise EoN.EoNError("must define at least one of rho and Y0")
    if X0 is not None and Y0 is None:
        raise EoN.EoNError("cannot define X0 without defining Y0") 

    
    if nodelist is None:
        if Y0 is not None:
            raise EoN.EoNError("cannot define Y0 without defining nodelist")
        nodelist = G.nodes()
           
    #nodelist is now guaranteed to exist.
    if rho is not None:
        if Y0 is not None:
            raise EoN.EoNError("cannot define both rho and Y0")
        else:
            Y0 = rho*scipy.ones(len(nodelist))
    
    if X0 is None:
        X0 = 1- Y0
        
    trans_rate_fxn, rec_rate_fxn = EoN._get_rate_functions_(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    index_of_node = {}
    for i, node in enumerate(nodelist):
        index_of_node[node] = i

    N = len(X0)
    times = scipy.linspace(tmin, tmax, tcount)
    V0 = scipy.concatenate((X0,Y0), axis=0)
    V = integrate.odeint(_dSIR_individual_based_, V0, times, 
                            args = (G, nodelist, index_of_node, 
                                    trans_rate_fxn, rec_rate_fxn))
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


def SIS_individual_based_pure_IC(G, tau, gamma, initial_infecteds, nodelist = None,
                                    tmin = 0, tmax = 100, tcount = 1001, 
                                    transmission_weight=None, 
                                    recovery_weight=None, 
                                    return_full_data = False):
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIS_individual_based is that this 
    one assumes a "pure initial condition", that is, we know exactly 
    what the statuses of the nodes are at the initial time.  
    
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate of disease
    **gamma** number      (default None)
        global recovery rate  

    **initial_infecteds**  list or set
        the set of nodes initially infected

    **nodelist**  list  (default None)
        list of nodes in G in desired order. (only matters if return_full_data==True)

    **tmin**  number       (default 0)
        minimum report time

    **tmax**  number       (default 100)
        maximum report time

    **tcount**  integer       (default 1001)
        number of reports

    **transmission_weight**  string       (default None)
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**  boolean      (default False)


    :Returns: 
    if return_full_data is True,
        returns **times, Ss, Is**
    if return_full_data is False,
        returns **times, S, I**

    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        import matplotlib.pyplot as plt
        
        G = nx.configuration_model([3,10]*1000)
        nodelist = G.nodes()
        initial_infecteds = range(100)
        t, S, I = EoN.SIS_individual_based(G, 0.3, 1, initial_infecteds, nodelist,
                                            tmax = 20)
        plt.plot(t,I)

    '''
    #make Y0[u] be 1 if infected 0 if not
    if nodelist is None:
        nodelist = G.nodes()
    initial_infecteds = set(initial_infecteds)
    Y0 = scipy.array([1 if u in initial_infecteds else 0 for u in nodelist])

    return SIS_individual_based(G, tau, gamma, Y0=Y0, nodelist=nodelist, 
                                tmin=tmin, tmax=tmax, tcount=tcount,
                                transmission_weight=transmission_weight, 
                                recovery_weight = recovery_weight, 
                                return_full_data = return_full_data)
        


def SIR_individual_based_pure_IC(G, tau, gamma, initial_infecteds, 
                                    initial_recovereds = None, nodelist=None, 
                                    tmin = 0, tmax = 100, tcount = 1001, 
                                    transmission_weight=None, 
                                    recovery_weight=None, 
                                    return_full_data = False):
    '''Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIR_individual_based is that this 
    one assumes a "pure initial condition", that is, we know exactly 
    what the statuses of the nodes are at the initial time.
    
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate of disease

    **gamma** number      (default None)
        global recovery rate  
 
    **initial_infecteds**  list or set
        the set of nodes initially infected
      
    **nodelist**  list  (default None)
        list of nodes in G in desired order. (only matters if return_full_data==True)
   
    **initial_recovereds** list or set  (default None)
        initially recovered nodes
        if equal to None, then all non-index nodes are initially 
        susceptible.

    **tmin**  number       (default 0)
        minimum report time

    **tmax**  number       (default 100)
        maximum report time

    **tcount**  integer       (default 1001)
        number of reports

    **transmission_weight**  string       (default None)
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**  boolean      (default False)

    :Returns: 
    if return_full_data is True,
        returns **times, S, I, R, Ss, Is, Rs**
    if return_full_data is False,
        returns **times, S, I, R**
    
    '''
    if nodelist is None:
        nodelist = G.nodes()
    N = len(nodelist)
    initial_infecteds = set(initial_infecteds)
    #make Y0[u] be 1 if infected 0 if not
    Y0 = scipy.array([1 if u in initial_infecteds else 0 for u in nodelist])
    if initial_recovereds is None:
        X0 = 1 - Y0
    else:
        initial_recovereds = set(initial_recovereds) #for fast membership test
        non_susceptibles = initial_recovereds.union(initial_infecteds)
        X0 = scipy.array([0 if u in non_susceptibles  else 1
                            for u in nodelist])
    
    return SIR_individual_based(G, tau, gamma, nodelist, X0, Y0, tmin, 
                                tmax, tcount, transmission_weight, 
                                recovery_weight, return_full_data)

########   PAIR BASED

def _dSIS_pair_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, rec_rate_fxn):
    '''
    <\dot{Y}_i> = tau \sum_j g_{ij} <XiYj>  -  gamma_i <Yi>
    <\dot{XY}_ij> = tau sum_{k \neq i} g_{jk} <XiXj><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XiYj>/<Xi>
                   - tau g_{ij}<XiYj> - gamma_j <XiYj>  
                   + ***gamma_i <YiYj>***
    <\dot{XX}_ij> = - tau sum_{k\neq i} g_{jk} <XiXj><XjYk>/<Xj>
                    - tau sum_{k neq j} g_{ik} <YkXi><XiXj>/<Xi>
                    + **** \gamma_i <YiXj> + gamma_j <XiYj>****

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
    All derivatives are initialized to 0, and then the loop only makes 
    changes for those terms where an edge exists.
    '''
    #print(t)
    N=G.order()
    Y = V[0:N] #infecteds
    X = 1-Y    #susceptibles
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
            #there are places where we divide by X[i] which may = 0.
            #In those cases the numerator is (very) 0, so it's easier
            #to set this up as mult by inverse with a dummy value when
            #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[N: N+N**2]
    XX = V[N+N**2:]

    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T   #not really needed, 
                #but helps keep consistent with equations as written.

    YY = 1 - XY-XX-YX

    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, but I think this sequence of for 
    # loops is easier to read, or at least understand.
    #I expect this isn't the bottleneck.  
    #Will avoid (premature) optimization for now.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate_fxn(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dY[i] += trans_rate_fxn(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*XY[i,j] \
                            + rec_rate_fxn(u)*YY[i,j]
            dXX[i,j] +=  rec_rate_fxn(u)*YX[i,j] + rec_rate_fxn(v)*XY[i,j]
            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(v):
                if w == u: #skip these
                    continue
                #so w != v. 
                k= index_of_node[w]

                dXY[i,j] += trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j] 
                dXX[i,j] += -trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j] 
            for w in G.neighbors(u):
                if w == v:
                    continue #skip these
                k = index_of_node[w]
                dXY[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XX[i,j]*Xinv[i]


    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dY[:,None], dXY, dXX), axis=0).T[0]
    return dV

def _dSIR_pair_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, rec_rate_fxn):
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
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
            #there are places where we divide by X[i] which may = 0.
            #In those cases the numerator is (very) 0, so it's easier
            #to set this up as mult by inverse with a dummy value when
            #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[2*N: 2*N+N**2]
    XX = V[2*N+N**2:]

    #print X.shape, Y.shape, XY.shape, XX.shape, N
    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T   #not really needed, 
                #but helps keep consistent with equations as written.

    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, 
    #but I think this sequence of for loops is easier to read, 
    #or at least understand.
    #I expect it to run quickly regardless.  Will avoid (premature) 
    #optimization for now.
    
    #Okay -this does not appear to run quickly (at least for a complete graph). 
    #I'll need to do some profiling.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate_fxn(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dX[i] += -trans_rate_fxn(u,v)*XY[i,j]
            dY[i] += trans_rate_fxn(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*XY[i,j] 

            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(v):
                if w == u: #skip these
                    continue
                #so w != u.  
                k= index_of_node[w]
                #i corresponds to u, j to v and k to w.
                dXY[i,j] += trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  
                dXX[i,j] += -trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j] 
            for w in G.neighbors(u):
                if w == v:
                    continue #skip these
                k = index_of_node[w]
                dXY[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XX[i,j]*Xinv[i]

    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dX[:, None], dY[:,None], dXY, dXX), axis=0).T[0]
    return dV


def SIS_pair_based(G, tau, gamma, rho = None, nodelist = None,
                    Y0=None, XY0=None, XX0 = None, tmin = 0, tmax = 100,  
                    tcount = 1001, transmission_weight=None, 
                    recovery_weight=None, return_full_data = False):
    r'''
    Encodes System (3.26) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    WARNING: this does NOT solve the pairwise equations.  Look at
    SIS_homogeneous_pairwise and SIS_heterogeneous_pairwise for that.
    
    This system solves equations for an SIS disease model spreading on a 
    given graph.  
    
    It captures the dependence with pairs, but not triples.

    It does not include corrections for triangles (or any other cycles).  
    
    The corrections for triangles are provided in the text, but not 
    implemented here.

    There are some inefficiencies in the implementation:
        we track all pairs, rather than just those pairs in edges, but 
        this is unlikely to significantly affect the calculation time.  
        
        This makes it much easier to vectorize things.
        
        We track pairs in both directions: e.g., XX[1,2] and XX[2,1].


    :Arguments: 
    
    **G** networkx Graph
        
    **tau** positive float
        transmission rate of disease

    **gamma** number
        global recovery rate  
            
    **rho**  (float between 0 and 1, default None)
        proportion assumed initially infected.  If None, then Y0 is used
        if Y0 is also None, then rho = 1./N

    **nodelist**  list
        list of nodes in G in some prescribed order (just since there is no 
        guarantee that G returns nodes in the same order if things change 
        a bit.)
            
    **Y0**  scipy array
        the array of initial infection probabilities for each node in 
        the same order as in nodelist

    **XY0** 2D scipy array (default None)
        (each dimension has length number of nodes of G)
        XY0[i,j] is probability node i is susceptible and j is 
        infected.
        if None, then assumes that infections are introduced 
        randomly according to Y0.

    **XX0** 2D scipy array (default None)
        (each dimension has length number of nodes of G)
        XX0[i,j] is probability nodes i and j are susceptible.
        if None, then assumes that infections are introduced 
        randomly according to Y0.

    **tmin**  number (default 0)
        minimum report time

    **tmax**  number (default 100)
        maximum report time 

    **tcount**  integer (default 1001)
        number of reports

    **transmission_weight**  string
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**  boolean      (default False)
        if True:
            returns times, S, I, R, Xs, Ys, Zs, XY, XX
        if False:
            returns times, S, I, R

    :Returns: 
    if return_full_data is True:
        returns **times, S, I, Xs, Ys, XY, XX**
    if False:
        returns times, **S, I**
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
        G = nx.fast_gnp_random_graph(1000,0.004)
        nodelist = G.nodes()
        Y0 = scipy.array([1 if node<10 else 0 for node in nodelist]) #infect first 10
        t, S, I = EoN.SIS_pair_based(G, 2, 0.5, nodelist, Y0, tmax = 4, tcount = 101)
        plt.plot(t,I)
        
'''

    N = G.order()
        
    if Y0 is None and rho is None:
        rho = 1./N
    if Y0 is not None and rho is not None:
        raise EoN.EoNError("either Y0 or rho must be defined")
    if Y0 is not None and  nodelist is None:
        raise EoN.EoNError("cannot define Y0 without nodelist")
        
    if  nodelist is None: #only get here if Y0 is None
        nodelist = G.nodes()
        Y0 = scipy.array([rho]*N)
    if len(Y0) != N:
        raise EoN.EoNError("incompatible length for Y0")            

    trans_rate_fxn, rec_rate_fxn = EoN._get_rate_functions_(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin,tmax,tcount)


    
    X0=1-Y0

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    elif XY0.shape != (N,N):
        raise EoN.EoNError("incompatible lengths for XY0 and Y0")

    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    elif XX0.shape != (N,N):
        raise EoN.EoNError("incompatible lengths for XX0 and Y0")
        
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs not
    XX0 = XX0*A  #in an edge, but this led to an error with odeint.  
                 #Multiplying by A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((Y0[:,None], XY0, XX0), axis=0).T[0]
    index_of_node = {node:i for i, node in enumerate(nodelist)}

    V = _my_odeint_(_dSIS_pair_based_, V0, times, 
                            args = (G, nodelist, index_of_node, trans_rate_fxn, 
                                    rec_rate_fxn))
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


def SIS_pair_based_pure_IC(G, tau, gamma, initial_infecteds, nodelist = None,
                            tmin = 0, tmax = 100, tcount = 1001,
                            transmission_weight = None, recovery_weight=None,
                            return_full_data = False):
    r''' 
    Encodes System (3.26) of Kiss, Miller, & Simon, using a "pure initial 
    condition".  That is, we can specify the exact status of all nodes at tmin
    
    Please cite the book if using this algorithm.
    
    :Arguments: 
    
    **G** networkx Graph
        
    **tau** positive float
        transmission rate of disease

    **gamma** number
        global recovery rate  
            
    **initial_infecteds**  list or set
        the set of nodes initially infected
        
    **nodelist**  list
        list of nodes in G in some prescribed order (just since there is no 
        guarantee that G returns nodes in the same order if things change 
        a bit.)
            
    **tmin**  number (default 0)
        minimum report time

    **tmax**  number (default 100)
        maximum report time 

    **tcount**  integer (default 1001)
        number of reports

    **transmission_weight**  string
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**  boolean      (default False)
        if True:
            returns times, S, I, R, Xs, Ys, Zs, XY, XX
        if False:
            returns times, S, I, R

    :Returns: 
    if return_full_data is True:
        returns **times, S, I, Xs, Ys, XY, XX**
    if False:
        returns times, **S, I**
    
    '''
    if nodelist is None:
        nodelist = G.nodes()
    N = len(nodelist)
    #make Y0[u] be 1 if infected 0 if not
    initial_infecteds = set(initial_infecteds) #for fast membership test
    Y0 = scipy.array([1 if u in initial_infecteds else 0 for u in nodelist])    
    
    return SIS_pair_based(G, tau, gamma, nodelist = nodelist,
                            Y0=Y0, tmin=tmin, tmax=tmax, tcount=tcount,
                            transmission_weight=transmission_weight, 
                            recovery_weight=recovery_weight,
                            return_full_data=return_full_data)

def _SIR_pair_based_initialize_node_data(G, rho, nodelist, X0, Y0):
    #inputs must define either rho or Y0.  In first case nodelist is optional.
    if (rho and Y0 is not None) or (not rho and Y0 is None):
        raise EoN.EoNError("need rho or Y0 defined for initial condition, \
                        but not both")
    if Y0 is not None and nodelist is None:
        raise EoN.EoNError("order in Y0 is ambiguous if nodelist is not given.")

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

def _SIR_pair_based_initialize_edge_data(G, edgelist, nodelist, XY0, YX0, 
                                            XX0, X0, Y0, index_of_node):
    if (not XY0 is None or YX0 is None or XX0 is None) \
            and (XY0 is not None or YX0 !=None  or XX0 is not None):  
            #at least one defined and one not defined
        raise EoN.EoNError("must define all of XY0, YX0, and XX0 or none of them")
    if not edgelist:
        if XY0:
            raise EoN.EoNError("order in XY0, YX0, and XX0 is ambiguous if \
                            edgelist is not given.")
        else:
            edgelist = list(G.edges())

    if XY0:
        #test that XY0 <= X0Y0, same for  YX0 and XX0
        for index,(u,v) in enumerate(edgelist):
           i_u = index_of_node[u]
           i_v = index_of_node[v]
           if XY0[index] >X0[i_u]*Y0[i_v] or YX0[index]>Y0[i_u]*X0[I_v] \
                                or XX0[index]>X0[i_u]*X0[i_v]:
               raise EoN.EoNError("edge probabilities inconsistent with node \
                                probabilities")
    else:
        XY0 = scipy.array([X0[index_of_node[u]]*Y0[index_of_node[v]] 
                            for u,v in edgelist])
        YX0 = scipy.array([Y0[index_of_node[u]]*X0[index_of_node[v]] 
                            for u,v in edgelist])
        XX0 = scipy.array([X0[index_of_node[u]]*X0[index_of_node[v]] 
                            for u,v in edgelist])
    return edgelist, XY0, YX0, XX0


def SIR_pair_based(G, tau, gamma, rho = None, nodelist=None, Y0=None, 
                    X0 = None, XY0=None, 
                    XX0 = None, tmin = 0, tmax = 100, tcount = 1001, 
                    transmission_weight=None, recovery_weight=None, 
                    return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    WARNING This does not solve the pairwise equations.  Look at 
    SIR_homogeneous_pairwise and SIR_heterogeneous_pairwise for that.
    
    This system solves equations for an SIR disease model spreading on a
    given graph.  It captures the dependence with pairs, but not 
    triples.
    
    It will be exact for a tree.

    There are NO CORRECTIONS for the existence of TRIANGLES or any other
    CYCLES.
    
    Some corrections for triangles are provided in the text, but not 
    implemented here.
    
    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention 
        strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology
    

    :Arguments: 
    **G** networkx Graph

    **tau** positive float
            transmission rate of disease

    **gamma** number
            global recovery rate  
    
    **rho**  float between 0 and 1   (default None)
            proportion assumed initially infected.  If None, then Y0 is used
            if Y0 is also None, then rho = 1./N

    **nodelist**  list
            list of nodes in G in the some prescribed order (just since there is no 
            guarantee that G returns nodes in the same order if things change 
            a bit.)

    **Y0**  scipy array
            the array of initial infection probabilities for each node in 
            order as in nodelist/

    **X0**  scipy array (default None)
            probability a random node is initially susceptible.
            the probability of initially recovered will be 1-X0-Y0.  By 
            default we assume no initial recoveries, so X0=1-Y0 will be 
            assumed unless both Y0 and X0 are given.

    **XY0** 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is 
                infected.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    **XX0** 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    **tmin**  number (default 0)
            minimum report time

    **tmax**  number (default 100)
            maximum report time 

    **tcount**  integer (default 1001)
            number of reports

    **transmission_weight**  string
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    **return_full_data**  boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    :Returns: 

    if return_full_data is True:
        returns times, S, I, R, Xs, Ys, Zs, XY, XX
    if ... is False:
        returns times, S, I, R


    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
        G = nx.fast_gnp_random_graph(1000,0.004)
        nodelist = G.nodes()
        Y0 = scipy.array([1 if node<10 else 0 for node in nodelist]) #infect first 10
        t, S, I, R = EoN.SIR_pair_based(G, nodelist, Y0, 2, 0.5, tmax = 4, tcount = 101)
        plt.plot(t,I)
    '''

    N = G.order()
        
    if Y0 is None and rho is None:
        rho = 1./N
    if Y0 is not None and rho is not None:
        raise EoN.EoNError("either Y0 or rho must be defined")
    if Y0 is not None and  nodelist is None:
        raise EoN.EoNError("cannot define Y0 without nodelist")
        
    if  nodelist is None: #only get here if Y0 is None
        nodelist = G.nodes()
        Y0 = scipy.array([rho]*N)
    if len(Y0) != N:
        raise EoN.EoNError("incompatible length for Y0")            

    trans_rate_fxn, rec_rate_fxn = EoN._get_rate_functions_(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin,tmax,tcount)
    if X0 is None:
        X0 = 1-Y0

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    else:
        if XY0.shape != (N,N):
            raise EoN.EoNError("incompatible lengths for XY0 and Y0")
    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    else:
        if XX0.shape != (N,N):
            raise EoN.EoNError("incompatible lengths for XX0 and Y0")
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs not
    XX0 = XX0*A  #in an edge, but this led to the error with odeint.  
                 #Multiplying by A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((X0[:,None], Y0[:,None], XY0, XX0), axis=0).T[0]
    #print V0.shape
    index_of_node = {node:i for i, node in enumerate(nodelist)}
    #    index_of_node = {}
    #    for i, node in enumerate(nodelist):
    #        index_of_node[node] = i


    V = integrate.odeint(_dSIR_pair_based_, V0, times, 
                            args = (G, nodelist, index_of_node, trans_rate_fxn, 
                                    rec_rate_fxn))
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


def SIR_pair_based_pure_IC(G, tau, gamma, initial_infecteds, 
                            initial_recovereds = None, nodelist=None, 
                            tmin = 0, tmax = 100, tcount = 1001, 
                            transmission_weight=None, recovery_weight=None, 
                            return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon, using a "pure initial 
    condition".  Please cite the
    book if using this algorithm.

    uses SIR_pair_based after finding the appropriate initial conditions.
    
    :Arguments: 
    
    **G** networkx Graph
        
    **tau** positive float
        transmission rate of disease

    **gamma** number
        global recovery rate  
            
    **initial_infecteds**  list or set
        the set of nodes initially infected
        
    **initial_recovereds** list or set (default None)
        the set of nodes that are already recovered.
        
    **nodelist**  list
        list of nodes in G in a prescribed order (just since there is no 
        guarantee that G returns nodes in the same order if things change 
        a bit.)
            
    **tmin**  number (default 0)
        minimum report time

    **tmax**  number (default 100)
        maximum report time 

    **tcount**  integer (default 1001)
        number of reports

    **transmission_weight**  string
        the label for a weight given to the edges.
        G.edge[i][j][transmission_weight] = g_{ij}

    **recovery_weight**  string       (default None)
        a label for a weight given to the nodes to scale their 
        recovery rates
            gamma_i = G.node[i][recovery_weight]*gamma
    **return_full_data**  boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    :Returns: 

    if return_full_data is True:
        returns **times, S, I, R, Xs, Ys, Zs, XY, XX**
    if ... is False:
        returns times, S, I, R
    '''
    
    if nodelist is None:
        nodelist = G.nodes()
    N = len(nodelist)
    #make Y0[u] be 1 if infected 0 if not
    initial_infecteds = set(initial_infecteds) #for fast membership test
    Y0 = scipy.array([1 if u in initial_infecteds else 0 for u in nodelist])
    if initial_recovereds is None:
        X0 = 1 - Y0
    else:
        initial_recovereds = set(initial_recovereds) #for fast membership test
        non_susceptibles = initial_recovereds.union(initial_infecteds)
        X0 = scipy.array([0 if u in non_susceptibles  else 1
                            for u in nodelist])
    
    return SIR_pair_based(G, tau, gamma, nodelist=nodelist, Y0=Y0, X0=X0,
                            tmin=tmin, tmax=tmax, tcount=tcount,
                            transmission_weight=transmission_weight, 
                            recovery_weight=recovery_weight,
                            return_full_data=return_full_data)
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

           
def SIS_homogeneous_meanfield(S0, I0, n, tau, gamma, tmin=0, tmax=100, 
                                tcount=1001):
    '''Encodes System (4.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    In the text this is often referred to as the 
    "mean-field model closed at the level of pairs"

       [\dot{S}] = \gamma [I] - tau n[S][I]/N
       
       [\dot{I}] = \tau n[S][I]/N - \gamma [I]

    This is the SIS version of the "Kermack-McKendrick equations".
    
    :Arguments: 

    **S0** number
        initial number susceptible
    **I0** number
        initial number infected
    **n** integer
        (average) degree of all nodes.
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports

    :Returns: 
    **times, S, I**    all scipy arrays
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 999
        I0 = 1
        n = 4 #degree
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_homogeneous_meanfield(S0, I0, n, tau, gamma)
    '''

    N=S0+I0
    X0=scipy.array([S0,I0])
    times=scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_meanfield_, X0,times, 
                            args=(float(n)/N, tau, gamma))
    S, I= X.T
    return times, S, I

def SIR_homogeneous_meanfield(S0, I0, R0, n, tau, gamma, tmin=0, tmax=100, 
                                tcount=1001):
    '''Encodes System (4.9) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of pairs"
    
    These are often referred to as the "Kermack-McKendrick equations"

    [\dot{S}] = - tau n[S][I]/N
    [\dot{I}] = \tau n[S][I]/N - \gamma [I]
    [\dot{R}] = \gamma [I]


    :Arguments: 

    **S0** number
        initial number susceptible
    **I0** number
        initial number infected
    **R0** number
        initial number recovered
    **n** integer
        degree of each node
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    
    :Returns: 
    **times, S, I, R** all scipy arrays
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 999
        I0 = 1
        n = 4 #degree
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_homogeneous_meanfield(S0, I0, 0, n, tau, gamma)
        
    '''

    N=S0+I0+R0
    X0= scipy.array([S0,I0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_meanfield_, X0, times, 
                            args=(float(n)/N, tau, gamma))
    S, I= X.T
    R = N-S-I
    return times, S, I, R

def SIS_homogeneous_meanfield_from_graph(G, tau, gamma, 
                                        initial_infecteds=None, rho = None, 
                                        tmin = 0, tmax=100, tcount=1001):
    r'''
    Calls SIR_homogeneous_pairwise after calculating S0, I0, & n, 
    based on the graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports

    :Returns: 
    **times, S, I**    all scipy arrays
'''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    kave = G.size()*2.0/G.order()
    if initial_infecteds is not None:
        I0 = len(initial_infecteds)
    elif rho is not None:
        I0 = rho*G.order()
    else:
        I0 = 1.
    S0 = G.order()-I0
    return SIS_homogeneous_meanfield(S0, I0, kave, tau, gamma, tmin=tmin, tmax=tmax, 
                                tcount=tcount)
                                
def SIR_homogeneous_meanfield_from_graph(G, tau, gamma, initial_infecteds=None, 
                                    initial_recovereds = None, rho = None, 
                                    tmin = 0, tmax=100, tcount=1001):
    r'''
    Calls SIR_homogeneous_pairwise after calculating S0, I0, R0, & n, 
    based on the graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports

    :Returns: 
    **times, S, I, R** all scipy arrays
'''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    kave = G.size()*2.0/G.order()
    if initial_infecteds is not None:
        I0 = len(initial_infecteds)
        if initial_recovereds is None:
            initial_recovereds = []
    elif rho is not None:
        I0 = rho*G.order()
    else:
        I0 = 1.
    R0 = len(initial_recovereds)
        
    S0 = G.order()-I0 - R0
    return SIR_homogeneous_meanfield(S0, I0, R0, kave, tau, gamma, tmin=tmin, tmax=tmax, 
                                tcount=tcount)

####   HOMOGENEOUS PAIRWISE

def _dSIS_homogeneous_pairwise_(X, t, N, n, tau, gamma):
    r'''
    [\dot{S}] = gamma [I] - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I] = -[\dot{S}]
    [\dot{SI}] = \gamma([II]-[SI])+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] 
                 - \tau [SI]
    [\dot{SS}] = 2\gamma[SI] - 2\tau ((n-1)/n) [SI][SS]/[S]
    [\dot{II}] = -2\gamma[II] + 2\tau((n-1)/n) [SI]^2/[S] + 2\tau[SI]

    conserved quantities: [S]+[I]
                          [SS]+2[SI]+[II]

    n([S]+[I]) should equal [SS]+2[SI]+[II], so II will be calculated 
    based on this.
    '''
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

def SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin = 0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    r'''Encodes System (4.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of triples"

    :Arguments: 

    **S0** number
        initial number susceptible
    **I0** number
        initial number infected
    **SI0** number
        initial number of SI edges
    **SS0** number
        initial number of SS edges
    **n** integer
        (common) degree of nodes.
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I or all calculated data.
    
    :Returns: 

    if return_full_data is True:
        **t, S, I, SI, SS, II**
    if return_full_data is False:
        **t, S, I**
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 990
        I0 = 10
        SI0 = 50
        SS0 = 4900
        n = 5
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, 
                                                tmax = 20)
    '''
    N = S0+I0

    if SS0 + SI0*2>n*N:
        raise EoN.EoNError('Initial condition has more SS, SI, and IS edges than allowed')

    X0 = scipy.array([S0, SI0, SS0])
    
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_pairwise_, X0, times, 
                            args=(N, n, tau, gamma))
    S, SI, SS= X.T
    I = N-S
    
    if return_full_data:
        II = N*n - SS-2*SI
        return times, S, I, SI, SS, II
    else:
        return times, S, I
    
    
def SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin = 0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    '''Encodes System (4.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of triples"

    [\dot{S}] = - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I]
    [\dot{R}] = \gamma [I]    ;    [R] = N-[S]-[I]
    [\dot{SI}] = -\gamma [SI]+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] 
                 - \tau [SI]
    [\dot{SS}] = - 2\tau ((n-1)/n) [SI][SS]/[S]

    conserved quantities: [S]+[I]+[R]  also 
                          [SS]+[II]+[RR] + 2([SI] + [SR] + [IR])

    :Arguments: 

    **S0** float
        Initial number susceptible
    **I0** float
        Initial number infected
    **R0** float
        Initial number recovered
    **SI0** float
        Initial number of SI edges
    **SS0** float
        Initial number of SS edges
    **n** float
        Degree of nodes
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        tells whether to just return times, S, I, R or all calculated data.
        if True, then returns times, S, I, R, SI, SS
    :Returns: 
    if return_full_data is True:
        **times, S, I, R, SI, SS**
    if return_full_data is False:
        **times, S, I, R **

    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        S0 = 990
        I0 = 10
        R0 = 1
        SI0 = 45
        SS0 = 4900
        n = 5
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, 
                                                    tmax = 20)

    '''
    N = S0+I0+R0
    if SS0 + 2*SI0 > n*N:
        raise EoN.EoNError('Initial condition has more SS, SI, and IS edges than allowed')
    X0 = scipy.array([S0, I0, SI0, SS0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_pairwise_, X0, times, 
                            args=(n, tau, gamma))
    S, I, SI, SS = X.T
    R = N-S-I
    if return_full_data:
        return times, S, I, R, SI, SS
    else:
        return times, S, I, R


def SIS_homogeneous_pairwise_from_graph(G, tau, gamma, initial_infecteds=None, 
                                        rho = None, tmin = 0, tmax=100, 
                                        tcount=1001, return_full_data=False):
    r'''
    Calls SIS_homogeneous_pairwise after calculating S0, I0, SI0, SS0, n based
    on the graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
        the contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        tells whether to just return times, S, I, or all calculated data.
        if True, then returns times, S, I, SI, SS

    :Returns: 
        
    if return_full_data is True:
        **t, S, I, SI, SS, II**
    if return_full_data is False:
        **t, S, I**
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.fast_gnp_random_graph(10000,0.0005)
        tau = 1
        gamma = 3
        rho = 0.02
        t, S, I = EoN.SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho, tmax = 20)
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")

    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds)
        I0= len(initial_infecteds)
        S0 = N-I0
        SS0=0
        II0=0
        SI0=0
        for edge in G.edges():
            if status[edge[0]]==status[edge[1]]:
                if status[edge[0]]=='S':
                    SS0 += 2
                else:
                    II0 += 2
            else:
                SI0 += 1
    else:
        if rho is None:
            rho = 1./G.order()
        S0 = (1-rho)*N
        I0 = rho*N
        SI0 = (1-rho)*N*n*rho
        SS0 = (1-rho)*N*n*(1-rho)
    return SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin, 
                                        tmax, tcount, return_full_data)

def SIR_homogeneous_pairwise_from_graph(G, tau, gamma, initial_infecteds=None, 
                                        initial_recovereds = None,
                                        rho = None, tmin = 0, 
                                        tmax=100, tcount=1001, 
                                        return_full_data=False):
    r'''
    Calls SIR_homogeneous_pairwise after calculating S0, I0, R0, SI0, SS0, n, 
    based on the graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
        the contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        tells whether to just return times, S, I, R or all calculated data.
        if True, then returns times, S, I, R, SI, SS

    :Returns: 
        
    if return_full_data is True:
        **t, S, I, SI, SS, II**
    if return_full_data is False:
        **t, S, I**
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.fast_gnp_random_graph(10000,0.0005)
        tau = 1
        gamma = 3
        rho = 0.02
        t, S, I, R = EoN.SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho, 
                                                                tmax = 20)
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    
    
    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()

    if initial_infecteds is not None:
        if initial_recovereds is None:
            initial_recovereds = []
        status = _initialize_node_status_(G, initial_infecteds, initial_recovereds)
        I0 = len(initial_infecteds)
        R0 = len(initial_recovereds)
        S0 = N-I0-R0
        SS0 = 0
        SI0 = 0
        for edge in G.edges():
            if status[edge[0]] == 'S' and status[edge[1]] == 'S':
                SS0 += 2
            else:
                if ((status[edge[0]] == 'S' and status[edge[1]] == 'I') or
                    (status[edge[0]] == 'I' and status[edge[1]] == 'S')):
                    SI0 += 1
    else:
        if rho is None:
            rho = 1./G.order()
        S0 = (1-rho)*N
        I0 = rho*N
        SI0 = (1-rho)*N*n*rho
        SS0 = (1-rho)*N*n*(1-rho)
        R0=0
    return SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin,
                                    tmax, tcount, return_full_data)




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
    Rk = scipy.array(X[1:])
    ks = scipy.arange(len(Rk))
    Sk = S0*(theta**ks)
    Ik = Nk - Sk - Rk
    pi_I = ks.dot(Ik)/ks.dot(Nk)
    dRkdt = gamma*Ik
    dThetadt = - tau *pi_I * theta

    dX = scipy.concatenate(([dThetadt],dRkdt), axis=0)
    return dX




def SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin = 0, tmax=100, 
                                tcount=1001, return_full_data=False):
    '''Encodes System (5.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of pairs"

    This is also called Degree-baded Mean Field or Mean Field Social 
    Heterogeneity
    
    a few notes on the inputs:
    Sk0 is an array (or a list). 
    
    It is not a dict.  
    
    Sk0[k] is the *number* of nodes that are susceptible and have degree 
    k (even if some degrees missing).  
    
    A dict like this can be converted into an array by
    Sk0 = scipy.array([Sk0dict.get(k,0) 
                       for k in range(max(Sk0dict.keys())+1)])

    Ik0 is similar to Sk0.


    [\dot{S}_k] = \gamma [I_k] - \tau k [S_k] \pi_I
    [\dot{I}_k] = -(above)
    \pi_I = \sum_k k [I_k] / \sum_k  k [N_k]



    :Arguments: 

    **Sk0** scipy array
        number susceptible for each k
    **Ik0** scipy array
        number infected for each k
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        tells whether to just return times, S, I or all calculated data.
        if True, returns t, S, I, Sk, Ik

    :Returns: 
        
    if return_full_data is True:
        **times, S, I, Sk**  (Sk is scipy 2D arrays)
    if return_full_data is False:
        **times, S, I**      (all scipy arrays)
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        Sk0 = [995, 995, 995, 995, 995]
        Ik0 = [5, 5, 5, 5, 5]
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmax = 10)
    '''
    if len(Sk0) != len(Ik0):
        raise EoN.EoNError('length of Sk0 not equal to length of Ik0')
    Sk0 = scipy.array(Sk0)
    Ik0 = scipy.array(Ik0)

    kcount = len(Sk0)
    
    X0 = scipy.concatenate((Sk0,Ik0), axis=0)
    Nk = Sk0+Ik0
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIS_heterogeneous_meanfield_, X0, times, 
                            args=(kcount,tau,gamma))
    Sk = scipy.array(X.T[:kcount])
    Ik = scipy.array(X.T[kcount:])
    S = Sk.sum(axis=0)
    I = Ik.sum(axis=0)
    if return_full_data:
        return times, S, I, Sk, Ik
    else:
        return times, S, I
	
    
def SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin = 0, tmax=100, 
                                    tcount=1001, return_full_data=False):
    '''
    Encodes System (5.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of pairs"

    This is also called Degree-baded Mean Field or Mean Field Social 
    Heterogeneity
    
    Ik0 and Rk0 are similar to Sk0.

    [S_k] = [S_k](0) theta^k
    [I_k] = [N_k] - [S_k] - [R_k]
    [\dot{R}_k] = \gamma [I_k]
    pi_I = \sum_k k[I_k]


    :Arguments: 

    **Sk0** array
        Sk0[k] is the number of
        nodes that are susceptible and have degree k (even if some degrees 
        missing).
    **Ik0** array
        as in Sk0
    **Rk0** array
        as in Sk0
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
        
    if return_full_data is True:
        **times, S, I, R, Sk, Ik, Rk** (the Xk are scipy 2D arrays)
    if return_full_data is False:
        **times, S, I, R**          (all scipy arrays)
        
    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        Sk0 = [995, 995, 995, 995, 995]
        Ik0 = [5, 5, 5, 5, 5]
        Rk0 = [0,0,0,0,0]
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, 
                                                        tmax = 10)
    '''
    if len(Sk0) != len(Ik0) or len(Sk0) != len(Rk0):
        raise EoN.EoNError('length of Sk0, Ik0, and Rk0 must be the same')

    theta0=1
    Sk0 = scipy.array(Sk0)
    Ik0 = scipy.array(Ik0)
    Rk0 = scipy.array(Rk0)
    Nk = Sk0+Ik0 +Rk0
    X0 = scipy.concatenate(([theta0],Rk0), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_heterogeneous_meanfield_, X0, times, 
                            args = (Sk0, Nk, tau, gamma))
    
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
    
def SIS_heterogeneous_meanfield_from_graph(G, tau, gamma,  initial_infecteds=None, 
                                            rho = None, tmin = 0, tmax=100, 
                                            tcount=1001, return_full_data=False):
    r'''
    Calls SIS_heterogeneous_meanfield after calculating Sk0, Ik0 based on 
    the graph G and random fraction infected rho.

    :Arguments: 
    **G** networkx Graph
        the contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
        
    if return_full_data is True:
        **times, S, I, Sk, Ik**  (the Xk are scipy 2D arrays)
    if return_full_data is False:
        **times, S, I**         (all scipy arrays)
        
    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        G = nx.configuration_model([1,2,3,4]*1000)
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, 
                                                                tmax = 15)

    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")

    Nk, Sk0, Ik0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds=initial_infecteds,
                                            rho=rho, SIR=False)
    
    return SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin, tmax, tcount, return_full_data)

def SIR_heterogeneous_meanfield_from_graph(G, tau, gamma,  initial_infecteds=None, 
                                            initial_recovereds = None, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''
    Calls SIR_heterogeneous_meanfield after calculating Sk0, Ik0, Rk0 based
    on a graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
        the contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all calculated data.

    :Returns: 

    if return_full_data is True
        **times, Sk, Ik, Rk**     (the Xk are scipy 2D arrays)
    if False,
        **times, S, I, R**     (all scipy arrays)
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.configuration_model([1,2,3,4]*1000)
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, 
                                                                tmax = 10)
    
    '''
    Nk, Sk0, Ik0, Rk0 = _get_Nk_and_IC_as_arrays_(G, 
                                                initial_infecteds = initial_infecteds, 
                                                initial_recovereds = initial_recovereds,
                                                rho=rho, SIR=True)
    
    return SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin, tmax, 
                                        tcount, return_full_data=False)


#######      HETEROGENEOUS PAIRWISE
def _dSIS_heterogeneous_pairwise_(X, t, Nk, NkNl, tau, gamma, Ks):
    '''
    Gives the derivatives


    [\dot{Sk}] = gamma [Ik] - tau [SkI]
    [\dot{Ik}] = tau [SkI] - gamma [Ik]   = - [\dot{Sk}]
    [\dot{SkIl}] = gamma([IkIl] - [SkIl]) + tau([SkSlI] - [ISkIl] 
                                                 - [SkIl])
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


    :Arguments: 

        X : current values of variables
        t : current time
        Nk : The number of nodes of each degree; Nk[i] is the number of 
            nodes with degree Ks[i]
        NkNl : number of edges of various types.  NkNl[i,j] corresponds to 
            Ks[i] and Ks[j].
    **tau** transmission rate
    **gamma** recovery rate
    **Ks** a scipy array --- gives the observed degrees in increasing 
            order.
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
    kxSk = Ks*tmpSk  # k [S_k]
    lxSl = Ls*tmpSk  # l [S_l]
    
    kxSk[kxSk==0] = 1#avoid division by 0. This is okay since those places
    lxSl[lxSl==0] = 1#where it is 0 will have 0 numerator, and the whole
                     #expression should evaluate to 0.
    SkSlI = SkSl * ((Ls-1))*SlI/ (lxSl)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (kxSk)).T
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

    #the tmp variables defined below are used because in some places
    #the denominator becomes zero, but in those places the numerator
    #also becomes zero in such a way that the ratio should be 0.
    #making the terms in the denominator become 1 make the division go okay.
    
    tmpSk = 1.*Sk #the 1.* is to make it a copy, not the same object
    tmpSk[tmpSk==0] = 1  
    tmpSl = tmpSk
    
    tmpKs = 1.*Ks
    tmpKs[tmpKs==0] = 1
    tmpLs = tmpKs

    
    SkSlI = SkSl * ((Ls-1))*SlI/ (tmpLs*tmpSl)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (tmpKs*tmpSk)).T
    ISkSl = SkSlI.T 
    IkSlI = ISkIl.T 

    dSk = -tau*SkI
    dIk = tau*SkI - gamma*Ik
    dSkIl = -gamma*SkIl + tau*(SkSlI - ISkIl - SkIl)
    dSkSl = -tau*(SkSlI + ISkSl)

    dSkIl.shape=(kcount**2,1)
    dSkSl.shape = (kcount**2,1)

    dX = scipy.concatenate((dSk[:,None], dIk[:,None], 
                                dSkSl, dSkIl),axis=0).T[0]
    return dX



def SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, gamma, 
                                tmin = 0, tmax=100, tcount=1001, 
                                return_full_data = False, Ks = None):
                                
    r'''Encodes System (5.13) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    heterogeneous mean-field model closed at the level of triples

    :Arguments: 
    **Sk0** array.  
        (if Ks is defined, the definition changes slightly, see below)
        
        Sk0[k] is the number of nodes that are susceptible and have degree 
        k.  If one is empty, it becomes 0.

    **Ik0** array
        (if Ks is defined, the definition changes slightly, see below)
        
        similar to Sk0, but for infected.
    
    **SkSl0**  2D numpy array
        (if Ks is defined, the definition changes slightly, see below)

        SkSl0[k][l] is [S_kS_l] at 0
        see below for constraints these should satisfy related to 
        Sk0 and Ik0.  
        The code does not enforce these constraints.

    **SkIl0**  2D numpy array
        as in SkSl0
        
    **IkIl0**  2D numpy array
        as in SkSl0
        
    **tau** positive float
        transmission rate

    **gamma** number
        recovery rate
    
    **tmin**  number (default 0)
        minimum report time

    **tmax**  number (default 100)
        maximum report time 

    **tcount**  integer (default 1001)
        number of reports

    **return_full_data**  boolean (default False)
        If True, return times, Sk, Ik, SkIl, SkSl, IkIl
        If False, return times, S, I

    **Ks** scipy array. (default None)
        (helps prevent memory errors) if some degrees are not
        observed, then the corresponding entries of these arrays are
        zero.  This can lead to memory errors in the case of a
        network with many missing degrees.  So Ks is an (assumed)
        ordered vector stating which Ks are actually observed.  Then
        the Sk0[i] is the number of nodes that are susceptible and
        have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.

    In principle, there are constraints relating Sk with SkSl and SkIl 
    and similarly relating Ik with IkIl and SkIl.T.  
    
    No attempt is made to enforce these.  
    
    It is assumed the user will ensure acceptible inputs.

    We could also remove Sk0 and Ik0 as inputs and infer them from the 
    others, but for consistency with elsewhere, this is not done here.
    
    
    :Returns: 
        
    if return_full_data is True:
        returns **times, S, I, Sk, Ik, SkIl, SkSl, IkIl**
    if return_full_data is False:
        returns **times, S, I**

    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        Sk0 = 100 * scipy.ones(4)
        Ik0 = scipy.zeros(4)
        Ik0[3]=1
        SkSl0 = scipy.matrix([[0, 0,0,0],[0,100,0,0],[0,0,200,0],[0,0,0,294]])
        #only interact within a degree class, so the deg 1 and 2 are safe.
        SkIl0 = scipy.zeros((4,4))
        SkIl0[3,3] = 3
        IkIl0 = scipy.zeros((4,4))
        tau = 1
        gamma = 1
        
        t, S, I = EoN.SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, 
                                                    gamma)
    
    '''

    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
    times = scipy.linspace(tmin,tmax,tcount)
    Nk = Sk0+Ik0
    kcount = len(Nk)
    NkNl = SkSl0 + SkIl0 + IkIl0 + SkIl0.T
    SkSl0 = SkSl0.copy()
    SkIl0 = SkIl0.copy()
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)

    X0 = scipy.concatenate((Sk0[:,None], SkSl0, SkIl0), axis=0).T[0]

    X = _my_odeint_(_dSIS_heterogeneous_pairwise_, X0, times, 
                    args = (Nk, NkNl, tau, gamma, Ks))

    kcount = len(Nk)
    Sk = X.T[:kcount]
    #print Sk.size
    S = Sk.sum(axis=0)
    #print S.size
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

def SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma, 
                                tmin = 0, tmax=100, tcount=1001, 
                                return_full_data=False, Ks = None):
    '''Encodes System (5.15) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of triples"

    [\dot{S}_k] = -tau [S_k I]
    [\dot{I}_k] = tau [S_k I] - gamma [I_k]
    [\dot{R}_k] = gamma [I_k]  (but using Rk=Nk-Sk-Ik for this equation)
    [\dot{S_kI_l}] = -gamma[S_k I_l] + tau([S_k S_l I] - [I S_k I_l] 
                                           - [S_k I_l])
    [\dot{S_kS_l}] = -tau([S_k S_l I] + [I S_k S_l])

    [A_l S_k I] = ((k-1)/k) [A_l S_k] [S_k I]/ [S_k]
    [I S_k A_l] = ((k-1)/k) [I S_k] [S_k A_l]/ [S_k]

    :Arguments: 

    **Sk0** scipy array, 
        (if Ks is defined, the definition changes slightly, see below)
        
        Sk0[k] is number of degree k susceptible at 
        time 0.
        
    **Ik0** scipy array
        (if Ks is defined, the definition changes slightly, see below)
        
        as in Sk0

    **Rk0** scipy array
        (if Ks is defined, the definition changes slightly, see below)
        
        as in Sk0
        
    **SkSl0** scipy 2D array
        (if Ks is defined, the definition changes slightly, see below)
        
        SkSl0[k][l] is [S_kS_l] at 0

    **SkIl0** scipy 2D array
        (if Ks is defined, the definition changes slightly, see below)
        
        as in SkSl0
        
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        If True, return times, Sk, Ik, Rk, SkIl, SkSl
        If False, return times, S, I, R
    **Ks** scipy array. (default None)
        (helps prevent memory errors) if some degrees are not
        observed, then the corresponding entries of these arrays are
        zero.  This can lead to memory errors in the case of a
        network with many missing degrees.  So Ks is an (assumed)
        ordered vector stating which Ks are actually observed.  Then
        the Sk0[i] is the number of nodes that are susceptible and
        have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.        

    :Returns: 

    if return_full_data is True
        returns **times, S, I, R, Sk, Ik, Rk, SkIl, SkSl**
    if return_full_data is False
        return **times, S, I, R**
    
    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        
'''
    
    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
#    print "Ks is ", Ks
    times = scipy.linspace(tmin,tmax,tcount)

    Nk = Sk0+Ik0+Rk0
    kcount = len(Ks)
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)

    X0 = scipy.concatenate((Sk0[:,None], Ik0[:,None], SkSl0, SkIl0), 
                                axis=0).T[0]
    X = integrate.odeint(_dSIR_heterogeneous_pairwise_, X0, times, 
                            args = (tau, gamma, Nk, Ks))

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


def SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, initial_infecteds = None,
                                            rho = None, tmin = 0, 
                                            tmax=100, tcount=1001, 
                                            return_full_data = False):
    r'''
    Calls SIS_heterogeneous_pairwise after calculating Sk0, Ik0, SkSl0, SkIl0, IkIl0
    from a graph G and initial fraction infected rho.

    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        tells whether to just return times, S, I, or all calculated data.
        if True, then returns times, S, I, SI, SS

    
    :Returns: 
    if return_full_data is True:
        returns **times, S, I, Sk, Ik, SkIl, SkSl, IkIl**
    if return_full_data is False:
        returns **times, S, I**
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.fast_gnp_random_graph(10000,0.0005)
        tau = 1
        gamma = 3
        rho = 0.02
        t, S, I = EoN.SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho, tmax = 20)
    
    :WARNING:
        
    ::
    
        This can have segmentation faults if there are too many degrees in the
        graph.  This appears to happen because of trouble in numpy, and I have
        not been able to find a way around it.
    
    '''
    
    Nk, Sk0, Ik0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds = initial_infecteds,
                                            rho=rho, SIR=False)

    NkNl, SkSl0, SkIl0, IkIl0, Ks = \
                _get_NkNl_and_IC_as_arrays_(G, initial_infecteds = initial_infecteds,
                                            rho=rho, withKs = True, SIR=False)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    return SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, 
                                        gamma, tmin, tmax, tcount, 
                                        return_full_data, 
                                        Ks = scipy.array(Ks))

    
def SIR_heterogeneous_pairwise_from_graph(G, tau, gamma, initial_infecteds=None, 
                                        initial_recovereds = None, rho = None, tmin=0, 
                                            tmax=100, tcount = 1001, 
                                            return_full_data=False):
    r'''Calls SIR_heterogeneous_pairwise after calculating Sk0, Ik0, Rk0, SkSl0, SkIl0
    from a graph G and initial fraction infected rho. 
    

    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        If True, return times, Sk, Ik, Rk, SkIl, SkSl
        If False, return times, S, I, R

    :Returns: 

    if return_full_data is True
        returns **times, S, I, R, Sk, Ik, Rk, SkIl, SkSl**
    if return_full_data is False
        return **times, S, I, R**
        
    :WARNING:
    
        This can have segmentation faults if there are too many degrees in the
        graph.  This appears to happen because of trouble in numpy, and I have
        not been able to find a way around it.

    '''
    
    Nk, Sk0, Ik0, Rk0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds=initial_infecteds, 
                                        initial_recovereds = initial_recovereds, rho=rho, SIR=True)
    NkNl, SkSl0, SkIl0, Ks = _get_NkNl_and_IC_as_arrays_(G, initial_infecteds=initial_infecteds, 
                                        initial_recovereds = initial_recovereds, rho=rho, withKs = True, 
                                                        SIR = True)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    Rk0 = scipy.array([Rk0[k] for k in Ks])
    return SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma,
                                        tmin = tmin, tmax=tmax, tcount=tcount,
                                        return_full_data=return_full_data,
                                        Ks = Ks)




#######    COMPACT PAIRWISE


def _dSIS_compact_pairwise_(X, t, Nk, twoM, tau, gamma):
    Sk = X[:-2]
    SI, SS= X[-2:]

    Ik = Nk - Sk
    II = twoM - SS - 2*SI

    ks = scipy.arange(len(Sk))
    SX = ks.dot(Sk)#float(SI + SS)
    Q = (1./SX**2) * (ks*(ks-1)).dot(Sk)

    dSk = gamma *Ik - tau * ks *Sk *SI/SX
    dSI = gamma*(II-SI) + tau*(SS-SI)*(SI)*Q - tau*SI
    dSS = 2*gamma*SI - 2*tau*SS*SI*Q

    dX = scipy.concatenate((dSk, [dSI, dSS]), axis=0)
    return dX

def _dSIR_compact_pairwise_(X, t, N, tau, gamma):
    Sk = X[:-3]
    SS, SI, R = X[-3:]
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


def SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin = 0, 
                            tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.18) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [dot{S}_k] = gamma[I_k] - tau k [S_k] [SI]/[SX]
    
    [dot{I}_k] = tau * k *[S_k] [SI]/SX - gamma [I_k] = -[dot{S}_k]

    [dot{SI}] = gamma([II]-[SI]) + tau([SS]-[SI])[SI]Q - tau[SI]

    [dot{SS}] = 2 gamma[SI] - 2 tau [SS] [SI] Q

    [dot{II}] = 2 tau[SI] = 2 gamma[II] + 2 tau[SI]^2Q

    [SX] = sum_k k [S_k]

    Q = (1/[SX]^2) sum_k (k-1)k[S_k]

    :conserved quantities:  
    [Sk]+[Ik]
    
    SS + II + 2SI

    :Arguments: 
    **Sk0** scipy array
        number susceptible for each k
    **Ik0** scipy array
        number infected for each k
    **SI0** number
        number of SI edges
    **SS0** number
        number of SS edges
    **II0** number
        number of II edges
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        if True, return times, S, I, Sk, Ik, SI, SS, II
        if False,  return times, S, I
    
    
    :Returns: 
    
    All scipy arrays
    
    if return_full_data
        return **times, S, I, Sk, Ik, SI, SS, II**
    else
        return **times, S, I**


    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        
    '''
    Nk = Sk0+Ik0
    twoM = SS0+II0+2*SI0 #each of the M edges counted twice
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.concatenate((Sk0, [SI0,SS0]), axis=0)
    X = integrate.odeint(_dSIS_compact_pairwise_, X0, times, 
                            args = (Nk, twoM, tau, gamma))

    Sk = X.T[:-2]
    S = Sk.sum(axis=0)

    Ik = Nk[:,None] - Sk
    I = Ik.sum(axis=0)


    if return_full_data:
        SI, SS = X.T[-2:]
        II = twoM - SS - 2*SI
        return times, S, I, Sk, Ik, SI, SS, II
    else:
        return times, S, I

def SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=0, tmax=100,
                            tcount=1001, return_full_data=False):
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

    :Arguments: 

    **Sk0** scipy array
        initial number of suscetibles of each degree k
    **I0** number
        initial number infected
    **R0** number
        initial number recovered
    **SS0** number
        initial number of SS edges
    **SI0** number
        initial number of SI edges
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
    
    if return_full_data:
       **times, Sk, I, R, SS, SI**
    else:
       **times, S, I, R**


    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
    '''
    
    times = scipy.linspace(tmin,tmax,tcount)
    N = I0+R0+sum(Sk0)
    X0 = scipy.concatenate((Sk0, [SS0, SI0, R0]), axis=0)
    X = integrate.odeint(_dSIR_compact_pairwise_, X0, times, 
                            args = (N, tau, gamma))
    SI, SS, R = X.T[-3:]
    Sk = X.T[:-3]
    S = Sk.sum(axis=0)
    I = N - R - S
    if return_full_data:
        return times, Sk, I, R, SS, SI
    else:
        return times, S, I, R


def SIS_compact_pairwise_from_graph(G, tau, gamma, initial_infecteds=None, rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    r'''Calls SIS_compact_pairwise after calculating Sk0, Ik0, SI0, SS0, II0
    from the graph G and initial fraction infected rho.
    
    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean (default False)
        if True, return times, S, I, Sk, Ik, SI, SS, II
        if False,  return times, S, I
    
    
    :Returns: 
    
    All scipy arrays
    if return_full_data:
        return **times, S, I, Sk, Ik, SI, SS, II**
    else:
        return **times, S, I**
    
    '''
    
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is None and initial_infecteds is None:
        rho = 1./G.order()

    Nk, Sk0, Ik0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds=initial_infecteds, rho=rho, SIR=False)
    
    if initial_infecteds is not None:
        SS0, SI0, II0 = _count_edge_types_(G, initial_infecteds, SIR=False)   
    else:
        maxk = max(dict(G.degree()).values())
        SS0 = sum(Nk[k]*k*(1-rho)*(1-rho) for k in range(maxk+1))
        SI0 = sum(Nk[k]*k*(1-rho)*rho for k in range(maxk+1))
        II0 = sum(Nk[k]*k*rho*rho for k in range(maxk+1))
    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, 
                                    tmax, tcount, return_full_data)

    
def SIR_compact_pairwise_from_graph(G, tau, gamma,  initial_infecteds=None, 
                                    initial_recovereds = None, rho = None, 
                                    tmin = 0, tmax=100, tcount=1001,
                                    return_full_data=False):
    r'''Calls SIR_compact_pairwise after calculating Sk0, I0, R0, SS0, SI0
    from the graph G and initial fraction infected rho.
    
    
    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
    
    if return_full_data:
       **times, Sk, I, R, SS, SI**
    else:
       **times, S, I, R**
    '''
    
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is None and initial_infecteds is None:
        rho = 1./G.order()

    Nk, Sk0, Ik0, Rk0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds=initial_infecteds, 
                                                initial_recovereds=initial_recovereds, 
                                                rho=rho, SIR=True)
    
    I0 = sum(Ik0)
    R0 = sum(Rk0)


    if initial_infecteds is not None:
        SS0, SI0 = _count_edge_types_(G, initial_infecteds, initial_recovereds, SIR=True)   
    else: 
        ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]
        SX0 = scipy.dot(Sk0,ks)
        SS0 = (1-rho)*SX0
        SI0 = rho*SX0
    
    return SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=tmin,
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)




#######SUPER COMPACT PAIRWISE

def _dSIS_super_compact_pairwise_(X, t, tau, gamma, N, k_ave, ksquare_ave, 
                                    kcube_ave):
    '''    
    [\dot{I}] = tau [SI] - gamma [I]
    [\dot{SS}] = 2 gamma [SI] - 2 tau [SI] [SS] Q
    [\dot{SI}] = gamma ([II]-[SI]) + tau [SI] ([SS]-[SI])Q - tau [SI]
    [\dot{II}] = -2 gamma [II] + 2 tau [SI]^2 Q + 2 tau [SI]
    
    Q = ((<K^2>(<K^2>-n_S<K>) + <K^3>(n_S-<K>))
        /(n_S(<K^2>-<K>^2)) - 1)/n_S[S]
        
        
    n_S = ([SI] + [SS])/(N-[I])
    '''

    I, SS, SI, II = X
    S = N-I
    
    n_S = (SS+SI)/(S)
    
    Q = ((ksquare_ave*(ksquare_ave-n_S*k_ave) \
            + kcube_ave*(n_S-k_ave))/(n_S*(ksquare_ave-k_ave**2))-1)/(S*n_S)
    dIdt = tau*SI - gamma*I
    dSSdt = 2*gamma*SI - 2*tau*SI*SS*Q
    dSIdt = gamma*(II - SI) + tau*SI*(SS-SI)*Q - tau*SI
    dIIdt = -2*gamma*II + 2*tau*SI**2*Q + 2*tau*SI

    dX = scipy.array([dIdt, dSSdt, dSIdt, dIIdt])
    return dX

def _dSIR_super_compact_pairwise_(X, t, tau, gamma, psihat, psihatPrime, 
                                    psihatDPrime, N):
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


def SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, k_ave, 
                                ksquare_ave, kcube_ave, tmin = 0, tmax=100, 
                                tcount=1001, return_full_data=False):
    '''
    Encodes system (5.20) of Kiss, Miller, & Simon.  Please cite the 
    book if using this algorithm.

    :Arguments: 

    **S0** number
        initial number susceptible
    **I0** number
        initial number infected
    **SS0** number
        initial number of susceptible-susceptible edges
    **SI0** number
        initial number of susceptible-infected edges
    **II0** number
        initial number of infected-infected edges.
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **k_ave** number
        average value of degree k
    **ksquare_ave** number
        average value of k**2
    **kcube_ave** number
        average value of k**3
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all 
        calculated data.

    :Returns: 
        
    if return_full_data is True
        returns **times, S, I, SS, SI, II**
    if return_full_data is False
        returns **times, S, I**
        
    :SAMPLE USE:

    ::
    
        import networkx as nx
        import EoN
    
    '''
    X0 = [I0, SS0, SI0, II0]
    N = S0+I0
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_super_compact_pairwise_, X0, times, 
                            args = (tau, gamma, N, k_ave, ksquare_ave, 
                                    kcube_ave))
    I, SS, SI, II = X.T
    S = N-I
    if return_full_data:
        return times, S, I, SS, SI, II
    else:
        return times, S, I



def SIR_super_compact_pairwise(R0, SS0, SI0,  N, tau, gamma, psihat, 
                                psihatPrime, psihatDPrime, tmin = 0, 
                                tmax = 100, tcount = 1001, 
                                return_full_data = False):
    '''
    Encodes system (5.22) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{theta} = -tau [SI]/N*psihat(theta)

    [dot{SS}] = -2 tau [SS] [SI] Q

    [dot{SI}] = -gamma[SI] + tau([SS]-[SI])[SI]Q - tau*[SI]

    [dot{R}] = gamma*[I]

    [S] = N psihat(theta)

    [I] = N-[S]-[R]

    Q = psihat_xx(theta)/N(psihat_x(theta))^2


    
    :Arguments: 
    **R0** number
        initial number of R nodes.
    **SS0** number
        initial number of SS edges
    **SI0** number
        initial number of SI edges
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all 
        calculated data.

    :Returns: 
        
    if return_full_data:
        return **times, S, I, R, SS, SI**
    else:
        return **times, S, I, R**

'''
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.array([1., SS0, SI0, R0])
    X = integrate.odeint(_dSIR_super_compact_pairwise_, X0, times, 
                            args = (tau, gamma, psihat, psihatPrime, 
                                    psihatDPrime, N))
    theta, SS, SI, R = X.T
    S = N*psihat(theta)
    I = N-S-R
    if return_full_data:
        return times, S, I, R, SS, SI
    else:
        return times, S, I, R

def SIS_super_compact_pairwise_from_graph(G, tau, gamma, initial_infecteds=None, 
                                            rho = None, tmin = 0,
                                            tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''Calls SIS_super_compact_pairwise after calculating S0, I0, SS0, SI0, II0
    from the graph G and initial fraction infected rho'''

    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    
    Nk, Sk0, Ik0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds = initial_infecteds, 
                                            rho=rho, SIR=False)

    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]

    S0 = sum(Sk0)
    I0 = sum(Ik0)

    if initial_infecteds is not None:
        SS0, SI0, II0 = _count_edge_types_(G, initial_infecteds, SIR=False)   
    else:    
        if rho is None:
            rho = 1./G.order()
    

        SX0 = scipy.dot(Sk0,ks)
        SS0 = (1-rho)*SX0
        SI0 = rho*SX0
        II0 = scipy.dot(Nk,ks)-SX0
        
    Pk = get_Pk(G)
    Pks = scipy.array([Pk.get(k,0) for k in ks])
    k_ave = scipy.dot(Pks, ks)
    ksquare_ave = scipy.dot(Pks, ks*ks)
    kcube_ave = scipy.dot(Pks, ks*ks*ks)
    
    return SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, 
                                        k_ave, ksquare_ave, kcube_ave, 
                                        tmin=tmin, tmax=tmax, tcount=tcount, 
                                        return_full_data=return_full_data)

def SIR_super_compact_pairwise_from_graph(G, tau, gamma,  initial_infecteds=None, 
                                            initial_recovereds = None, 
                                            rho = None, tmin = 0,
                                            tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''Calls SIR_super_compact_pairwise after calculating R0, SS0, SI0
    from the graph G and initial fraction infected rho
    
    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or all 
        calculated data.

    :Returns: 
        
    if return_full_data:
        return **times, S, I, R, SS, SI**
    else:
        return **times, S, I, R**
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is None and initial_infecteds is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = _get_Nk_and_IC_as_arrays_(G, initial_infecteds = initial_infecteds, 
                                                initial_recovereds = initial_recovereds,
                                                rho=rho, SIR=True)

    Pk = get_Pk(G)
    
    N = G.order()
    
    R0 = sum(Rk0)

    if initial_infecteds is not None:
        SS0, SI0 = _count_edge_types_(G, initial_infecteds, initial_recovereds, SIR=True)
        def psihat(x): #probably faster if vectorized, 
                   #but need to be careful with broadcasting...
            return sum((Sk0[k]*(x**k)) for k in Pk)/N
        def psihatPrime(x):
            return sum(k*Sk0[k]*x**(k-1) for k in Pk)/N
        def psihatDPrime(x):
            return sum(k*(k-1)*Sk0[k]*x**(k-2) for k in Pk)/N

    else:       
        ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]
        SX0 = scipy.dot(Sk0,ks)
        SS0 = (1-rho)*SX0
        SI0 = rho*SX0

        
        def psihat(x): #probably faster if vectorized, 
                   #but need to be careful with broadcasting...
            return (1-rho)*sum(Pk[k]*x**k for k in Pk)
        def psihatPrime(x):
            return (1-rho)*sum(k*Pk[k]*x**(k-1) for k in Pk)
        def psihatDPrime(x):
            return (1-rho)*sum(k*(k-1)*Pk[k]*x**(k-2) for k in Pk)

    return  SIR_super_compact_pairwise(R0, SS0, SI0, N, tau, gamma, psihat, 
                                        psihatPrime, psihatDPrime, 
                                        tmin = tmin, tmax = tmax, 
                                        tcount = tcount, 
                                        return_full_data = return_full_data)






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
    #note that dSIR_effective_degree has some commands for vectorizing 
    #tentatively coded (and commented out)
    ksq = original_shape[0]*original_shape[1]
    Ssi = X[:ksq]
    Isi = X[ksq:]
    Ssi.shape = original_shape
    Isi.shape = original_shape

    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    ISI = sum([sum([i*(i-1)*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SI = sum([sum([i*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])

    g1 = scipy.zeros(original_shape)
    g2 = scipy.zeros(original_shape)
    t1 = scipy.zeros(original_shape)
    t2 = scipy.zeros(original_shape)
    
    dSsi = scipy.zeros(original_shape)
    dIsi = scipy.zeros(original_shape)
    for s in range(original_shape[0]):
        for i in range(original_shape[1]):
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
            
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*Isi[s,i] \
                        + gamma*((i+1)*Ssm1ip1 - i*Ssi[s,i]) \
                        + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
            dIsi[s,i] =  tau*i*Ssi[s,i] - gamma*Isi[s,i] \
                        + gamma*((i+1)*Ism1ip1 - i*Isi[s,i]) \
                        + tau*(ISI/SI + 1)*((s+1)*Isp1im1 - s*Isi[s,i])# 

    dSsi.shape = (original_shape[0]*original_shape[1])
    dIsi.shape = (original_shape[0]*original_shape[1])

    dX = scipy.concatenate((dSsi, dIsi), axis=0)
    return dX


def _dSIR_effective_degree_(X, t, N, original_shape, tau, gamma):
    R = X[-1]
    Ssi = X[:-1]
    Ssi.shape=original_shape
    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    
    #commenting out commands for vectorizing this.  
    #I should do this eventually, but not now.  Apply to SIS version as well.
    #ultimately I think right option is to pad with zeros and then 
    #cut out appropriate section.
    #ivec = scipy.array(range(original_shape[1]))
    #ivec.shape = (1,original_shape[1])
    #svec = scipy.array(range(original_shape[0]))
    #svec.shape = (original_shape[0],1)
    #
    #Ssip1 = Ssi[:,1:]
    #scipy.pad(Ssip1, pad_width=((0,0),(0,1)), mode = 'constant', 
    #          constant_values=0)
    #Ssp1im1 = Ssi[1:,:-1]
    #scipy.pad(Ssp1im1, pad_width=((0,1),(1,0)), mode = 'constant', 
    #           constant_values=0)
    #dSsi = - tau* ivec*Ssi + gamma*((ivec+1)*Ssip1 - i*Ssi) 
    #       + tau *ISS*((svec+1)*Sp1im1 - svec*Ssi)/SS

    dSsi = scipy.zeros(original_shape)
    for s in range(original_shape[0]):
        for i in range(original_shape[1]):
            if i+1 == original_shape[1]:
                Ssip1 = 0
            else:
                Ssip1 = Ssi[s,i+1]
            if s+1 == original_shape[0] or i == 0:
                Ssp1im1=0
            else:
                Ssp1im1 = Ssi[s+1,i-1]    
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*((i+1)*Ssip1 - i*Ssi[s,i]) \
                        + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
    S = Ssi.sum() 
    I = N-S-R
    dR = gamma*I

    dSsi.shape = (original_shape[0]*original_shape[1])
    dX = scipy.concatenate((dSsi, [dR]), axis=0)
    return dX



def SIS_effective_degree(Ssi0, Isi0, tau, gamma, tmin = 0, tmax=100, 
                            tcount=1001, return_full_data=False):
    '''Encodes system (5.36) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    :Arguments: 

        Ssi0 and Isi0 : (square) numpy 2D arrays of same shape.
            Entries are initial number susceptible or infected 
            with given initial number of susceptible/infected 
            neighbors.
    **tau** positive float
            transmission rate
    **gamma** number
            recovery rate
    **tmin**  number (default 0)
            minimum report time
    **tmax**  number (default 100)
            maximum report time 
    **tcount**  integer (default 1001)
            number of reports
    **return_full_data**  boolean (default False)
            tells whether to just return times, S, I, R or all calculated data.
            if True, 
                return times, S, I, Ssi, Isi
            if False, 
                return times, S, I
                   
   :Returns: 
       
    if return_full_data:
        return times, S, I, Ssi, Isi
    else:
        return times, S, I
    
    '''
    times = scipy.linspace(tmin,tmax,tcount) 
    original_shape = Ssi0.shape
    ksq = original_shape[0]*original_shape[1]
    Ssi0.shape = (1,ksq)
    Isi0.shape = (1,ksq)
    
    X0= scipy.concatenate((Ssi0[0], Isi0[0]), axis=0)
    X = integrate.odeint(_dSIS_effective_degree_, X0, times, 
                            args = (original_shape, tau, gamma))
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

def SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=0, tmax=100, 
                            tcount=1001, return_full_data=False):
    '''Encodes system (5.38) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_{s,i} = - tau i S_{s,i}  + gamma((i+1)S_{s,i+1} - i S_{s,i})
                    + tau [ISS]((s+1)S_{s+1,i-1} - sS_{s,i})/[SS]
    \dot{R} = gamma I
    S = \sum_{s,i} S_{s,i}
    I = N-S-R

    :Arguments: 

    **S_si0**  (square) numpy 2-D array
        S_{s,i} at time 0
    **I0** number
        number of infected individuals at time 0
    **R0** number
        number of recovered individuals at time 0
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or 
        all calculated data.

    :Returns: 
        
    if return_full_data==False
        **times** scipy.array of times
        
        **S** scipy.array of number susceptible
        
        **I** scipy.array of number infected
        
        **R** scipy.array of number recovered
    else
        **times** as before
    
        **S** number susceptible
            
        **I** number infected

        **R** number recovered

        **S_si** S_{s,i} at each time in times

    '''
    times = scipy.linspace(tmin,tmax, tcount)
    N = S_si0.sum()+I0+R0
    original_shape = S_si0.shape
    S_si0.shape = (original_shape[0]*original_shape[1]) 
    #note this makes it array([[values...]])
    R0=scipy.array([R0])
    R0.shape=(1)
    X0 = scipy.concatenate((S_si0, R0), axis=0)
    X = integrate.odeint(_dSIR_effective_degree_, X0, times, 
                            args = (N, original_shape, tau, gamma))

    R = X.T[-1]
    S_si = X.T[:-1]
    S = S_si.sum(axis=0)
    I = N - R - S
    if return_full_data:
        S_si.shape = (original_shape[0], original_shape[1], tcount)
        return  times, S, I, R, S_si
    else:
        return times, S, I, R


def SIS_effective_degree_from_graph(G, tau, gamma, initial_infecteds=None, 
                                    rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
                                    

    r'''Calls SIS_effective_degree after calculating Ssi0, Isi0 from
    the graph G and initialf fraction infected rho.
    
        
    :WARNING:
        
    ::
    
        This can have segmentation faults if there are too many degrees in the
        graph.  This appears to happen because of trouble in numpy, and I have
        not been able to find a way around it.
'''

    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")

    Nk = Counter(dict(G.degree()).values())
    maxk = max(Nk.keys())
    

    S_si0 = scipy.zeros((maxk+1,maxk+1))
    I_si0 = scipy.zeros((maxk+1,maxk+1))

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds)
        for node in G.nodes():
            s = sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
            i = G.degree(node)-s
            if status[node] == 'S':
                S_si0[s][i] += 1
            else:
                I_si0[s][i] += 1
    else:  
        if rho is None:
            rho = 1./G.order()
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])       
        for s in range(maxk+1):
            for i in range(maxk+1-s):
                binomial_result = scipy.special.binom(s+i,i)
                if binomial_result < float('Inf'):
                    #sometimes sp.special.binom() returns 'inf', this tries to avoid those cases
                    S_si0[s,i] = (1-rho)*Nk[s+i] * binomial_result \
                                * (rho**i) * (1-rho)**s
                    I_si0[s,i] = rho*Nk[s+i] * binomial_result \
                                * (rho**i) * (1-rho)**s
                else:
                    #implicitly assuming that those cases where it return 'Inf'
                    #have rho**i *(1-rho)**s is small enough to make this
                    #effectively 0.
                    S_si0[s,i] = 0
                    I_si0[s,i] = 0

    return SIS_effective_degree(S_si0, I_si0, tau, gamma, tmin = tmin, 
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)

def SIR_effective_degree_from_graph(G, tau, gamma, initial_infecteds=None, 
                                    initial_recovereds = None, rho = None, 
                                    tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    r'''Calls SIR_effective_degree after calculating S_si0, I0, R0 from the
    graph G and initial fraction infected rho
    
    
    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports    
        
    **return_full_data**  boolean
        tells whether to just return times, S, I, R or 
        all calculated data.

    :Returns: 
        
    if return_full_data==False
        **times** scipy.array of times
        
        **S** scipy.array of number susceptible
        
        **I** scipy.array of number infected
        
        **R** scipy.array of number recovered
    else
        **times** as before
    
        **S** number susceptible
            
        **I** number infected

        **R** number recovered

        **S_si** S_{s,i} at each time in times
        
    :WARNING:
        
    
    This can have segmentation faults if there are too many degrees in the
    graph.  This appears to happen because of trouble in numpy, and I have
    not been able to find a way around it.
'''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    

    Nk = Counter(dict(G.degree()).values())
    maxk = max(Nk.keys())
    S_si0 = scipy.zeros((maxk+1,maxk+1))
    I0 = 0
    R0 = 0

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds)
        for node in G.nodes():
            s = sum(1 for nbr in G.neighbors(node) if status[nbr] == 'S')
            i = sum(1 for nbr in G.neighbors(node) if status[nbr] == 'I')
            if status[node] == 'S':
                S_si0[s][i] += 1
            elif status[node] == 'I':
                I0 += 1
            else:
                R0 += 1
    else:  
        if rho is None:
            rho = 1./G.order()
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])       
        for s in range(maxk+1):
            for i in range(maxk+1-s):
                binomial_result = scipy.special.binom(s+i,i) 
                if binomial_result < float('Inf'):
                    S_si0[s,i] = (1-rho)*Nk[s+i] * binomial_result \
                                * (rho**i) * (1-rho)**s
                else:
                    #implicitly assuming that those cases where it return 'Inf'
                    #have rho**i *(1-rho)**s is small enough to make this
                    #effectively 0.
                    S_si0[s,i] = 0


        I0 = rho*sum(Nk)
        R0=0

    return SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=tmin, 
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)




#######     COMPACT EFFECTIVE DEGREE


def SIS_compact_effective_degree(Sk0, Ik0, SI0, SS0, II0, tau, gamma, 
                                    tmin = 0, tmax=100, tcount=1001, 
                                    return_full_data=False):
    r'''
    Encodes system (5.44) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    This model is identical to the SIS compact pairwise model, so it 
    simply calls SIS_compact_pairwise()'''

    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, 
                                    tmax, tcount, return_full_data)

def SIS_compact_effective_degree_from_graph(G, tau, gamma, initial_infecteds = None,
                                            rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''because the SIS compact effective degree model is identical to the
    compact pairwise model, simply calls SIS_compact_pairwise_from_graph'''
                                            
    return SIS_compact_pairwise_from_graph(G, tau, gamma, 
                                            initial_infecteds = initial_infecteds, 
                                            rho = rho, tmin = tmin, tmax=tmax, 
                                            tcount=tcount, 
                                            return_full_data=return_full_data)


def _dSIR_compact_effective_degree_(X, t, N, tau, gamma):
    Skappa = X[:-2]
    R, SI = X[-2:]
    I = N- R- Skappa.sum()
    kappas = scipy.arange(len(Skappa))
    effectiveI = float(SI) /Skappa.dot(kappas)
    dSkappa = effectiveI*(-(tau+gamma)*kappas*Skappa \
                + gamma*shift(kappas*Skappa,-1))
    dSI = -(tau+gamma)*SI \
            + tau*(effectiveI-2*effectiveI**2)*sum(kappas*(kappas-1)*Skappa)

    dR = gamma*I
    dX = scipy.concatenate((dSkappa, [dR, dSI]), axis=0) 
    return dX
    
def SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, tmin=0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    '''
    Encodes system (5.43) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_kappa = <I> [-(tau+gamma) kappa S_kappa 
                         + gamma(kappa+1)S_{kappa+1}
    [\dot{SI}] = -(tau+gamma)[SI] 
                + tau(<I> - 2 <I>^2) sum_{kappa} kappa(kappa-1) S_kappa
    \dot{R} = gamma I
    <I> = [SI]/sum_kappa kappa S_kappa
    S = sum_kappa S_kappa
    I = N - S - R

    :Arguments: 

        Skappa0 : scipy array
            from S_0(0) up to S_kappamax(0) of number susceptible with 
            each effective degree
            Skappa = number of nodes that are susceptible and have 
            kappa non-recovered neighbors
    **I0** number
            number of infected individuals at time 0
    **R0** number
            initial number recovered
    **SI0** number
            initial number of SI edges
    **tau** positive float
            transmission rate
    **gamma** number
            recovery rate
    **tmin**  number (default 0)
            minimum report time
    **tmax**  number (default 100)
            maximum report time 
    **tcount**  integer (default 1001)
             number of reports
    **return_full_data**  boolean
            tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
        
    if return_full_data==False
        **times** scipy.array of times
        
        **S** scipy.array of number susceptible
        
        **I** scipy.array of number infected
        
        **R** scipy.array of number recovered
    else
        **times** as before
    
        **S** number susceptible
            
        **I** number infected

        **R** number recovered

        **SI** S_{s,i} 
            number of SI edges
    '''

    N = Skappa0.sum()+I0+R0
    X0= scipy.concatenate((Skappa0,[R0,SI0]), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_compact_effective_degree_, X0, times, 
                            args = (N, tau, gamma))
    Skappa = X.T[:-2]
    S = Skappa.sum(axis=0)
    R, SI = X.T[-2:]
    I = N - S -R
    if return_full_data:
        return times, S, I, R, Skappa, SI
    else:
        return times, S, I, R

def SIR_compact_effective_degree_from_graph(G, tau, gamma, initial_infecteds=None, 
                                        initial_recovereds = None, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''Calls SIR_compact_effective_degree after calculating Skappa0, I0, R0, SI0
    from the graph G and initial fraction infected rho.
    :Arguments: 

    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmin**  number (default 0)
        minimum report time
    **tmax**  number (default 100)
        maximum report time 
    **tcount**  integer (default 1001)
        number of reports    
    **return_full_data**  boolean
            tells whether to just return times, S, I, R or all calculated data.

    :Returns: 
        
    if return_full_data==False
        **times** scipy.array of times
        
        **S** scipy.array of number susceptible
        
        **I** scipy.array of number infected
        
        **R** scipy.array of number recovered
    else
        **times** as before
    
        **S** number susceptible
            
        **I** number infected

        **R** number recovered

        **SI** S_{s,i} 
            number of SI edges
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    
    
    if initial_infecteds is not None:
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Skappa0 = scipy.zeros(maxk+1)
        I0=0
        R0=0
        SI0=0
        status = _initialize_node_status_(G, initial_infecteds, 
                                    initial_recovereds = initial_recovereds)
                                    
        for node in G.nodes():
            if status[node] == 'S':
                kappa = sum(1 for nbr in G.neighbors(node) if status[nbr] !='R')
                Skappa0[kappa] += 1
                SI0 += sum(1 for nbr in G.neighbors(node) if status[nbr] == 'I')
            elif status[node] == 'I':
                I0 += 1
            else: #status[node]=='R'
                R0 += 1            
    else:
        if rho is None:
            rho = 1./G.order()
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])
        Skappa0 = Nk*(1-rho)  #Skappa0 = Sk0
        I0 = rho*sum(Nk)
        R0=0
        SI0 = sum([k*Skappa0[k]*rho for k in range(maxk+1)])
    return SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, 
                                        tmin=tmin, tmax=tmax, tcount=tcount, 
                                        return_full_data=return_full_data)






#######################################
#    EBCM and other related results   #
#######################################

def Epi_Prob_discrete(Pk, p, number_its = 100):
    '''Encodes System (6.2) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    :Arguments: 

    **Pk**  dict
            Pk[k] is probability a node has degree k.

    **p**  number
        transmission probability

    **number_its**    int    
        number of iterations before assumed converged.
        default value is 100

    :Returns: 
    **PE**  float
        Calculated Epidemic probability in the limit of a negligible
        fraction initially infected (assuming configuration model)
    '''
    psi = get_PGF(Pk)
    psiPrime = get_PGFPrime(Pk)
    
    alpha = 1-p
    k_ave = psiPrime(1.)
    for counter in range(number_its):
        alpha = 1-p +p *psiPrime(alpha)/k_ave
    return 1- psi(alpha)


def Epi_Prob_cts_time(Pk, tau, gamma, umin=0, umax = 10, ucount = 1001, 
                        number_its = 100):
    r'''Encodes System (6.3) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The equations are rescaled by setting $u=\gamma T$.  Then it becomes

    P = 1- \int_0^\infty \psi(\alpha(u/\gamma)) e^{-u} du
    alpha_d(u/\gamma) = 1- p(u/\gamma)
                        + p(u/\gamma)
                          \int_0^\infty 
                                (\psiPrime(\alpha(\hat{u}/\gamma))/<K>)
                                e^{-u}du

    where p(u/\gamma) = 1 - e^{-\tau u/\gamma}

    Define 
        \hat{p}(u) = p(u/\gamma), and \hat{\alpha}(u) = \alpha(u/\gamma)
    and then drop hats to get

    P = 1-\int_0^\infty \psi(\alpha(u)) e^{-u} du
    \alpha(u) = 1-p(u) + p(u) 
                         \int_0^\infty 
                              (\psiPrime(\alpha(u))/<K>)e^{-u} du

    with initial guess 
        \alpha_1(u) = e^{-\tau u/\gamma} 
    and 
        p(u) = 1-e^{-\tau u/\gamma}
    
    :Arguments: 

    **Pk** dict
        Pk[k] is probability a node has degree k.

    **tau** float
        transmission rate

    **gamma** float
        recovery rate

    **umin** minimal value of \gamma T used in calculation
    **umax** maximum value of \gamma T used in calculation
    **ucount** number of points taken for integral.
             So this integrates from umin to umax using simple Riemann 
             sum.

    **number_its**    int    
        number of iterations before assumed converged.
        default value is 100

    :Returns: 
    
    **PE**  float
            Calculated Epidemic probability (assuming configuration model)
    '''
    psi = get_PGF(Pk)
    psiPrime = get_PGFPrime(Pk)

    us = scipy.linspace(umin, umax, ucount) 
    alpha = scipy.e**(-tau*us/gamma)  #initial guess for alpha(u)
    p = 1- scipy.e**(-tau*us/gamma)    #initial guess for p(u)
    exp_neg_u = scipy.e**(-us)        #e^{-u}
    for counter in range(number_its):
        alpha = 1 - p + p* (psiPrime(alpha)/kave).dot(exp_neg_u)/(ucount-1.)
    return 1 - psi(alpha).dot(exp_neg_u)/(ucount - 1)


def Epi_Prob_non_Markovian(Pk, Pxidxi, po, number_its = 100):
    r'''Encodes system (6.5) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    :Arguments: 

    **Pk**  dict Pk[k] is probability a node has degree k.

    **Pxidxi** dict
        Pxidxi[xi] is  P(xi)dxi for user-selected xi.  The 
        algorithm will replace the integral with
        \sum_{xi \in Pxidxi.keys()} \psi(\alpha(xi)) Pxidxi(xi)

    **po** a function.
        returns p_o(xi), the probability a node will transmit to a 
        random neighbor given xi.

    **number_its**    int    
        number of iterations before assumed converged.
        default value is 100


    :Returns: 
    **PE**  float
        Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    psi = get_PGF(Pk)
    psiPrime = get_PGFPrime(Pk)
    
    xis = Pxidxi.keys()
    alpha = {xi: 1-po(xi) for xi in xis}
    for counter in range(number_its):
        newalpha = {}
        for xi in xis:
            newalpha[xi] = 1 - po(xi)  \
                            + po(xi)*sum(psiPrime(alpha[xihat])*Pxidxi(xihat) 
                                                for xihat in xis)/kave
        alpha = newalpha
    return 1 - sum(psi(alpha[xi])*Pxidxi[xi] for xi in xis)

def Attack_rate_discrete(Pk, p, rho = None, Sk0=None, 
                            phiS0=None, phiR0=0, number_its=100):
    r'''
    Encodes systems (6.6) and (6.10) of Kiss, Miller, & Simon.  Please 
    cite the book if using this algorithm.

    To use system (6.6), leave rho and Sk0 as None.

    :Arguments: 

    **Pk**  dict
        Pk[k] is the probability a randomly selected node has degree k.
    **tau** positive float
        per-edge transmission rate.
    **gamma** number
        per-node recovery rate
    **number_its**    int    The solution is found iteratively, so this determines 
        the number of iterations.
    **rho**  Number (default None)
        proportion of the population to be randomly infected at time 0
        only one of rho and Sk0 can be defined.
        The other (or both) should remain None.
        if rho=0, then calculates the limiting attack rate as rho->0
        (assuming an epidemic happens)
    **Sk0** dict (default None)
        only one of rho and Sk0 can be defined.  
        The other (or both) should remain None.
        Sk0 is a dict such that Sk0[k] is the probability that a 
        degree k node is susceptible at start.
    **phiS0** number (default None)
        Should only be used if Sk0 is not None.  
        If it is None, then assumes that initial introduction is 
        randomly introduced
    **phiR0** number (default 0)
        As with phiS0, only used if Sk0 is not None.

    :Returns: 
        
    **AR**  float
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoN.EoNError("at most one of rho and Sk0 can be defined")
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
    for counter in range(number_its):
        theta = 1-p + p*(phiR0 +  phiS0*psihatPrime(theta)/psihatPrime(1))
    return 1 - psihat(theta)

def Attack_rate_discrete_from_graph(G, p, initial_infecteds=None, 
                                        initial_recovereds = None, 
                                        rho = None, number_its = 100 
                                    ):
    r''' if initial_infecteds and initial_recovereds is defined, then it
    will find Sk0, phiS0, and phiR0 and then call Attack_rate_discrete.  
    
    Otherwise it calls attack_rate_discrete with rho.
    '''
    
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    
    Pk = get_Pk(G)
    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, 
                                    initial_recovereds = initial_recovereds)
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])
        SS=0
        SX=0
        for node in G.nodes():
            if status[node] == 'S':
                k = G.degree(node)
                Sk0[k] += 1./Nk[k]
                SS += sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
                SR += sum(1 for nbr in G.neighbors(node) if status[nbr]=='R')
                SX += k
        phiS0 = SS*1./SX
        phiR0 = SR*1./SX
        
    else:
        Sk0=None
        phiS0 = None
        phiR0 = 0
        
    
    return Attack_rate_discrete(Pk, p, rho = rho, Sk0=Sk0, phiS0=PhiS0, 
                                phiR0=phiR0, number_its = number_its)

def Attack_rate_cts_time(Pk, tau, gamma, number_its =100, rho = None, 
                            Sk0 = None, phiS0=None, phiR0=0):
    #tested in test_SIR_final_sizes
    '''Encodes system (6.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    
    This system predicts the fraction of nodes infected if an epidemic 
    occurs in a Configuration Model network assuming a continuous-time 
    Markovian SIR disease.  
    
    This gives the limit of the attack rate of epidemics as the initial 
    fraction infected approaches 0.

    If we look for the limit of a nonzero initial fraction infected, we 
    introduce rho or Sk0

    :Arguments: 
        
    **Pk**  dict    
            the probability a randomly selected node has degree k.

    **tau** positive float
            per-edge transmission rate.

    **gamma** number
            per-node recovery rate

    **number_its**    int
            The solution is found iteratively, so this determines the number 
            of iterations.

    **rho**  number, optional
            The initial proportion infected (defaults to None).  If None, then
            result is limit of rho->0.

    **Sk0** dict (default None)
            only one of rho and Sk0 can be defined.  
            The other (or both) should remain None.
            rho gives the fraction of nodes randomly infected.
            Sk0 is a dict such that Sk0[k] is the probability that a 
            degree k node is susceptible at start.

    :Returns: 

    **AR**  float
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoN.EoNError("at most one of rho and Sk0 can be defined")
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
        #print omega
        omega = gamma/(gamma+tau) \
                + tau*phiS0*psihatPrime(omega)/(psihatPrime(1)*(gamma+tau)) \
                + tau*phiR0/(gamma+tau)
    return 1 - psihat(omega)

def Attack_rate_cts_time_from_graph(G,  tau, gamma, initial_infecteds=None, 
                                        initial_recovereds = None, rho=None, 
                                        number_its =100):
    r'''
    Given a graph, predicts the attack rate for Configuration Model 
    networks with the given degree distribution.  This does not account 
    for any structure in G beyond degree distribution.
    
    First calculates the degree distribution and then calls 
    `Attack_rate_cts_time`.
    
    See also:
    `estimate_SIR_prob_size(G, p)` which accounts for entire structure of G, not just
    degree distribution.
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")
    
    Pk = get_Pk(G)
    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, 
                                    initial_recovereds = initial_recovereds)
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])
        SS=0
        SX=0
        for node in G.nodes():
            if status[node] == 'S':
                k = G.degree(node)
                Sk0[k] += 1./Nk[k]
                SS += sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
                SR += sum(1 for nbr in G.neighbors(node) if status[nbr]=='R')
                SX += k
        phiS0 = SS*1./SX
        phiR0 = SR*1./SX
    else:
        Sk0=None
        phiS0 = None
        phiR0 = 0

    return Attack_rate_cts_time(Pk, tau, gamma, rho=rho, Sk0=Sk0, phiS0 = phiS0,
                                phiR0=phiR0, number_its=number_its)
    
def Attack_rate_non_Markovian(Pk, Pzetadzeta, pi, number_its = 100):
    r'''
    Encodes system (6.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    :Arguments: 
    **Pk**  dict
        Pk[k] is probability of degree k.
    **Pzetadzeta** a dict.  
        gives P(zeta)dzeta for user-selected zeta.  The 
        algorithm will replace the integral with
        \sum_{zeta \in Pzetadzeta.keys()} \psi(\alpha(zeta)) Pzetadzeta(zeta)

    **pi** a function.
        returns p_i(zeta), the probability a node will receive a transmission from 
        random infected neighbor given zeta.

    **number_its**    int   (default 100)  
        number of iterations before assumed converged.
                 default value is 100
    :Returns: 
        
    **AR** float
        attack rate
            
    :Comments:
        Because of the symmetry for epidemic probability, this works by simply
        calling Epi_Prob_non_Markovian.
    '''
    return Epi_Prob_non_Markovian(Pk, Pzetadzeta, pi, number_its = number_its)


def EBCM_discrete(N, psihat, psihatPrime, p, phiS0, phiR0=0, R0=0, tmin = 0, tmax = 100,
                    return_full_data = False):
    #tested in test_basic_discrete_SIR
    r'''
    Encodes system (6.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    theta(t) = (1-p) + p(phi_R(0) 
               + phi_S(0) psihatPrime(theta(t-1))/psihatPrime(1))
    R(t) = R(t-1) + I(t-1)
    S(t) = N psihat(theta(t))
    I(t) = N-S-R

    :Arguments: 

    **N**  positive integer
        size of population
    **psihat**   function
        psihat(x) = \sum_k S(k,0) x^k
    **psihatPrime**   function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    **p**  number
        per edge transmission probability
    **phiS0**   number
        initial proportion of edges (of susceptible nodes) 
        connecting to susceptible nodes
    **phiR0**   number
        initial proportion of edges (of susceptible nodes) 
        connecting to recovered nodes
    **R0** number
        number of recovered nodes at time 0
    **tmax**  number
        maximum time
    **return_full_data**  boolean
        if False, 
           return t, S, I, R
        if True 
           return t, S, I, R, and theta

    :Returns: 
        
    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, all scipy arrays
    '''
    times = [tmin]
    theta = [1]
    R = [R0]
    S = [N*psihat(1)]
    I = [N-S[-1]-R[-1]]

    for time in range(tmin+1,tmax+1):
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
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R)
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R), scipy.array(theta)

def EBCM_discrete_from_graph(G, p, initial_infecteds=None, 
                                initial_recovereds = None, rho = None, 
                                tmin = 0, tmax=100, 
                                return_full_data=False):
    #tested in test_basic_discrete_SIR
    '''
    Takes a given graph, finds the degree distribution 
    (from which it gets psi),
    assumes a constant proportion of the population is infected at time 
    0, 
    and then uses the discrete EBCM model.
    
    :Arguments: 

    **G** Networkx Graph
        the contact network
    **p**   number
        per edge transmission probability
    **initial infecteds**  node or iterable of nodes   (default None)
        if a single node, then this node is initially infected
        if an iterable, then whole set is initially infected
        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.
        If both initial_infecteds and rho are assigned, then there
        is an error.       
    **initial_recovereds**  iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho**  float between 0 and 1   (default None)
        the fraction to be randomly infected at time 0
        If None, then rho=1/N is used where N = G.order()
    **tmax**  number
        maximum time
    **return_full_data** boolean
        if False, 
            return t, S, I, R and if True return t, S, I, R, and theta

    :Returns: 
        
    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, all scipy arrays
            
    
 '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")

    Pk = get_Pk(G)
    N= G.order()

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, 
                                    initial_recovereds = initial_recovereds)
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])
        Sk0 = 0*Nk
        SS=0
        SR=0
        SX=0
        R0=0
        for node in G.nodes():
            if status[node] == 'S':
                k = G.degree(node)
                #print(k, Sk0[k], 1./Nk[k])
                Sk0[k] += 1
                #print(k, Sk0[k], 1./Nk[k])
                SS += sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
                SR += sum(1 for nbr in G.neighbors(node) if status[nbr]=='R')
                SX += k
                #print(SS, SR, SX)
            elif status[node] == 'R':
                R0+=1
        def psihat(x):
            return sum(Pk[k]*Sk0[k]*x**k/Nk[k] for k in Pk)
        def psihatPrime(x):
            return sum(k*Pk[k]*Sk0[k]*x**(k-1)/Nk[k] for k in Pk)
        phiS0 = SS*1./SX
        phiR0 = SR*1./SX
        #print('here',Sk0, len(initial_infecteds)/sum(Nk))
        

    else:
        if rho is None:
            rho = 1./N
        def psihat(x):
            return (1-rho)*sum(Pk[k]*x**k for k in Pk)
        def psihatPrime(x):
            return (1-rho)*sum(k*Pk[k]*x**(k-1) for k in Pk)
        phiS0 = 1-rho
        phiR0 = 0
        R0 = 0
    #print(phiS0, phiR0, R0, psihat(1), psihatPrime(1))
    return EBCM_discrete(N, psihat, psihatPrime, p, phiS0, phiR0=phiR0, R0=R0, 
                        tmin = tmin, tmax = tmax, 
                        return_full_data = return_full_data)
                        
def EBCM_discrete_uniform_introduction(N, psi, psiPrime, p, rho, tmax=100, 
                                        return_full_data=False):
    #tested in test_basic_discrete_SIR
    '''
    Handles the case that the disease is introduced uniformly as opposed
    to depending on degree.

    :Arguments: 

    **N** Positive integer
        number of nodes
    **psi** function
        psi(x) = \sum P(k) x^k
    **psiPrime** function
        psiPrime(x)=d psi(x)/dx = \sum kP(k) x^{k-1}
    **p**   float (between 0 and 1)
        per edge transmission probability
    **rho**  number
        initial proportion of infected nodes
    **tmax**  number
        maximum time
    **return_full_data** boolean
        if False, 
            return t, S, I, R
        if True 
            return t, S, I, R, and theta

    :Returns: 
    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, all scipy arrays
    '''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM_discrete(N, psihat, psihatPrime, p, 1-rho, tmax=tmax, 
                            return_full_data=return_full_data)




def _dEBCM_(X, t, N, tau, gamma, psihat, psihatPrime, phiS0, phiR0):
    theta = X[0]
    R = X[1]
    
    dtheta = -tau*theta + tau*phiS0*psihatPrime(theta)/psihatPrime(1) \
                + gamma*(1-theta) + tau*phiR0

    S = N*psihat(theta)
    I = N-S-R
    dR = gamma*I
    return scipy.array([dtheta, dR])
    
def EBCM(N, psihat, psihatPrime, tau, gamma, phiS0, phiR0=0, R0=0, tmin=0, 
            tmax=100, tcount=1001, return_full_data=False):
    '''
    Encodes system (6.12) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    note : R0 is R(0), not the reproductive number


    :Arguments: 

    **N**  positive integer
        size of population
    **psihat**   function
        psihat(x) = \sum_k S(k,0) x^k
    **psihatPrime**   function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    **tau** positive float
        per edge transmission rate
    **gamma** positive float
        per node recovery rate
    **phiS0**  positive float
        initial proportion of edges (of susceptible nodes) 
        connecting to susceptible nodes
    **phiR0**   positive float
        initial proportion of edges (of susceptible nodes) 
        connecting to recovered nodes
    **R0** positive integer
        number of recovered nodes at time 0
    **tmin**  number
        start time
    **tmax**  number
        stop time
    **tcount**  positive integer
        number of distinct times to calculate
    **return_full_data**  boolean
        if False, 
            return t, S, I, R
        if True 
            return t, S, I, R, and theta

    :Returns: 
        
    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, all scipy arrays
    '''
    times = scipy.linspace(tmin, tmax, tcount)
    X0 = scipy.array([1, R0])
    X = integrate.odeint(_dEBCM_, X0, times, 
                            args = (N, tau, gamma, psihat, psihatPrime, phiS0, 
                                        phiR0))
    theta = X.T[0]
    R = X.T[1]
    S = N*psihat(theta)
    I = N-S-R
    if not return_full_data:
        return times, S, I, R
    else:
        return times, S, I, R, theta


def EBCM_from_graph(G, tau, gamma, initial_infecteds=None, 
                    initial_recovereds = None, rho = None, tmin = 0, tmax=100, 
                    tcount=1001, return_full_data=False):
    r'''
    Given network G and rho, calculates N, psihat, psihatPrime, and calls EBCM.
    '''
    if rho is not None and initial_infecteds is not None:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho is not None and initial_recovereds is not None:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")

    Pk = get_Pk(G)
    N = G.order()

    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds, 
                                    initial_recovereds = initial_recovereds)
        Nk = Counter(dict(G.degree()).values())
        maxk = max(Nk.keys())
        Nk = scipy.array([Nk[k] for k in range(maxk+1)])
        SS=0
        SR=0
        SX=0
        R0=0
        Sk0 = scipy.zeros(maxk+1)
        for node in G.nodes():
            if status[node] == 'S':
                k = G.degree(node)
                Sk0[k] += 1./Nk[k]
                SS += sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
                SR += sum(1 for nbr in G.neighbors(node) if status[nbr]=='R')
                SX += k
            elif status[node] == 'R':
                R0+=1
        def psihat(x):
            return sum(Pk[k]*Sk0[k]*x**k for k in Pk)
        def psihatPrime(x):
            return sum(k*Pk[k]*Sk0[k]*x**(k-1) for k in Pk)
        phiS0 = SS*1./SX
        phiR0 = SR*1./SX

    else:
        if rho is None:
            rho = 1./N
        def psihat(x):
            return (1-rho)*sum(Pk[k]*x**k for k in Pk)
        def psihatPrime(x):
            return (1-rho)*sum(k*Pk[k]*x**(k-1) for k in Pk)
        phiS0 = 1-rho
        phiR0 = 0
        R0 = 0
    return EBCM(N, psihat, psihatPrime, tau, gamma, phiS0, phiR0=phiR0, R0=R0, 
                        tmin = tmin, tmax = tmax, tcount=tcount,
                        return_full_data = return_full_data)
                        

def EBCM_uniform_introduction(N, psi, psiPrime, tau, gamma, rho, tmin=0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    r'''
    Handles the case that the disease is introduced uniformly as opposed
    to depending on degree.

    :Arguments: 

    **N** positive integer
        size of population
    **psi** function
        psihat(x) = \sum_k S(k,0) x^k
    **psiPrime** function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    **tau** positive float
        per edge transmission rate
    **gamma** number
        per node recovery rate
    **rho**  number
        initial proportion infected
    **tmin**  number
        start time
    **tmax**  number
        stop time
    **tcount**  integer
        number of distinct times to calculate
    **return_full_data**  boolean
        if False,
            return t, S, I, R 
        if True 
            return t, S, I, R, and theta

    :Returns: 
        
    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, all scipy arrays
'''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM(N, psihat, psihatPrime, tau, gamma, 1-rho, tmin=tmin, 
                tmax=tmax, tcount=tcount, return_full_data=return_full_data)

    
def _dEBCM_pref_mix_(X, t, rho, tau, gamma, Pk, Pnk):
    #print t
    R= X[0]
    theta = {}
    phiR={}
    for index, k in enumerate(sorted(Pk.keys())):
        theta[k] = X[1+2*index]
        phiR[k] = X[2+2*index]
    S = (1-rho)*sum([Pk[k]*theta[k]**k for k in Pk.keys()])
    #print 'S', S
    I = 1 - S - R
    returnval = [gamma*I]#xidot, I2dot, Rdot
    phiS = {}
    phiI = {}
    for k1 in Pk.keys():
        phiS[k1] = (1-rho)*sum([Pnk[k1][k2]*theta[k2]**(k2-1) for k2 in Pnk[k1].keys()])
        phiI[k1] = theta[k1] - phiS[k1]  - phiR[k1]
    for k in sorted(Pk.keys()):
        dthetak_dt = -tau*phiI[k]
        dphiRk_dt = gamma*phiI[k]
        returnval.extend([dthetak_dt, dphiRk_dt])
    return scipy.array(returnval)

#N, psihat, psihatPrime, tau, gamma, phiS0, phiR0=0, R0=0, tmin=0, 
#            tmax=100, tcount=1001, return_full_data=False
def EBCM_pref_mix(N, Pk, Pnk, tau, gamma, rho = None, tmin = 0, tmax = 100, tcount = 1001, return_full_data=False):
    r'''
    
    Encodes the system derived in exercise 6.21 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    I anticipate eventually adding an option so that the initial condition is
    not uniformly distributed.  So could give rho_k
    
    :Arguments: 
    **N**  positive integer
        number of nodes.
    **Pk**  dict  (could also be an array or a list)
        Pk[k] is the probability a random node has degree k.
    **Pnk** dict of dicts (possibly array/list)
        Pnk[k1][k2] is the probability a neighbor of a degree k1 node has
        degree k2.
    **tau** positive float
        transmission rate
    **gamma** number
        recovery rate
    **rho**  number (optional)
        initial proportion infected.  Defaults to 1/N.
    **tmin**  number (default 0)
        minimum time
    **tmax**  number (default 100)
        maximum time
    **tcount**  integer (default 1001)
        number of time points for data (including end points)
    **return_full_data**  boolean (default False)
        whether to return theta or not
            

    :Returns: 

    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**
        where theta[k] is a scipy array giving theta for degree k
            
    '''
    if rho is None:
        rho = 1./N
        
    ts = scipy.linspace(tmin, tmax, tcount)
    IC = [0] #R(0)
    for k in sorted(Pk.keys()):
        IC.extend([1,0]) #theta_k(0), phiR_k(0)
    IC = scipy.array(IC)
    X = integrate.odeint(_dEBCM_pref_mix_, IC, ts, args = (rho, tau, gamma, Pk, Pnk))
    R =  X.T[0]
    theta = {}
    phiR={}
    for index, k in enumerate(sorted(Pk.keys())):
        theta[k] = X.T[1+2*index]
        phiR[k] = X.T[1+2*index]
    
    S = (1-rho)*sum([Pk[k]*theta[k]**k for k in Pk.keys()])
    I = 1-S-R
    if return_full_data:
        return ts, N*S, N*I, N*R, theta
    else:
        return ts, N*S, N*I, N*R

def EBCM_pref_mix_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax = 100, tcount = 1001, return_full_data=False):
    r'''
    Takes a given graph, finds degree correlations, and calls EBCM_pref_mix
    
    
    I anticipate eventually adding an option so that the initial condition is
    not uniformly distributed.  So could give rho_k
    
    :Arguments: 
    **G** networkx Graph
        The contact network
    **tau** positive float
        transmission rate
    **gamma** positive float
        recovery rate
    **rho**  positive float (default None)
        initial proportion infected.  Defaults to 1/N.
    **tmin**  number (default 0)
        minimum time
    **tmax**  number (default 100)
        maximum time
    **tcount**  integer (default 1001)
        number of time points for data (including end points)
    **return_full_data**  boolean (default False)
        whether to return theta or not
            

    :Returns: 

    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, 
        where theta[k] is a scipy array giving theta for degree k
    '''

    N=G.order()
    Pk = get_Pk(G)
    Pnk = get_Pnk(G)
    return EBCM_pref_mix(N, Pk, Pnk, tau, gamma, rho = rho, tmin = tmin, tmax = tmax, tcount = tcount, return_full_data=return_full_data)
    
    
def EBCM_pref_mix_discrete(N, Pk, Pnk, p, rho = None, tmin = 0, tmax = 100, return_full_data=False):
    r'''
    
    Encodes the discrete version of exercise 6.21 of Kiss, Miller, & Simon.  
    Please cite the book if using this algorithm.

    I anticipate eventually adding an option so that the initial condition is
    not uniformly distributed.  So could give rho_k
    
    :Arguments: 
    **N**  positive integer
        number of nodes.
    **Pk**  dict  (could also be an array or a list)
        Pk[k] is the probability a random node has degree k.
    **Pnk** dict of dicts (possibly array/list)
        Pnk[k1][k2] is the probability a neighbor of a degree k1 node has
        degree k2.  
    **p**   positive float (0 <= p <= 1)
        transmission probability
    **rho**  number (optional)
        initial proportion infected.  Defaults to 1/N.
    **tmin**  number (default 0)
        minimum time
    **tmax**  number (default 100)
        maximum time
    **tcount**  integer (default 1001)
        number of time points for data (including end points)
    **return_full_data**  boolean (default False)
        whether to return theta or not
            

    :Returns: 

    if return_full_data == False:
        returns **t, S, I, R**, all scipy arrays
    if ...== True
        returns **t, S, I, R** and **theta**, 
        where theta is a dict and theta[k] is the thetas for given k.
                
    '''
    if rho is None:
        rho = 1./N

    
    times = [0]
    theta = {k:[1] for k in Pk.keys()}
    R = [0]
    S = [N*(1-rho)]
    I = [N*rho]
    phiS={k: (1-rho) for k in Pk.keys()}
    phiI={k : rho for k in Pk.keys()}
    phiR = {k:0 for k in Pk.keys()}
    for time in range(1,tmax+1):
        times.append(time)
        newtheta = {k:theta[k][-1]-p*phiI[k] for k in Pk.keys()}
        newR = R[-1]+I[-1]
        newS = N*(1-rho)*sum(Pk[k]*newtheta[k]**k for k in Pk.keys())
        newI = N-newR-newS
        for k in newtheta:
            theta[k].append(newtheta[k])
        R.append(newR)
        S.append(newS)
        I.append(newI)
        #may be worth holding on to phiS and phiI
        phiS = {k1: (1-rho)*sum([Pnk[k1][k2]*theta[k2][-1]**(k2-1) for k2 in Pnk[k1].keys()]) for k1 in Pk.keys()}
        phiR = {k:phiR[k]+(1-p)*phiI[k] for k in Pk.keys()}
        phiI = {k: theta[k][-1] - phiS[k]  - phiR[k] for k in Pk.keys()}
    if return_full_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R), theta
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R)

def EBCM_pref_mix_discrete_from_graph(G, p, rho = None, tmin = 0, tmax = 100, return_full_data=False):
    
    '''
    Takes a given graph, finds degree correlations, and calls EBCM_pref_mix_discrete
    
    :SAMPLE USE:
        
    ::

        import networkx as nx
        import EoN
        import matplotlib.pyplot as plt
        G = nx.bipartite.configuration_model([5]*300000, [2]*750000)
        t, S, I, R = EoN.basic_discrete_SIR(G, 0.6, rho = 0.002)
        tx, Sx, Ix, Rx = EoN.EBCM_pref_mix_discrete_from_graph(G, 0.6, rho=0.002, tmax=t[-1])
        plt.plot(t, I, label = 'simulation')
        plt.plot(tx, Ix, '--', label = 'analytic prediction')
        plt.legend(loc='upper right')
        plt.show()
    '''
            
    N = G.order()
    Pk = get_Pk(G)
    Pnk = get_Pnk(G)
    return EBCM_pref_mix_discrete(N, Pk, Pnk, p, rho=rho, tmin=tmin, tmax=tmax, return_full_data=return_full_data)
    
