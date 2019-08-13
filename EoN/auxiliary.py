import networkx as nx
import EoN
import numpy as np
import random

def subsample(report_times, times, status1, status2=None, 
                status3 = None):
    r'''
    Given 
      S, I, and/or R as lists of numbers of nodes of the given status
      at given times

    returns them 
      subsampled at specific report_times.
    

    :Arguments: 

    **report_times** iterable (ordered)
        times at which we want to know state of system
                   
    **times** iterable (ordered)
        times at which we have the system state (assumed no change 
        between these times)
            
    **status1**  iterable 
        generally S, I, or R
        
        number of nodes in given status at corresponding time in times.
        
    **status2**  iterable  (optional, default None)
        generally S, I, or R
        
        number of nodes in given status at corresponding time in times.

    **status3**  iterable (optional, default None)
        generally S, I, or R
        
        number of nodes in given status at corresponding time in times.
                                
    :Returns:

    If only status1 is defined
        **report_status1** numpy array 
        gives ``status1`` subsampled just at ``report_times``.
                     
    If more are defined then it returns a list, either
        **[report_status1, report_status2]**
    or
        **[report_status1, report_status2, report_status3]**
    In each case, these are subsampled just at report_times.

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import numpy as np
        import matplotlib.pyplot as plt

        """ in this example we will run 100 stochastic simulations.
            Each simulation will produce output at a different set
            of times.  In order to calculate an average we will use
            subsample to find the epidemic sizes at a specific set
            of times given by report_times.
        """

        G = nx.fast_gnp_random_graph(10000,0.001)
        tau = 1.
        gamma = 1.
        report_times = np.linspace(0,5,101)
        Ssum = np.zeros(len(report_times))
        Isum = np.zeros(len(report_times))
        Rsum = np.zeros(len(report_times))
        iterations = 100
        for counter in range(iterations): 
            t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds = range(10))
            #t, S, I, and R have an entry for every single event.
            newS, newI, newR = EoN.subsample(report_times, t, S, I, R)
            #could also do: newI = EoN.subsample(report_times, t, I)
            plt.plot(report_times, newS, linewidth=1, alpha = 0.4)
            plt.plot(report_times, newI, linewidth=1, alpha = 0.4)
            plt.plot(report_times, newR, linewidth=1, alpha = 0.4)
            Ssum += newS
            Isum += newI
            Rsum += newR
        Save = Ssum / float(iterations)
        Iave = Isum / float(iterations)
        Rave = Rsum / float(iterations)
        plt.plot(report_times, Save, "--", linewidth = 5, label = "average")
        plt.plot(report_times, Iave, "--", linewidth = 5)
        plt.plot(report_times, Rave, "--", linewidth = 5)
        plt.legend(loc = "upper right")
        plt.savefig("tmp.pdf")

    If only one of the sample times is given then returns just that.

    If report_times goes longer than times, then this simply assumes the 
    system freezes in the final state.
    
    This uses a recursive approach if multiple arguments are defined.


    '''
    if report_times[0] < times[0]:
        raise EoN.EoNError("report_times[0]<times[0]")
        
    report_status1 = []
    next_report_index = 0
    next_observation_index = 0
    while next_report_index < len(report_times):
        while next_observation_index < len(times) and \
              times[next_observation_index]<= report_times[next_report_index]:
            candidate = status1[next_observation_index]
            next_observation_index += 1
        report_status1.append(candidate)
        next_report_index +=1
        
    report_status1= np.array(report_status1)
    
    if status2 is not None:
        if status3 is not None:
            report_status2, report_status3 = subsample(report_times, times, status2, status3)
            return report_status1, report_status2, report_status3
        else:
            report_status2 = subsample(report_times, times, status2)
            return report_status1, report_status2
    else:
        return report_status1



def get_time_shift(times, L, threshold):
    r'''
    Identifies the first time at which list/array L crosses a threshold.  
    Useful for shifting times.
    
    :Arguments: 
    **times** list or numpy array (ordered)
        the times we have observations
    **L** a list or numpy array
        order of L corresponds to times
    **threshold** number
        the threshold value

    :Returns:
        
    **t**  number
        the first time at which L reaches or exceeds threshold.

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import numpy as np
        import matplotlib.pyplot as plt

        """ in this example we will run 20 stochastic simulations.
            We plot the unshifted curves (grey) and the curves shifted 
            so that t=0 when 1% have been infected (I+R = 0.01N) (red)
        """
        plt.clf() # just clearing any previous plotting.
        
        N=100000
        kave = 10.
        G = nx.fast_gnp_random_graph(N,kave/(N-1.))
        tau = 0.2
        gamma = 1.
        report_times = np.linspace(0,5,101)
        Ssum = np.zeros(len(report_times))
        Isum = np.zeros(len(report_times))
        Rsum = np.zeros(len(report_times))
        iterations = 20
        for counter in range(iterations):
            R=[0]
            while R[-1]<1000: #if an epidemic doesn't happen, repeat
                t, S, I, R = EoN.fast_SIR(G, tau, gamma)
                print R[-1]
            plt.plot(t, I, linewidth = 1, color = 'gray', alpha=0.4)
            tshift = EoN.get_time_shift(t, I+R, 0.01*N)
            plt.plot(t-tshift, I, color = 'red', linewidth = 1, alpha = 0.4)
        plt.savefig("timeshift_demonstration.pdf")
    '''
    
    for index, t in enumerate(times):
        if L[index]>= threshold:
            break
    return t




def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, leaf_vs_root_factor = 0.5):

    '''
    Based on Joel's answer at https://stackoverflow.com/a/29597209/2966723,
    but with some modifications.  

    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    We include this because it may be useful for plotting transmission trees,
    and there is currently no networkx equivalent (though it may be coming soon).
    
    There are two basic approaches we think of to allocate the horizontal 
    location of a node.  
    
    - Top down: we allocate horizontal space to a node.  Then its ``k`` 
      descendants split up that horizontal space equally.  This tends to result
      in overlapping nodes when some have many descendants.
    - Bottom up: we allocate horizontal space to each leaf node.  A node at a 
      higher level gets the entire space allocated to its descendant leaves.
      Based on this, leaf nodes at higher levels get the same space as leaf
      nodes very deep in the tree.  
      
    We use use both of these approaches simultaneously with ``leaf_vs_root_factor`` 
    determining how much of the horizontal space is based on the bottom up 
    or top down approaches.  ``0`` gives pure bottom up, while 1 gives pure top
    down.   
    
    
    :Arguments: 
    
    **G** the graph (must be a tree)

    **root** the root node of the tree 
    - if the tree is directed and this is not given, the root will be found and used
    - if the tree is directed and this is given, then the positions will be 
      just for the descendants of this node.
    - if the tree is undirected and not given, then a random choice will be used.

    **width** horizontal space allocated for this branch - avoids overlap with other branches

    **vert_gap** gap between levels of hierarchy

    **vert_loc** vertical location of root
    
    **leaf_vs_root_factor**

    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, leftmost, width, leafdx = 0.2, vert_gap = 0.2, vert_loc = 0, 
                    xcenter = 0.5, rootpos = None, 
                    leafpos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if rootpos is None:
            rootpos = {root:(xcenter,vert_loc)}
        else:
            rootpos[root] = (xcenter, vert_loc)
        if leafpos is None:
            leafpos = {}
        children = list(G.neighbors(root))
        leaf_count = 0
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            rootdx = width/len(children)
            nextx = xcenter - width/2 - rootdx/2
            for child in children:
                nextx += rootdx
                rootpos, leafpos, newleaves = _hierarchy_pos(G,child, leftmost+leaf_count*leafdx, 
                                    width=rootdx, leafdx=leafdx,
                                    vert_gap = vert_gap, vert_loc = vert_loc-vert_gap, 
                                    xcenter=nextx, rootpos=rootpos, leafpos=leafpos, parent = root)
                leaf_count += newleaves

            leftmostchild = min((x for x,y in [leafpos[child] for child in children]))
            rightmostchild = max((x for x,y in [leafpos[child] for child in children]))
            leafpos[root] = ((leftmostchild+rightmostchild)/2, vert_loc)
        else:
            leaf_count = 1
            leafpos[root]  = (leftmost, vert_loc)
#        pos[root] = (leftmost + (leaf_count-1)*dx/2., vert_loc)
        print(leaf_count)
        return rootpos, leafpos, leaf_count

    xcenter = width/2.
    if isinstance(G, nx.DiGraph):
        leafcount = len([node for node in nx.descendants(G, root) if G.out_degree(node)==0])
    elif isinstance(G, nx.Graph):
        leafcount = len([node for node in nx.node_connected_component(G, root) if G.degree(node)==1 and node != root])
    rootpos, leafpos, leaf_count = _hierarchy_pos(G, root, 0, width, 
                                                    leafdx=width*1./leafcount, 
                                                    vert_gap=vert_gap, 
                                                    vert_loc = vert_loc, 
                                                    xcenter = xcenter)
    pos = {}
    for node in rootpos:
        pos[node] = (leaf_vs_root_factor*leafpos[node][0] + (1-leaf_vs_root_factor)*rootpos[node][0], leafpos[node][1]) 
#    pos = {node:(leaf_vs_root_factor*x1+(1-leaf_vs_root_factor)*x2, y1) for ((x1,y1), (x2,y2)) in (leafpos[node], rootpos[node]) for node in rootpos}
    xmax = max(x for x,y in pos.values())
    for node in pos:
        pos[node]= (pos[node][0]*width/xmax, pos[node][1])
    return pos