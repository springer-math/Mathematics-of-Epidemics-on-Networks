import networkx as nx
import EoN
import scipy

def subsample(report_times, times, status1, status2=None, 
                status3 = None):
    r'''
    Given 
      S, I, and/or R as lists of numbers of nodes of the given status
      at given times

    returns them 
      subsampled at specific report_times.
    

    Arguments:

        report_times : iterable (ordered)
            times at which we want to know state of system
                   
        times : iterable (ordered)
            times at which we have the system state (assumed no change 
            between these times)
            
        statusX (X one of 1, 2 or 3) : iterable (order corresponds to times)
                          generally S, I, or R
                          number of nodes in given status.
    Returns:

        :
        If only status1 is defined
            report_status1 : scipy array gives status1 subsampled just at 
                 report_times.
                     
        If more are defined then it returns a list, either
            [report_status1, report_status2]
        or
            [report_status1, report_status2, report_status3]

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import scipy
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
        report_times = scipy.linspace(0,5,101)
        Ssum = scipy.zeros(len(report_times))
        Isum = scipy.zeros(len(report_times))
        Rsum = scipy.zeros(len(report_times))
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
        
    report_status1= scipy.array(report_status1)
    
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
    Identifies the first time at which L crosses a threshold.  
    Useful for shifting times.
    
    Arguments:
        times : list or scipy array (ordered)
            the times we have observations
        L : a list or scipy array
            order of L corresponds to times
        threshold : number
            a threshold value

    Returns:
        :
            t  (number)
                the first time at which L reaches or exceeds a threshold.

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import scipy
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
        report_times = scipy.linspace(0,5,101)
        Ssum = scipy.zeros(len(report_times))
        Isum = scipy.zeros(len(report_times))
        Rsum = scipy.zeros(len(report_times))
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



def visualize(G, plot_times, node_history, pos = None, 
                SIR = True, filetype = 'png', filenamebase = 'tmp', 
                colorS = '#009a80', colorI = '#ff2020', 
                colorR = 'gray', show_edges = True, plot_args = ()):
    r''' 
    Creates a set of plots showing statuses of nodes at different times.  By 
    default, the plot for t = 1.3 would be put into "tmp1p3.png"
    
    Arguments:
        G : NetworkX Graph
        
        plot_times : list (or array, maybe even a set)
            collection of times to output plot
        
        node_history : dict
            node_history[node] = (times, statuses)
            where times is a list of times of status change
            and statuses is the new status at that time.
             
        recovery_times : dict
            see infection_times, except this has time(s) of recovery.
            
        pos : dict, optional
            the positions to plot nodes of G.  By default spring_layout is used.
        
        SIR : boolean, default True
            True if the simulation is SIR, False if it is SIS.
            
        filetype : string (default 'png')
            the type of figure to make.
            
        filenamebase : string (default 'tmp')
            base name for output plot file names.
            
        colorS : default '#009a80'
            something that will be interpreted as a color by matplotlib.  
            Used to plot susceptible nodes.  Note, if using RGB as a tuple
            of length 3, if there are 3 susceptible nodes, the tuple will
            be interpreted as 3 different colors to use.
            
        colorI : default '#ff2020'
            see colorS
            
        colorR : default 'gray'
            see colorS
            
        show_edges : Boolean, default True
            whether the edges should be plotted
            
        plot_args : tuple, default ()
            arguments to be passed to the networkx drawing commands
            draw_networkx_nodes and draw_networkx_edges.
         
         
    :SAMPLE USE:

    ::

        import EoN
        import networkx as nx
        import matplotlib.pyplot as plt
        
        G = nx.fast_gnp_random_graph(100,0.06)
        plot_times = scipy.linspace(0,10,101) #0, 0.1, 0.2, ..., 10
        
        #let's create 101 figures for an SIS epidemic
        times, S, I, inf_times, rec_times = EoN.fast_SIS(G, 1., 1., return_full_data=True)
        EoN.visualize(G, plot_times, inf_times, rec_times, filenamebase = 'tmpSIS', SIR = False)
        
        #let's create 101 figures for an SIR epidemic
        times, S, I, R, inf_time, rec_time = EoN.fast_SIR(G, 1., 1., return_full_data=True)
        EoN.visualize(G, plot_times, inf_time, filenamebase = 'tmpSIR', rec_time)
        
    '''
    
    if pos is None:
        pos = nx.spring_layout(G)

        
    for time in plot_times:
        plt.clf()
        S = set()
        I = set()
        R = set()

        for node in G.nodes():
            if node not in node_history:
                S.add(node)
            else:
                changetimes = node_history[node][0]
                number_swaps = len([changetime for changetime in changetimes if changetime<= time])
                status = node_history[node][1][number_swaps-1]
                if status == 'S':
                    S.add(node)
                elif status == 'I':
                    I.add(node)
                else:  #won't happen in SIS case
                    R.add(node)    

        nx.draw_networkx_nodes(G, pos = pos, node_color = colorS, nodelist = list(S), *plot_args)            
        nx.draw_networkx_nodes(G, pos = pos, node_color = colorI, nodelist = list(I), *plot_args)            
        nx.draw_networkx_nodes(G, pos = pos, node_color = colorR, nodelist = list(R), *plot_args)
        if show_edges:
            nx.draw_networkx_edges(G, pos, *plot_args)
        plt.savefig(filenamebase+str(time).replace('.', 'p')+'.'+filetype)
                
