import networkx as nx
import EoN
import scipy
import matplotlib.pyplot as plt

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
                   
    **times** : iterable (ordered)
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
        **report_status1** scipy array 
        gives status1 subsampled just at report_times.
                     
    If more are defined then it returns a list, either
        **[report_status1, report_status2]**
    or
        **[report_status1, report_status2, report_status3]**
    In each case, these are subsampled just at report_times.

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
    Identifies the first time at which list/array L crosses a threshold.  
    Useful for shifting times.
    
    :Arguments: 
    **times** list or scipy array (ordered)
        the times we have observations
    **L** a list or scipy array
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


