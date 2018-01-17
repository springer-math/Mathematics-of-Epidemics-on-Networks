


def visualize(G, plot_times, node_history, pos = None, 
                filetype = 'png', filenamebase = 'tmp', 
                colorS = '#009a80', colorI = '#ff2020', colorR = 'gray',
                show_edges = True, highlightSI = False,
                plot_args = (), number_by_time = False):
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
            
        highlight_SI : Boolean, default False
            whether the SI edges should be plotted in a different color.
            
        plot_args : tuple, default ()
            arguments to be passed to the networkx drawing commands
            draw_networkx_nodes and draw_networkx_edges.
         
        number_by_time : Boolean
            if True, then the filename gets the actual calculated time.  Else
            it is numbered sequentially.
         
    :SAMPLE USE:

    ::

        import EoN
        import networkx as nx
        import matplotlib.pyplot as plt
        
        G = nx.fast_gnp_random_graph(100,0.06)
        plot_times = scipy.linspace(0,10,101) #0, 0.1, 0.2, ..., 10
        
        #let's create 101 figures for an SIS epidemic
        times, S, I, node_history = EoN.fast_SIS(G, 1., 1., return_full_data=True)
        EoN.visualize(G, plot_times, node_history, filenamebase = 'tmpSIS')
        
        #let's create 101 figures for an SIR epidemic
        times, S, I, R, node_history = EoN.fast_SIR(G, 1., 1., return_full_data=True)
        EoN.visualize(G, plot_times, node_history, filenamebase = 'tmpSIR')
        
    '''
    
    if pos is None:
        pos = nx.spring_layout(G)

        
    for index, time in enumerate(plot_times):
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

        nx.draw(G, pos = pos, node_color = colorS, nodelist = list(S), edgelist = [], *plot_args)            
        nx.draw(G, pos = pos, node_color = colorI, nodelist = list(I), edgelist = [], *plot_args)            
        nx.draw(G, pos = pos, node_color = colorR, nodelist = list(R), edgelist = [], *plot_args)
        if show_edges:
            if not highlightSI:
                nx.draw(G, pos, nodelist=[], edgelist = G.edges(), *plot_args)
            else:
                SIedges = {(u,v) for u, v in G.edges() if (u in S and v in I) or (u in I and v in S)}
                nonSIedges = {edge for edge in G.edges() if edge not in SIedges}
                nx.draw(G, pos, nodelist = [], edgelist = list(SIedges), edge_color = colorI)
                nx.draw(G, pos, nodelist = [], edgelist = list(nonSIedges))
        if number_by_time:
            plt.savefig(filenamebase+str(time).replace('.', 'p')+'.'+filetype, bbox_inches='tight')
        else:
            plt.savefig(filenamebase+str(index).replace('.', 'p')+'.'+filetype, bbox_inches='tight')
                
