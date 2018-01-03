import networkx as nx
import EoN
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation

def _draw_edges_(G, pos, edgelist, ax, **nx_kwargs):
    drawn_edges = nx.draw_networkx_edges(G, pos, edgelist=edgelist, ax = ax, **nx_kwargs)
    return drawn_edges

def _draw_nodes_(G, pos, nodelist, colorlist, ax, **nx_kwargs):
    drawn_nodes = nx.draw_networkx_nodes(G, pos, nodelist = nodelist, node_color=colorlist, ax=ax, **nx_kwargs)
    return drawn_nodes

def _get_status_(nodelist, node_history, time):
    status = {}
    for node in nodelist:
        if node not in node_history:
                status[node] = 'S'
        else:
            changetimes = node_history[node][0]
            number_swaps = len([changetime for changetime in changetimes if changetime<= time])
            status[node] = node_history[node][1][number_swaps-1]
    return status

def _order_nodes_(nodelist, status, status_order):
    if status_order[0] is 'random':
        random.shuffle(nodelist) #random order so that there isn't an apparent bias from putting one color on top
    else: #but maybe you want to highlight, say, infected nodes.  So status_order = ['I', 'random']
        firstlist = [node for node in nodelist if status[node] == status_order[0]]
        random.shuffle(firstlist)
        if status_order[1] is 'random':
            secondlist = [node for node in nodelist if status[node] != status_order[1]]
            random.shuffle(secondlist)
            nodelist = firstlist + secondlist
        else:
            secondlist = [node for node in nodelist if status[node] == status_order[1]]
            used_states = set(status_order[0:2])
            thirdlist = [node for node in nodelist if status[node] not in used_states]
            random.shuffle(thirdlist)
            nodelist = firstlist+secondlist+thirdlist
    nodelist.reverse()

    pass

def plot_graph(G, node_history, time, pos = None, colorS = '#009a80', colorI = '#ff2020', colorR = 'gray', nodelist = None, calculated_times = None, S=None, I=None, R=None, timeseries_times = None, timelabel = '$t$', savefig = True, filename = 'tmp.png', **nx_kwargs):
    r'''
    Useful for visualizing the graph at some time point.  

    The plots come in two types. 
    
    The first option just plots the graph at a given time with the color denoting the status of the nodes.

    The second option does the above plot and can also plot S, I, and R as functions of t in the same figure as the graph plot.
    The time series plots show a vertical line marking the time observed.

    Arguments : 

        G (networkx Graph)
            The graph

        node_history (The node_history returned by EoN functions such as fast_SIS or fast_SIR when return_full_datae is True)

        time (float)
            The time at which we want a snapshot.

        pos (dict, default None)
            The positions of the nodes as used by networkx drawing functions.  If None, then uses nx.spring_layout

        colorS (color default greenish - chosen to be colorblind friendly)
            for susceptible nodes.  

        colorI (color default reddish - chosen to be colorblind friendly)
            for infected nodes.

        colorR (color default gray)
            for recovered nodes.

        nodelist (list  [if None defaults to list(G.nodes())])
            The nodes to plot.  No edges plotted that do not have both nodes in this list

        calculated_times (numpy array or list)
            The times that the EoN simulation returned.

        S (numpy array or list, default None)
            The number of susceptible nodes at each time in calculated_times.  If none of
            S, I, and R are given (and calculated times also not given), then ignores these.
            If only one, only 2 or all 3 are given (and calculated_times also) then plots 
            just those that are given

        I (numpy array or list, default None)
            see S

        R (numpy array or list, default None)
            see S

        timeseries_times (numpy array or list)
            the times to be used in the plot of the timeseries

        timelabel (string, default '$t$')
            the label to use for the horizontal axis of the timeseries.


        savefig (boolean, default True)
            Whether or not to save the figure.

        filename : (string, default 'tmp.png')
            name for saved figure.  extension tells which type to use.

        **nx_kwargs  (arguments for nx.draw commands)
            If wanting to set additional arguments such as node_size, this allows the use of any
            networkx drawing argument.  To enter an argument like this, simply do something like 
            (..., node_size = 5, ...) and this will determine the node_size the networkx drawing
            commands use.


        output : 

            If set to save a figure, saves the figure.

            Additionally returns the handle of the figure as well as for the parts which are modified 
            in the animation.
        
    '''
    timeseries_count = sum(x is not None for x in [S, I, R])

    include_timeseries = (timeseries_count>0 and calculated_times is not None) #both must be true

        
    if not include_timeseries and (timeseries_count > 0 or calculated_times is not None):
        raise EoN.EoNError("if plotting time series must include calculated_times and at least one of S, I, R if not plotting time series, do not include any of these")

            
    if pos is None:
        pos = nx.spring_layout(G)

    if include_timeseries:
        fig = plt.figure(figsize=(10,4))
        graph_ax = fig.add_subplot(121)
    else:
        fig = plt.figure()
        graph_ax = fig.add_subplot(111)
    if nodelist is None:
        nodelist = list(G.nodes())
        edgelist = list(G.edges())#[edge for edge in G.edges() if edge[0] in nodeset and edge[1] in nodeset]
    else:
        nodeset = set(nodelist)
        edgelist = [edge for edge in G.edges() if edge[0] in nodeset and edge[1] in nodeset]
        
    edge_drawing = _draw_edges_(G, pos, edgelist, graph_ax, **nx_kwargs)

    status2color = {'S': colorS, 'I':colorI, 'R':colorR}
    status = _get_status_(nodelist, node_history, time)
    colorlist = [status2color[status[node]] for node in nodelist]

    node_drawing = _draw_nodes_(G, pos, nodelist, colorlist, graph_ax, **nx_kwargs)

    cnt = 0
    if S is not None:
        cnt += 1
        S = EoN.subsample(timeseries_times, calculated_times, S)
        S_marker = _plt_series_(fig, S, cnt, colorS, timeseries_times, timeseries_count, time)
    else:
        S_marker = None
    if I is not None:
        cnt += 1
        I = EoN.subsample(timeseries_times, calculated_times, I)
        I_marker = _plt_series_(fig, I, cnt, colorI, timeseries_times, timeseries_count, time)
    else:
        I_marker = None
    if R is not None:
        cnt += 1
        R = EoN.subsample(timeseries_times, calculated_times, R)
        R_marker = _plt_series_(fig, R, cnt, colorR, timeseries_times, timeseries_count, time)
    else:
        R_marker = None
    if timeseries_count>1:
        for cnt in range(1, timeseries_count):
            ax = fig.add_subplot(timeseries_count, 2, 2*cnt)
            ax.set_xticks([])
    if include_timeseries:
        ax = fig.add_subplot(timeseries_count, 2, 2*timeseries_count)
        ax.set_xlabel(timelabel)
                                
    plt.tight_layout()
    if savefig:
        plt.savefig(basefilename + '.' + figure_type)
    return (fig, node_drawing, S_marker, I_marker, R_marker)

def _animation_update_(time, node_drawing, node_history, nodelist, status2color, saveframes, basefilename, figure_type, timeseries_data = None):
    print('tsd', timeseries_data)
    status = _get_status_(nodelist, node_history, time)
    print("nd",node_drawing)
    node_drawing.set_color([status2color[status[node]] for node in nodelist])
    if timeseries_data is not None:
        for X in timeseries_data:
            X.set_xdata([time, time])
        if saveframes:
            plt.savefig(basefilename+str(time).replace('.', 'p')+'.'+figure_type)            
        return (node_drawing, *timeseries_data)
        
    else:
        if saveframes:
            #See this answer: https://stackoverflow.com/a/41230791/2966723
            plt.savefig(basefilename+str(time).replace('.', 'p')+'.'+figure_type)            
        return node_drawing

def _plt_series_(fig, X, cnt, color, timeseries_times, timeseries_count, time0):
    r'''Plots a given time series X by adding a supblot to fig and and then
        plots it.  It puts a vertical dotted line at the first time.  It returns
        this vertical dotted line, which will be shifted in later plots.'''
    
    if X is not None:
        X_ax = fig.add_subplot(timeseries_count, 2, 2*cnt)
        Xmark_max = 1.1*max(X)
        X_ax.plot(timeseries_times, X, color = color)
        X_marker, = X_ax.plot([time0, time0], [0, Xmark_max], 'k:')
        return X_marker
    else:
        return None

def animate(G, node_history, frame_times, pos = None, colorS = '#009a80', colorI = '#ff2020', colorR = 'gray', nodelist = None,  calculated_times = None, S=None, I=None, R=None, timeseries_times = None, timelabel = '$t$', saveanimation = False, saveframes = False, basefilename = 'tmp', animation_type = 'mp4', figure_type = 'png', **nx_kwargs):
    r'''
    Animates the epidemic dynamics.

    if just the mandatory arguments are given, this will create an animation having just the graph.

    if the timeseries are given (S, I, R, etc), then the animation will have the graph and the time series.

    I want to update this to take arbitrary arguments for the networkx draw functions (in particular node size).

    Arguments :
       
        G (networkx Graph)
            the graph

        node_history (the node_history returned by sending 'return_full_data = True' to an EoN simulation)

        frame_times (list or numpy array)
            The times shown in the animation.

        pos (dict [optional])
            By default uses nx.spring_layout

        colorS (color default greenish - chosen to be colorblind friendly)
            for susceptible nodes.  

        colorI (color default reddish - chosen to be colorblind friendly)
            for infected nodes.

        colorR (color default gray)
            for recovered nodes.

        nodelist (list  [if None defaults to list(G.nodes())])
            The nodes to plot.  No edges plotted that do not have both nodes in this list

        calculated_times (numpy array or list)
            The times that the EoN simulation returned.

        S (numpy array or list)
            the timeseries of S returned by EoN simulation

        I (numpy array or list)
            the timeseries of I returned by EoN simulation

        R (numpy array or list)
            the timeseries of R returned by EoN simulation

        timeseries_times (numpy array or list)
            the times to be used in the plot of the timeseries

        timelabel (string)
            the label to use for the horizontal axis of the timeseries.

        saveanimation (boolean, default False)
            whether or not to save the animation.  Only tested on a mac, with mp4.
            likely requires installation of software to create the animation.
        
        saveframes (boolean, default False)
            whether or not to save the individual frames as images.

        basefilename (string, default 'tmp')
            What to use as the start of the file name. for frames or animation.

        animation_type (string, default 'mp4')
            What file type to save the animation.  Likely requires software installation.

        figure_type (string, default 'png')
            The file type to use if saving the frames.

        **nx_kwargs  (arguments for nx.draw commands)
            If wanting to set additional arguments such as node_size, this allows the use of any
            networkx drawing argument.  To enter an argument like this, simply do something like 
            (..., node_size = 5, ...) and this will determine the node_size the networkx drawing
            commands use.
    
    '''

    timeseries_count = sum(x is not None for x in [S, I, R])

    include_timeseries = (timeseries_count>0 and calculated_times is not None) #both must be true
    if nodelist is None:
        nodelist = list(G.nodes())
        
    if not include_timeseries and (timeseries_count > 0 or calculated_times is not None):
        raise EoN.EoNError("if plotting time series must include calculated_times and at least one of S, I, R if not plotting time series, do not include any of these")

    if include_timeseries and timeseries_times is None:
        timeseries_times = frame_times

    status2color = {'S':colorS, 'I':colorI, 'R':colorR}
    print('nx_kwargs', nx_kwargs)
    (fig, node_drawing, S_marker, I_marker, R_marker) = plot_graph(G, node_history, frame_times[0], pos=pos, colorS=colorS, colorI=colorI, colorR=colorR, nodelist=nodelist, calculated_times=calculated_times, S=S, I=I, R=R, timeseries_times=timeseries_times, timelabel = timelabel, savefig=saveframes, filename=basefilename+'.'+figure_type, **nx_kwargs)
    
    timeseries_data=[X for X in [S_marker, I_marker, R_marker] if X is not None]
    fargs = (node_drawing, node_history, nodelist, status2color, saveframes, basefilename, figure_type, timeseries_data)
    
    if saveframes: #does not repeat, because otherwise would keep rewriting as explained in this answer https://stackoverflow.com/a/41230791/2966723
        ani = FuncAnimation(fig, _animation_update_, frames = frame_times, fargs = fargs, repeat=False)
    else:
        ani = FuncAnimation(fig, _animation_update_, frames = frame_times, fargs = fargs)
        
    if saveanimation:
        ani.save(basefilename+'.' +animation_type, fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()
