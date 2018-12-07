import networkx as nx
import EoN
import matplotlib.pyplot as plt
import random
import scipy
from matplotlib.animation import FuncAnimation

from collections import defaultdict


        

class Simulation_Investigation():

    r'''Simulation_Display is a class which is used for creating a particular
    type of plot or an animation.
    
    The plot shows an image of the network at a snapshot in time.  In addition
    to the right of these plots it can show various timeseries from the simulation
    or from some other calculation.
    
    A longer term goal is to have the *_from_graph methods be directly callable and
    read in the IC and then get   
    '''
    #I want to improve how the labels, colors, linetypes etc are passed through
    #here.
    
    class _time_series_():
        def __init__(self, t, S, I, R=None, colordict=None, label=None, **kwargs):
            if colordict is None:
                raise EoN.EoNError("colordict must be defined")
            self._S_ = S
            self._I_ = I
            self._R_ = R
            self._t_ = t
            self.colordict = colordict
            self.label=label
            self.plt_kwargs = kwargs

        def _plot_(self, ax, series_string = ''):
            if self.label:
                if 'S' in series_string:
                    ax.plot(self._t_, self._S_, color = self.colordict['S'], label=self.label+': $S$', **self.plt_kwargs)
                if 'I' in series_string:
                    ax.plot(self._t_, self._I_, color = self.colordict['I'], label=self.label+': $I$', **self.plt_kwargs)
                if 'R' in series_string:
                    ax.plot(self._t_, self._R_, color = self.colordict['R'], label=self.label+': $R$', **self.plt_kwargs)
            else:
                if 'S' in series_string:
                    ax.plot(self._t_, self._S_, color = self.colordict['S'], **self.plt_kwargs)
                if 'I' in series_string:
                    ax.plot(self._t_, self._I_, color = self.colordict['I'], **self.plt_kwargs)
                if 'R' in series_string:
                    ax.plot(self._t_, self._R_, color = self.colordict['R'], **self.plt_kwargs)
        def update_kwargs(self, **kwargs):
            self.plt_kwargs.update(kwargs)
                            
 
    def __init__(self, G, node_history, transmissions, SIR = True, pos = None, 
                    colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}):
        self.G = G
        self._node_history_ = node_history
        self._transmissions_ = transmissions
        self.SIR = SIR
        self.sim_colordict = colordict
        self.pos = pos
        self.summary() #defines self._t_, self._S_, self._I_, and self._R_
        self._time_series_list_ = []
        self._simulation_time_series_ = self._time_series_(self._t_, self._S_, self._I_, self._R_, 
                                    colordict=self.sim_colordict, label = 'Simulation')
        self._time_series_list_.append(self._simulation_time_series_)
        
    def node_history(self, node):
        r'''
        returns the history of a node.
        
        :Arguments:
        **node**
            the node
        
        :Returns:
             
        **timelist, statuslist** lists
               
            the times at which the node changes status and what the new status is at each time.          
                
            '''
            
        return self._node_history_[node]
        
    def node_status(self, node, time):
        r'''
        returns the status of a given node at a given time.
    
        :Arguments:
    
        **node**
            the node
        **time** float
            the time of interest.
    
        :Returns:
    
        **status** string ('S', 'I', or 'R')
            status of node at time.
        '''
    
        changetimes = self._node_history_[node][0]
        number_swaps = len([changetime for changetime in changetimes if changetime<= time])
        status = self._node_history_[node][1][number_swaps-1]
        return status

    def get_statuses(self, nodelist=None, time=None):
        r'''
        returns the status of nodes at a given time.  
    
        :Arguments:
    
        **nodelist** iterable (default None):
            Some sort of iterable of nodes.
            If default value, then returns statuses of all nodes.
        **time** float (default None)
                the time of interest.
                if default value, then returns initial time
    
        :Returns: 
        **status** dict
            A dict whose keys are the nodes in nodelist giving their status at time.

        '''
        if nodelist is None:
            nodelist = self.G
        if time is None:
            time = self._t_[0]
        status = {}
        for node in nodelist:
            changetimes = self._node_history_[node][0]
            number_swaps = len([changetime for changetime in changetimes if changetime<= time])
            status[node] = self._node_history_[node][1][number_swaps-1]
        return status
    

    def summary(self, nodelist = None):
        r'''
        Provides the population-scale summary of the dynamics: t, S, I, and R
        
        :Arguments:
        **nodelist** (default None)
                The nodes that we want to focus on.  By default this is all nodes.
                If you want all nodes, the most efficient thing to do is to
                not include 'nodelist'.  Otherwise it will recalculate everything.
                    
        :Returns:
            
        if self.SIR is True then returns
            **t, S, I, R** --- scipy arrays 
        if self.SIR is False then returns
            **t, S, I** --- scipy arrays.
                
        Assumes that all entries in node_history start with same tmin'''
        if nodelist is None:  #calculate everything.
            nodelist =self.G
        if nodelist is self.G:
            try:
                self._t_  #after first time through, don't recalculate.
                if self.SIR:
                    return self._t_, self._S_, self._I_, self._R_
                else:
                    return self._t_, self._S_, self._I_
            except AttributeError:
                pass

        times = set()
        delta = {'S':defaultdict(int), 'I':defaultdict(int), 'R':defaultdict(int)}
        for node in nodelist:
            node_times = self._node_history_[node][0]
            node_statuses = self._node_history_[node][1]
            tmin = node_times[0] #should be the same for each node, but hard to choose a single node at start.
            times.add(tmin)
            delta[node_statuses[0]][tmin]+=1
            for new_status, old_status, time in zip(node_statuses[1:], node_statuses[:-1], node_times[1:]):
                delta[new_status][time] = delta[new_status][time]+1
                delta[old_status][time] = delta[old_status][time]-1
                times.add(time)
        t = scipy.array(sorted(list(times)))
        
        tmin = t[0]
        S=[delta['S'][tmin]]
        I=[delta['I'][tmin]]
        R=[delta['R'][tmin]]
        for time in t[1:]:
            S.append(S[-1]+delta['S'][time])
            I.append(I[-1]+delta['I'][time])
            R.append(R[-1]+delta['R'][time])

        if nodelist == self.G:   #we're going to save these to avoid recalculating 
            self._t_ = t
            self._S_ = scipy.array(S)
            self._I_ = scipy.array(I)
            if self.SIR:
                self._R_ = scipy.array(R)
            else:
                self._R_ = None
        
        else:
            S = scipy.array(S)
            I = scipy.array(I)
            if self.SIR:
                R = scipy.array(R)
            else:
                R = None
        
            if self.SIR:
                return t, S, I, R
            else:
                return t, S, I
                
    def t(self):
        r''' Returns the times of events
        Generally better to get these all through summary()'''
        return self._t_
    
    def S(self):
        r''' Returns the number susceptible at each time.
        Generally better to get these all through summary()'''
        return self._S_

    def I(self):
        r''' Returns the number infected at each time
        Generally better to get these all through summary()'''
        return self._I_

    def R(self):
        r''' Returns the number recovered at each time
        Generally better to get these all through summary()'''
        return self._R_
                
    def transmissions(self):
        r'''Returns a list of tuples (t,u,v) stating that node u infected node
        v at time t.  If v was infected at time tmin, then u is None
        
        Note - this only includes successful transmissions.  So if u tries
        to infect v, but fails because v is already infected this is not
        recorded.'''
        
        return self._transmissions_
        
    def transmission_tree(self):
        r'''
        
        Produces a MultiDigraph whose edges correspond to transmission events.  
        If SIR, then this is a tree (or a forest).
        
        :Returns: 
        
        **T** a directed Multi graph 
            T has all the information in `transmissions`.
            An edge from u to v with time t means u transmitted to v at time t.
        
        :Warning:
            
        Although we refer to this as a "tree", if the disease is SIS, there
        are likely to be cycles and/or repeated edges.  If the disease is SIR
        but there are multiple initial infections, then this will be a "forest".
        
        If it's an SIR, then this is a tree (or forest).
        
        The graph contains only those nodes that are infected at some point.
        '''
        
        T = nx.MultiDiGraph()
        
        for t, u, v in self._transmissions_:
            if u is not None:
                T.add_edge(u, v, time=t)
        return T
        
    def add_timeseries(self, t, S, I, R=None, colordict = None, label = None, **kwargs):
        r'''This allows us to include some additional timeseries for comparision
        with the simulation.  So for example, if we perform a simulation and 
        want to plot the simulation but also a prediction, this is what we 
        would use.
        
        :Arguments: 
        **t** list
            the times
        **S** list
            the number susceptible at each time
        **I** list
            the number infected at each time
        **R** list (default None)
            the number recovered at each time
        **colordict** dict  (default None)
            a dictionary mapping 'S', 'I', and (if SIR) 'R' to the color
            desired for their plots.  Defaults to the same as the simulation
        **label** (string)
            The label to be used for these plots in the legend.
        ****kwargs**
            any matplotlib key word args to affect how the curve is shown.
                
        :Returns:
        **ts** timeseries object
            
        :Modifies:
        This adds the timeseries object `ts` to the internal _time_series_list_
        
        '''
        if (R is not None and not self.SIR):
            raise EoN.EoNError("cannot define R if SIR isn't True")
        if (R is None and self.SIR):
            raise EoN.EoNError("cannot have SIR True if no R defined")
        if colordict is None:
            colordict = self.colordict
        ts = self._time_series_(t, S, I, R, self.SIR, colordict = colordict, label=label, **kwargs)
        self._time_series_list_.append(ts)
        return ts
        
    def update_ts_kwargs(self, ts, **kwargs):
        r'''Allows us to change some of the matplotlib key word arguments
        for a timeseries object
        
        :Arguments:
        **ts** (timeseries object)
                the timeseries object whose key word args we are updating.
        ****kwargs**
                the new matplotlib key word arguments
        '''
        ts.update_kwargs(**kwargs)
        
    def update_ts_label(self, ts, label):
        r'''updates the label for time series plots
        
        :Arguments:
        **ts** timeseries object
            the timeseries object whose key word args we are updating.
        **label** string
            the new label
        '''
        
        ts.label=label
    
    def update_ts_colordict(self, ts, colordict):
        r'''updates the colordict for time series plots
        
        :Arguments:
        **ts** timeseries object
            the timeseries object whose key word args we are updating.
        **colordict** dict
            the new colordict
        '''
        if not (colordict.has_key('S') and colordict.has_key('I')):
            raise EoN.EoNError("colordict must have keys 'S' and 'I'")
        if self.SIR and not colordict.has_key('R'):
            raise EoN.EoNError("if SIR, then colordict must have key 'R'")
        ts.colordict=colordict        
    
    def sim_update_kwargs(self, **kwargs):
        r'''Allows us to change some of the matplotlib key word arguments
        for the simulation.  This is identical to update_ts_kwargs except
        we don't need to tell it which time series to use.
        
        :Arguments:
        ****kwargs**
            the new matplotlib key word arguments
        '''
        self._simulation_time_series_.update_kwargs(**kwargs)

    def sim_update_label(self, label):
        r'''updates the label for the simulation in the time series plots
        
        :Arguments:
        **label** string
            the new label
        '''
        self.label=label

    def sim_update_colordict(self, colordict):
        r'''updates the colordict for the simulation 
        
        :Arguments:
        **colordict** dict
            the new colordict
        '''
        if not (colordict.has_key('S') and colordict.has_key('I')):
            raise EoN.EoNError("colordict must have keys 'S' and 'I'")
        if self.SIR and not colordict.has_key('R'):
            raise EoN.EoNError("if SIR, then colordict must have key 'R'")
        self.colordict = colordict
        
        
#    def add_timeseries_from_analytic_model(self, tau, gamma, model = EoN.EBCM_from_graph, tmin = 0, tmax = 10, tcount = 1001, SIR = True, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}, label = None):
#        r''' Uses one of the analytic models to predict the curve.
#        The analytic model needs to be one of the *_from_graph models.
#        (currently the pref mixing EBCM models do not work with this either)
#        only works for cts time models.
#        
#        Arguments:
#            tau (float)
#                transmission rate
#            gamma (float)
#                recovery rate
#            model (function)
#                A function like the *_from_graph models
#
#            
#        '''
#        node = self.G.nodes()[0]
#        tmin = self._node_history_[node][0][0]
#        initial_status = self.get_statuses(tmin)
#        
#        initial_infecteds = [node for node in self.G.nodes() if initial_status[node]=='I']
#
#        if SIR:
#            initial_recovereds = [node for node in self.G.nodes() if initial_status[node]=='R']
#            
#            t, S, I, R = model(self.G, tau, gamma, initial_infecteds=initial_infecteds, 
#                    initial_recovereds = initial_recovereds, tmin = tmin, tmax=tmax,
#                    tcount=tcount, return_full_data=False)
#            
#        else:
#            t, S, I = model(self.G, tau, gamma, initial_infecteds=initial_infecteds, 
#                    tmin = tmin, tmax=tmax,
#                    tcount=tcount, return_full_data=False)
#            R=None
#
#        self.add_timeseries(t, S=S, I=I, R=R, colordict=colordict, label=label)
#        def calculate_approximate_time_series(self, function, rho = None, **params):
#        r'''calls function to estimate time series.  If using one of the
#        networkx functions, it will use one of the X_from_graph functions'''
#        print("this function isn't up yet")
#        if rho is None:
#            if self.SIR:
#                initialS, initialI, initialR = get_IC(self.G, node_history, SIR=self.SIR)
#            else:
#                initialS, initialI = get_IC(self.G, node_history, SIR=self.SIR)
#
    def set_pos(self, pos):
        r'''Set the position of the nodes.
        
        :Arguments: 
        **pos** (dict)
            as in nx.draw_networkx
        '''
        self.pos = pos
        
    def _display_graph_(self, pos, status, nodelist, IonTop, ax, **nx_kwargs):
        #print(nodelist)
        if nodelist:
            #print('in if')
            nodeset = set(nodelist) #containment test in next line is quicker with set
            edgelist = [edge for edge in self.G.edges() if edge[0] in nodeset and edge[1] in nodeset]
        else:
            nodelist = list(self.G.nodes())
            random.shuffle(nodelist)#assume no order desired unless sent in
            edgelist = list(self.G.edges())
        #print(nodelist)
        if IonTop: #redefine nodelist order so that infected nodes on top
            I_nodes = [node for node in nodelist if status[node] == 'I']
            other_nodes = [node for node in nodelist if status[node]!='I']
            nodelist = other_nodes + I_nodes

        colorlist = [self.sim_colordict[status[node]] for node in nodelist]

        nx.draw_networkx_edges(self.G, pos, edgelist=edgelist, ax=ax, **nx_kwargs)
        drawn_nodes = nx.draw_networkx_nodes(self.G, pos, nodelist = nodelist, node_color=colorlist, ax=ax, **nx_kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        
        fakeLineS = plt.Line2D([0,0],[0,1], color=self.sim_colordict['S'], marker='o', linestyle='')
        fakeLineI = plt.Line2D([0,0],[0,1], color=self.sim_colordict['I'], marker='o', linestyle='')
        if self.SIR:
            fakeLineR = plt.Line2D([0,0],[0,1], color=self.sim_colordict['R'], marker='o', linestyle='')
            ax.legend([fakeLineS,fakeLineI,fakeLineR], ["$S$", "$I$", "$R$"])
        else:
            ax.legend([fakeLineS,fakeLineI], ["S", "I"])
            
        return drawn_nodes

    def _display_time_series_(self, fig, t, ts_plots, ts_list, timelabel):
        #the handling of the final element separately is ugly.
        #should figure out how to put it all into a single loop.
        if ts_list is None:
            ts_list = self._time_series_list_
        elif self._simulation_time_series_ not in ts_list:
            ts_list.append(self._simulation_time_series_)
        
        ts_axes = []
        time_markers = []
        ts_plot_count = len(ts_plots)        
        for cnt, ts_plot in enumerate(ts_plots[:-1]):
            if not self.SIR:
                ts_plot = "".join([x for x in ts_plot if x != 'R'])
            ax = fig.add_subplot(ts_plot_count, 2, 2*(cnt+1))
            ax.set_xticks([])
            for ts in reversed(ts_list):
                ts._plot_(ax, ts_plot)                
            ax.legend()
            ax.set_title(", ".join(ts_plot))
            tm = ax.axvline(x=t, linestyle='--', color='k')
            ts_axes.append(ax)
            time_markers.append(tm)
        ax = fig.add_subplot(ts_plot_count, 2, 2*ts_plot_count)
        ax.set_xlabel(timelabel)
        ts_plot = ts_plots[-1]
        if not self.SIR:
            ts_plot = "".join([x for x in ts_plot if x != 'R'])
        for ts in reversed(ts_list):
            ts._plot_(ax, ts_plot)                
        ax.legend()
        ax.set_title(", ".join(ts_plot))
        tm = ax.axvline(x=t, linestyle='--', color='k')
        ts_axes.append(ax)
        time_markers.append(tm)
        return ts_axes, time_markers
        
    def display(self, time, ts_plots = ['S', 'I', 'R'], ts_list = None, nodelist=None, IonTop=True, timelabel=r'$t$', pos=None, **nx_kwargs):
        r'''
        Provides a plot of the network at a specific time and (optionally) 
        some of the time series
        
        By default it plots the network and all time series.  The time series
        are plotted in 3 (for SIR) or 2 (for SIS) different plots to the right
        of the network.  There are options to control how many plots appear
        and which time series objects are plotted in it.
        
        
        :Arguments:
        **time** float
                the time for the snapshot of the network.
                
        **ts_plots** (list of strings, default ['S', 'I', 'R'])
                denotes what should appear in the timeseries plots.  The
                length of the list determines how many plots there are.  If
                entry i is 'AB' then plot i has both A and B plotted.
                .
                So the default has a plot with 'S', a plot with 'I' and another
                with 'R'.
            
        **ts_list** (list of timeseries objects - default None)
                If multiple time series have been added, we might want to plot
                only some of them.  This says which ones to plot.
                The simulation is always included.
            
        **nodelist** (list, default None)
                which nodes should be included in the network plot.  By default
                this is the entire network.  
                This also determines which nodes are on top of each other 
                (particularly if IonTop is False).
            
        **IonTop** (boolean, default True)
                In the network plot we put infected nodes on top.
            
        **timelabel** (string, default '$t$')
                the horizontal label to be used on the time series plots
                
        **pos**
                overrides self.pos for this display (but does not overwrite 
                self.pos.  Use set_pos if you want to do this)
                
        ****nx_kwargs**
                any networkx keyword arguments to go into the network plot.
            
        :Returns:
            
        **network_ax, ts_ax_list** (axis, list of axises)
            The axes for the network plot and a list of all the axes for the
            timeseries plots

        
        Notes : 
            
        If you only want to plot the graph, set ts_plots equal to [].  
         
        If you want S, I, and R on a single plot, set ts_plots equal to ['SIR']
        
        If you only want some of the timeseries objects, set ts_list to be those
        (the simulation time series will always be plotted).
        
        Examples :
            
        To show a plot where sim is the Simulation_Investigation object
        simply do
        
        ::
        
            sim.display()
            plt.show()
        
        To save it,
        
        ::
        
            sim.display()
            plt.savefig(filename).
        
        If you want to do more detailed modifications of the plots, this 
        returns the axes:
            
        ::
        
            network_ax, timeseries_axes = sim.display()
        
        '''
            
        if not self.SIR and ts_plots:
            ts_plots = [x for x in ts_plots if x != 'R']
        if ts_plots:
            fig = plt.figure(figsize=(10,4))
            graph_ax = fig.add_subplot(121)
        else:
            fig = plt.figure()
            graph_ax = fig.add_subplot(111)
        
        status = self.get_statuses(self.G, time)

        if pos is None:
            if self.pos is None:
                pos = nx.spring_layout(self.G)
            else:
                pos = self.pos

        self._display_graph_(pos, status, nodelist, IonTop, graph_ax, **nx_kwargs)
                
        if ts_plots:
            ts_ax_list, time_markers = self._display_time_series_(fig, time, ts_plots, ts_list, timelabel)
        else:
            ts_ax_list, time_markers = [], []
        plt.tight_layout()
        return graph_ax, ts_ax_list
                
    def _draw_infected_(self, pos, infected_nodes, ax, **nx_kwargs):
        drawn_infected = nx.draw_networkx_nodes(self.G, pos, nodelist = infected_nodes, node_color = self.sim_colordict['I'], **nx_kwargs)
        return drawn_infected
        
    def _update_ani_(self, time, pos, nodelist, drawn_nodes, drawn_infected, graph_ax, ts_axes, time_markers, nx_kwargs):
        status = self.get_statuses(self.G, time)
        infected_nodes = [node for node in nodelist if status[node] == 'I']
        drawn_nodes.set_color([self.sim_colordict[status[node]] for node in nodelist])
        drawn_infected[0].remove()
        drawn_infected[0] = nx.draw_networkx_nodes(self.G, pos, nodelist=infected_nodes, color = self.sim_colordict['I'], ax = graph_ax, **nx_kwargs)
        #print(len(time_markers),len(ts_axes))
        for index, ax in enumerate(ts_axes):
            time_markers[index].remove()
            time_markers[index] = ax.axvline(x=time, linestyle='--', color='k')
        return         


    def animate(self, frame_times=None, ts_plots = ['S', 'I', 'R'], 
                ts_list = None, nodelist=None, IonTop=True, timelabel=r'$t$',  
                pos = None, **nx_kwargs):
        r'''As in display, but this produces an animation.  
        
        To display an animation where sim is the Simulation_Investigation object
        simply do
        
        sim.animate()
        plt.show()
        
        To save an animation [on a mac with appropriate additional libraries
        installed], you can do
        
        ani = sim.animate()
        ani.save(filename, fps=5, extra_args=['-vcodec', 'libx264'])
        
        here ani is a matplotlib animation.
        See https://matplotlib.org/api/_as_gen/matplotlib.animation.Animation.save.html
        for more about the save command for matplotlib animations.
        
        :Arguments:
        The same as in display, except that time is replaced by frame_times
            
        **frame_times** (list/scipy array)
            The times for animation frames.  If nothing is given, then it
            uses 101 times between 0 and t[-1]
                
        **ts_plots** (list of strings, default 'S', 'I', 'R')
            The default means that there will be 3 plots showing time series
            with the first showing S, the second I, and the third R.
            
            If one of these is not wanted, it can simply not be included in
            the list. 
                
            Alternately if we want more than one to appear on the same plot
            the entry should be something like 'SI' or 'IR' or 'SIR' and the
            time series plot will show all of the plots.
                
            If this is an empty list, then only the network is shown, but 
            with a larger figure.

            
        **ts_list** list of timeseries objects  (default None)
            If multiple time series have been added, we might want to plot
            only some of them.  This says which ones to plot.
            The simulation is always included.
            
        **nodelist** list (default None)
            which nodes should be included in the network plot.  By default
            this is the entire network.  
            This also determines which nodes are on top of each other 
            (particularly if IonTop is False).
            
        **IonTop** boolean (default True)
            In the network plot we put infected nodes on top.
            
        **timelabel** string (default '$t$')
            the horizontal label to be used on the time series plots
            
        **pos** dict   (default None)
            overrides self.pos for this display (but does not overwrite 
            self.pos.  Use set_pos if you want to do this)
                
        ****nx_kwargs**
            any networkx keyword arguments to go into the network plot.
                
            
        '''
        
        if frame_times is None:
            frame_times = scipy.linspace(0,self._t_[-1], 101)
        if not self.SIR and ts_plots:
            ts_plots = [x for x in ts_plots if x != 'R']
        if ts_plots:
            fig = plt.figure(figsize=(10,4))
            graph_ax = fig.add_subplot(121)
        else:
            fig = plt.figure()
            graph_ax = fig.add_subplot(111)
        
        initial_status = self.get_statuses(self.G, frame_times[0])
    
        if pos is None:
            if self.pos is None:
                pos = nx.spring_layout(self.G)
            else:
                pos = self.pos
    
        if nodelist is None:
            nodelist = list(self.G.nodes())
            random.shuffle(nodelist)
        drawn_nodes = self._display_graph_(pos, initial_status, nodelist, False, graph_ax, **nx_kwargs)
        infected_nodes = [node for node in self.G if initial_status[node]=='I']
        drawn_infected = [self._draw_infected_(pos, infected_nodes, graph_ax, **nx_kwargs)] #making it a list so that I can change the entry in the list while still passing the same object
        
        
        if ts_plots:
            ts_axes, time_markers = self._display_time_series_(fig, frame_times[0], ts_plots, ts_list, timelabel)
        else:
            ts_axes, time_markers = [], []
        plt.tight_layout()
        
        fargs = (pos, nodelist, drawn_nodes, drawn_infected, graph_ax, ts_axes, time_markers, nx_kwargs)

        ani = FuncAnimation(fig, self._update_ani_, frames = frame_times, fargs = fargs, repeat=False)

        return ani
                
        #to show, do
        #simulation.animation()
        #plt.show()
        #to save do
        #ani = simulation.animation()
        #ani.save(filename, fps=5, extra_args=['-vcodec', 'libx264'])

