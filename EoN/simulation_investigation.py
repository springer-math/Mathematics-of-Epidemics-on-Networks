import networkx as nx
import EoN
import matplotlib.pyplot as plt
import random
import numpy as np
from matplotlib.animation import FuncAnimation

from collections import defaultdict
        

class Simulation_Investigation():

    r'''Simulation_Display is a class which is used for creating a particular
    type of plot or an animation.
    
    The plot shows an image of the network at a snapshot in time.  In addition
    to the right of these plots it can show various timeseries from the simulation
    or from some other calculation.
    
    A longer term goal is to have the *_from_graph methods be directly callable and
    read in the IC and then get the appropriate time series.
    '''
    #I want to improve how the labels, colors, linetypes etc are passed through
    #here.
    
    
    #within the Simulation_Investigation class there sits a _time_series_ class
    
    class _time_series_():
        
        
        def __init__(self, ts_data, color_dict, label=None, tex = True,
                    **kwargs):
            r'''
            
            :Arguments:
                
            **ts_data** a pair (t, D)
                where ``t`` is a numpy array of times and ``D`` is a dict such that
                ``D[status]`` is a numpy array giving the number of individuals
                of given status at corresponding time in ``t``.
                
            **color_dict** a dict
                ``color_dict[status]`` is the color to be used for status in plots
                
            **label** a string
                The label to be used for these plots
                
            **tex** A boolean
                tells whether the status should be rendered as tex math mode
                or not in the labels.
                
            **kwargs** key word arguments to be passed along to the plotting command
               
                '''
            
            self._t_ = ts_data[0]
            self._D_ = ts_data[1]
            self._tex_ = tex
            self.color_dict = color_dict
            self.label=label
            self.plt_kwargs = kwargs

        def _plot_(self, ax, statuses_to_plot):
            if self.label:
                if self._tex_:
                    for status in statuses_to_plot:
                        if status in self._D_:
                            ax.plot(self._t_, self._D_[status], color = self.color_dict[status], label=self.label+': ${}$'.format(status), **self.plt_kwargs)
                else:
                    for status in statuses_to_plot:
                        if status in self._D_:
                            ax.plot(self._t_, self._D_[status], color = self.color_dict[status], label=self.label+': {}'.format(status), **self.plt_kwargs)
            else:
                for status in statuses_to_plot:
                    if status in self._D_:
                        ax.plot(self._t_, self._D_[status], color = self.color_dict[status], **self.plt_kwargs)

        def update_kwargs(self, **kwargs):
            self.plt_kwargs.update(kwargs)
                            
 
    def __init__(self, G, node_history, transmissions = None, 
                    possible_statuses = None, pos = None, color_dict=None,
                    tex = True):
                    
        r'''
                                
        :Arguments:
            
        **G** The graph
        **node_history** (dict)
            ``node_history[node]`` is a tuple (times, statuses) where 
            - ``times`` is a list of the times at which ``node`` changes status.  
               the first entry is the initial time.
            - ``statuses`` is a list giving the status the new status of the node
               at each corresponding time.
        **transmissions** (list)
            Each event which is induced by a neighbor appears (in order) in 
            ``transmissions``.  It appears as a triple ``(time, source, target)``
            where 
            - ``time`` is the time of the event
            - ``source`` is the neighbor inducing the transition.
            - ``target`` is the node undergoing the transition.
        **possible_statuses** list (default None)
            a list of the statuses to be considered.
            If not given, then defaults to the values in node_history.
        **pos** (dict - default None)
            The ``pos`` argument to be given to the networkx plotting commands.
        **color_dict** (dict - default None)
            A dictionary stating for each status what color is to be used when
            plotting.
            If not given and your statuses are ``'S'``, and ``'I'`` or they are 
            ``'S'``, ``'I'``, and ``'R'``, it will attempt to use a greenish color for
            ``'S'`` and a reddish color for ``'I'`` and gray for ``'R'``.  These
            should be color-blind friendly, despite appearing green/red to
            me. 
            Otherwise if not given, it will cycle through a set of 7 colors
            which I believe are color-blind friendly.  If you have more than
            7 statuses, you probably want to set your own color_dict.
            
        **tex** Boolean (default `True`)
            If 'True`, then labels for statuses will be in tex's math mode
            If ``False``, just plain text.
        '''
        
        if possible_statuses is None:
            ps = set()
            for node in node_history:
                ps = ps.union(set(node_history[node]))
            possible_statuses = list(ps)
            
        if color_dict is None:
            if set(possible_statuses) == set(['S', 'I', 'R']):
                color_dict = {'S':'#009a80','I':'#ff2000', 'R':'gray'}
            elif set(possible_statuses) == set(['S', 'I']):
                color_dict = {'S':'#009a80','I':'#ff2000'}
            else:
                colors = ['#FF2000', '#009A80', '#5AB3E6', '#E69A00', '#CD9AB3', '#0073B3','#F0E442']
                color_dict = {status:colors[index%len(colors)] for index, status in enumerate(possible_statuses)}
                
                
        self.G = G
        self._node_history_ = node_history
        self._transmissions_ = transmissions
        self._tex_ = tex
        self._possible_statuses_ = possible_statuses
        self.sim_color_dict = color_dict
        self.pos = pos #don't go through the effort to define this until a plot
                       #is made
        self.summary() #defines self._t_, self._D_
        self._time_series_list_ = []
        self._simulation_time_series_ = self._time_series_(self._summary_,
                                    color_dict=self.sim_color_dict, 
                                    label = 'Simulation', tex = tex
                                    )
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
    
        **status** string (such as `'S'`, `'I'`, or `'R'`)
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
        
        
        
        Provides the population-scale summary of the dynamics.  It returns
        a numpy array t as well as numpy arrays for each of the ``possible_statuses``
        giving how many nodes had that status at the corresponding time.
        
        
        Assumes that all entries in node_history start with same tmin

        :Arguments:
        **nodelist** (default None)
                The nodes that we want to focus on.  By default this is all nodes.
                If you want all nodes, the most efficient thing to do is to
                not include ``'nodelist'``.  Otherwise it will recalculate everything.
                    
        :Returns:
           
        **summary** tuple
            a pair (t, D) where 
            - t is a numpy array of times and
            - D is a dict whose keys are the possible statuses and whose values
                are numpy arrays giving the count of each status at the specific
                times.
            If nodelist is empty, this is for the entire graph.  Otherwise
            it is just for the node in nodelist.
        '''
        if nodelist is None:  #calculate everything.
            nodelist =self.G
        if nodelist is self.G:
            try:
                self._summary_  #after first time through, don't recalculate.
                return self._summary_
                
            except AttributeError: #hey, it's the first time through, let's calculate
                pass

        times = set()
        delta = {status:defaultdict(int) for status in self._possible_statuses_}
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
        t = np.array(sorted(list(times)))
        tmin = t[0]
        
        mysummary = (t, {status:[delta[status][tmin]] for status in self._possible_statuses_})
        for time in t[1:]:
            for status in self._possible_statuses_:
                mysummary[1][status].append(mysummary[1][status][-1]+delta[status][time])
        
        for status in self._possible_statuses_:
            mysummary[1][status] = np.array(mysummary[1][status])
            
        if nodelist == self.G:   #
            self._summary_ = mysummary
            self._t_ = t
            self._D_ = mysummary[1]

        return mysummary

    
    def t(self):
        r''' Returns the times of events
        Generally better to get these all through summary()'''
        return self._summary_[0]
    
    def S(self):
        r''' 
        
        If ``'S'`` is a state, then this will return the number susceptible at each time. 
        
        Else it raises an error

        Generally better to get these all through ``summary()``  '''

        if 'S' in self._possible_statuses_:
            return self._summary_[1]['S']
        else:
            raise EoN.EoNError("'S' is not a possible status")
            
    def I(self):
        r''' 
        See notes for S
        
        Returns the number infected at each time
        Generally better to get these all through summary()'''
        if 'I' in self._possible_statuses_:
            return self._summary_[1]['I']
        else:
            raise EoN.EoNError("'I' is not a possible status")

    def R(self):
        r''' 
        See notes for S
        
        Returns the number recovered at each time
        Generally better to get these all through summary()'''
        if 'R' in self._possible_statuses_:
            return self._summary_[1]['R']
        else:
            raise EoN.EoNError("'R' is not a possible status")
                
    def transmissions(self):
        r'''Returns a list of tuples (t,u,v) stating that node u infected node
        v at time t.  In the standard code, if v was already infected at tmin, then 
        the source is None
        
        Note - this only includes successful transmissions.  So if u tries
        to infect v, but fails because v is already infected this is not
        recorded.'''
        
        if self._transmissions_ is None:
            raise EoN.EoNError("transmissions were not provided when created")
        return self._transmissions_
        
    def transmission_tree(self):
        r'''
        
        Produces a MultiDigraph whose edges correspond to transmission events.  
        If SIR, then this is a tree (or a forest).
        
        :Returns: 
        
        **T** a directed Multi graph 
            T has all the information in ``transmissions``.
            An edge from u to v with time t means u transmitted to v at time t.
        
        :Warning:
            
        Although we refer to this as a "tree", if the disease is SIS, there
        are likely to be cycles and/or repeated edges.  If the disease is SIR
        but there are multiple initial infections, then this will be a "forest".
        
        If it's an SIR, then this is a tree (or forest).
        
        The graph contains only those nodes that are infected at some point.
        
        '''
        
        if self._transmissions_ is None:
            raise EoN.EoNError("transmissions were not provided when created")

        T = nx.MultiDiGraph()
        
        for t, u, v in self._transmissions_:
            if u is not None:
                T.add_edge(u, v, time=t)
        return T
        
    def add_timeseries(self, ts_data, color_dict = None, label = None, tex = None,
                        **kwargs):
        r'''
        
        This allows us to include some additional timeseries for comparision
        with the simulation.  So for example, if we perform a simulation and 
        want to plot the simulation but also a prediction, this is what we 
        would use.
        
        :Arguments: 
        **ts_data** a pair (t, D)
            where t is a numpy array of times
            and D is a dict
            where D[status] is the number of individuals of given status at 
            corresponding time.
        **color_dict** dict  (default None)
            a dictionary mapping statuses to the color
            desired for their plots.  Defaults to the same as the simulation
        **label** (string)
            The label to be used for these plots in the legend.
        **tex** (boolean)
            Tells whether status should be rendered in tex's math mode in 
            labels.  Defaults to whatever was done for creation of this
            simulation_investigation object.
        ****kwargs**
            any matplotlib key word args to affect how the curve is shown.
                
        :Returns:
        **ts** timeseries object
            
        :Modifies:
        This adds the timeseries object ``ts`` to the internal ``_time_series_list_``
        
        '''
        if color_dict is None:
            color_dict = self.color_dict
        ts = self._time_series_(ts_data, color_dict = color_dict, label=label, 
                                tex=self._tex_, **kwargs)
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
        
    def update_ts_tex(self, ts, tex):
        r'''updates the tex flag for time series plots
        
        :Arguments:
        **ts** (timeseries object)
                the timeseries object whose key word args we are updating.
        **tex**
                the new value for ``tex``
        '''
        ts._tex_=tex
        
    def update_ts_label(self, ts, label):
        r'''updates the label for time series plots
        
        :Arguments:
        **ts** timeseries object
            the timeseries object whose key word args we are updating.
        **label** string
            the new label
        '''
        
        ts.label=label
    
    def update_ts_color_dict(self, ts, color_dict):
        r'''
        
        updates the color_dict for time series plots
        
        :Arguments:
        **ts** timeseries object
            the timeseries object whose key word args we are updating.
        **color_dict** dict
            the new color_dict
        '''

        for status in ts._D_:
            if status not in color_dict:
                raise EoN.EoNError("Status {} is not in color_dict".format(status))
        ts.color_dict=color_dict        
    
    def sim_update_kwargs(self, **kwargs):
        r'''Allows us to change some of the matplotlib key word arguments
        for the simulation.  This is identical to update_ts_kwargs except
        we don't need to tell it which time series to use.
        
        :Arguments:
        ****kwargs**
            the new matplotlib key word arguments
        '''
        self._simulation_time_series_.update_kwargs(**kwargs)
    
    def sim_update_tex(self, tex):
        r'''updates the tex flag for the simulation in the time series plots
        and in the network plots
        
        :Arguments:
        **tex** string
            the new value of ``tex``
        '''
        self._tex_=tex
        self._simulation_time_series_._tex_=tex
        
    def sim_update_label(self, label):
        r'''updates the label for the simulation in the time series plots
        
        :Arguments:
        **label** string
            the new ``label``
        '''
        self.label=label

    def sim_update_color_dict(self, color_dict):
        r'''
        
        updates the color_dict for the simulation 
        
        :Arguments:
        **color_dict** dict
            the new color_dict
        '''
        for status in self._possible_statuses_:
            if status not in color_dict:
                raise EoN.EoNError("Status {} is not in color_dict".format(status))

        self.sim_color_dict = color_dict
        self._simulation_time_series_.color_dict=color_dict
        
        
#    def add_timeseries_from_analytic_model(self, tau, gamma, model = EoN.EBCM_from_graph, tmin = 0, tmax = 10, tcount = 1001, SIR = True, color_dict={'S':'#009a80','I':'#ff2000', 'R':'gray'}, label = None):
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
#        self.add_timeseries(t, S=S, I=I, R=R, color_dict=color_dict, label=label)
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
            as in ``nx.draw_networkx``
        '''
        self.pos = pos
        
    def _display_graph_(self, pos, nodestatus, nodelist, status_order, statuses_to_plot, ax, **nx_kwargs):

        '''
        
        :Arguments:
            
        **pos** (dict)
            position as for networkx
        
        **nodestatus** (dict)
            status of all nodes at given time
            
        **nodelist** (list)
            a list of the nodes to plot.  This partially determines which
            nodes appear on top
            
        **status_order**  list of statuses  
            Each status will appear on top of all later statuses.  If list 
            empty or ``False``, will ignore.
            Any statuses not appearing in list will simply be below those on the
            list and will not have priority by status.
            
        **statuses_to_plot** list of statuses to plot.
            If given, then the other nodes will be left invisible when plotting
            but I think this requires networkx v2.3 or later.
            
            
        **ax** axis
        
        **nx_kwargs**
            
            
        '''

        
        if nodelist:
            nodeset = set(nodelist) #containment test in next line is quicker with set
            edgelist = [edge for edge in self.G.edges() if edge[0] in nodeset and edge[1] in nodeset]
        else:
            nodelist = list(self.G.nodes())
            random.shuffle(nodelist)#assume no order desired unless sent in
            edgelist = list(self.G.edges())
        if status_order: #redefine nodelist order so that particular status on top
            nodes_by_status = [[node for node in nodelist if nodestatus[node] == status] for status in reversed(status_order)]
            # I_nodes = [node for node in nodelist if nodestatus[node] == 'I']
            other_nodes = [node for node in nodelist if nodestatus[node] not in status_order]
            nodelist = other_nodes
            for L in nodes_by_status:
                nodelist.extend(L)

        color_list = [self.sim_color_dict[nodestatus[node]] if nodestatus[node] in statuses_to_plot else "None" for node in nodelist]

        nx.draw_networkx_edges(self.G, pos, edgelist=edgelist, ax=ax, **nx_kwargs)
        drawn_nodes = nx.draw_networkx_nodes(self.G, pos, nodelist = nodelist, node_color=color_list, ax=ax, **nx_kwargs)
        if "with_labels" in nx_kwargs and nx_kwargs['with_labels']==True:
            nx.draw_networkx_labels(self.G, pos)
        if "with_edge_labels" in nx_kwargs and nx_kwargs['with_edge_labels'] == True:
            nx.draw_networkx_edge_labels(self.G, pos)
        ax.set_xticks([])
        ax.set_yticks([])
        
        fakelines = []
        for status in statuses_to_plot:
            fakelines.append(plt.Line2D([0,0],[0,1], color=self.sim_color_dict[status], marker = 'o', linestyle = ''))
        
        if self._tex_:
            ax.legend(fakelines, ['${}$'.format(status) for status in statuses_to_plot])
        else:
            ax.legend(fakelines, statuses_to_plot)

        return drawn_nodes

    def _display_time_series_(self, fig, t, ts_plots, ts_list, timelabel):
        
        '''        
        :ARGUMENTS:
            
        **fig** a matplotlib figure
        
        **t** float
                the time for the snapshot of the network.
                
        **ts_plots** (list of lists or list of strings)
                lists such as ``[['S'], ['I'], ['R']]``  or ``[['S', 'I'], ['R']]``
                
                equivalently ``['S', 'I', 'R']`` and ``['SI', 'R']`` will do the same
                but is problematic if a status has a string longer than 1.
                
                denotes what should appear in the timeseries plots.  The
                length of the list determines how many plots there are.  If
                entry i is ``['A', 'B']`` then plot i has both ``'A'`` and ``'B'`` plotted.
                .
                So ``[['S'], ['I'], ['R']]``  or ``['SIR']`` will result in 
                3 plots, one with just ``'S'``, one with just ``'I'`` and one with just ``'R'``
                
                while ``[['S', 'I'], ['R']]`` or ``['SI', 'R']`` will result in 
                2 plots, one with both ``'S'`` and ``'I'`` and one with just ``'R'``.

        
        **ts_list** (list of timeseries objects - default ``None``)
                If multiple time series have been added, we might want to plot
                only some of them.  This says which ones to plot.
                The simulation is always included.
        
        **timelabel** (string, default ``'$t$'``)
                the horizontal label to be used on the time series plots
        '''
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
            ax = fig.add_subplot(ts_plot_count, 2, 2*(cnt+1))
            ax.set_xticks([])
            
            for ts in reversed(ts_list):
                ts._plot_(ax, ts_plot)  
                              
            ax.legend()
            
            if self._tex_:
                ax.set_title(", ".join(['${}$'.format(status) for status in ts_plot]))
            else:
                ax.set_title(", ".join(ts_plot))
                
            tm = ax.axvline(x=t, linestyle='--', color='k')
            ts_axes.append(ax)
            time_markers.append(tm)
            
        ax = fig.add_subplot(ts_plot_count, 2, 2*ts_plot_count)
        ax.set_xlabel(timelabel)
        ts_plot = ts_plots[-1]
        
        for ts in reversed(ts_list):
            ts._plot_(ax, ts_plot)                
        ax.legend()
        
        if self._tex_:
            ax.set_title(", ".join(['${}$'.format(status) for status in ts_plot]))
        else:
            ax.set_title(", ".join(ts_plot))

        tm = ax.axvline(x=t, linestyle='--', color='k')
        ts_axes.append(ax)
        time_markers.append(tm)
        return ts_axes, time_markers
        
    def display(self, time, ts_plots = None, ts_list = None, nodelist=None, 
                status_order=False, timelabel=r'$t$', pos=None, statuses_to_plot = None,
                **nx_kwargs):
        r'''
            
        Provides a plot of the network at a specific time and (optionally) 
        some of the time series
        
        By default it plots the network and all time series.  The time series
        are plotted in 3 (for SIR) or 2 (for SIS) different plots to the right
        of the network.  There are options to control how many plots appear
        and which time series objects are plotted in it.
        
        We can make the number of time series plots to the right be zero
        by setting ts_plots to be an empty list.
        
        
        :Arguments:
        **time** float
                the time for the snapshot of the network.
                
        **ts_plots** (list of strings, defaults to ``statuses_to_plot``, which defaults 
                      to ``self._possible_statuses_``)
        
                if ``[]`` or ``False`` then the display only shows the network.  
                
                lists such as ``[['S'], ['I'], ['R']]``  or ``[['S', 'I'], ['R']]``
                
                equivalently ``['S', 'I', 'R']`` and ``['SI', 'R']`` will do the same
                but is problematic if a status has a string longer than 1.
                
                denotes what should appear in the timeseries plots.  The
                length of the list determines how many plots there are.  If
                entry i is ``['A', 'B']`` then plot i has both ``'A'`` and ``'B'`` plotted.
                .
                So ``[['S'], ['I'], ['R']]``  or ``['SIR']`` will result in 
                3 plots, one with just ``'S'``, one with just ``'I'`` and one with just ``'R'``
                
                while ``[['S', 'I'], ['R']]`` or ``['SI', 'R']`` will result in 
                2 plots, one with both ``'S'`` and ``'I'`` and one with just ``'R'``.

                Defaults to the possible_statuses
                            
        **ts_list** (list of timeseries objects - default None)
                If multiple time series have been added, we might want to plot
                only some of them.  This says which ones to plot.
                The simulation is always included.
            
        **nodelist** (list, default None)
                which nodes should be included in the network plot.  By default
                this is the entire network.  
                This also determines which nodes are on top of each other 
                (particularly if ``status_order`` is ``False``).
            
        **status_order**  list of statuses  default ``False``
            Each status will appear on top of all later statuses.  If list 
            empty or ``False``, will ignore.
            Any statuses not appearing in list will simply be below those on the
            list and will not have priority by status.
            
        **timelabel** (string, default ``'$t$'``)
                the horizontal label to be used on the time series plots
                
        **pos**
                overrides self.pos for this display (but does not overwrite 
                self.pos.  Use set_pos if you want to do this)
                
        **statuses_to_plot** list of statuses to plot.
            If given, then the other nodes will be left invisible when plotting
            but I think this requires networkx v2.3 or later.
                
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
                      
        if statuses_to_plot is None:
            statuses_to_plot = self._possible_statuses_  
        if ts_plots is None:
            ts_plots = [[x] for x in statuses_to_plot]
        if ts_plots:
            fig = plt.figure(figsize=(10,4))
            graph_ax = fig.add_subplot(121)
        else:
            fig = plt.figure()
            graph_ax = fig.add_subplot(111)
        
        nodestatus = self.get_statuses(self.G, time)

        if pos is None:
            if self.pos is None:
                pos = nx.spring_layout(self.G)
            else:
                pos = self.pos
                
        self._display_graph_(pos, nodestatus, nodelist, status_order, statuses_to_plot, graph_ax, **nx_kwargs)
                
        if ts_plots:
            ts_ax_list, time_markers = self._display_time_series_(fig, time, ts_plots, ts_list, timelabel)
        else:
            ts_ax_list, time_markers = [], []
        plt.tight_layout()
        return graph_ax, ts_ax_list
       
    def _draw_specific_status(self, pos, nodes, status, ax, **nx_kwargs):
        drawn = nx.draw_networkx_nodes(self.G, pos, nodelist = nodes, node_color = self.sim_color_dict[status], **nx_kwargs) 
        return drawn
        
    def _update_ani_(self, time, pos, nodelist, drawn_nodes, drawn_elevated, status_order, graph_ax, ts_axes, time_markers, nx_kwargs):
        '''
        '''
        nodestatus = self.get_statuses(self.G, time)
        
        drawn_nodes.set_color([self.sim_color_dict[nodestatus[node]] for node in nodelist])
        for status in reversed(status_order):
            nodes_with_status = [node for node in nodelist if nodestatus[node] == status]
            drawn_elevated[status][0].remove()
            drawn_elevated[status][0] =  nx.draw_networkx_nodes(self.G, pos, nodelist=nodes_with_status, color = self.sim_color_dict[status], ax = graph_ax, **nx_kwargs)
        for index, ax in enumerate(ts_axes):
            time_markers[index].remove()
            time_markers[index] = ax.axvline(x=time, linestyle='--', color='k')
        return         


    def animate(self, frame_times=None, ts_plots = None, 
                ts_list = None, nodelist=None, status_order=False, timelabel=r'$t$',  
                pos = None, statuses_to_plot = None, **nx_kwargs):
        r'''
        
        As in display, but this produces an animation.  
        
        To display an animation where sim is the Simulation_Investigation object
        simply do
        
        ::
        
            sim.animate()
            plt.show()
        
        To save an animation [on a mac with appropriate additional libraries
        installed], you can do
        
        ::
        
            ani = sim.animate()
            ani.save(filename, fps=5, extra_args=['-vcodec', 'libx264'])
        
        here ``ani`` is a matplotlib animation.
        See 
        
        https://matplotlib.org/api/_as_gen/matplotlib.animation.Animation.save.html
        
        for more about the save command for matplotlib animations.
        
        :Arguments:
        The same as in display, except that time is replaced by frame_times
            
        **frame_times** (list/numpy array)
            The times for animation frames.  If nothing is given, then it
            uses 101 times between 0 and t[-1]
                
        **ts_plots** (list of strings, defaults to ``statuses_to_plot``, which defaults 
                      to ``self._possible_statuses_``)
        
                if ``[]`` or ``False`` then the display only shows the network.  
                
                lists such as ``[['S'], ['I'], ['R']]``  or ``[['S', 'I'], ['R']]``
                
                equivalently ``['S', 'I', 'R']`` and ``['SI', 'R']`` will do the same
                but is problematic if a status has a string longer than 1.
                
                denotes what should appear in the timeseries plots.  The
                length of the list determines how many plots there are.  If
                entry i is ``['A', 'B']`` then plot i has both ``'A'`` and ``'B'`` plotted.
                .
                So ``[['S'], ['I'], ['R']]``  or ``['SIR']`` will result in 
                3 plots, one with just ``'S'``, one with just ``'I'`` and one with just ``'R'``
                
                while ``[['S', 'I'], ['R']]`` or ``['SI', 'R']`` will result in 
                2 plots, one with both ``'S'`` and ``'I'`` and one with just ``'R'``.

                Defaults to the possible_statuses
            
        **ts_list** list of timeseries objects  (default None)
            If multiple time series have been added, we might want to plot
            only some of them.  This says which ones to plot.
            The simulation is always included.
            
        **nodelist** list (default None)
            which nodes should be included in the network plot.  By default
            this is the entire network.  
            This also determines which nodes are on top of each other 
            (particularly if status_order is ``False``).
            
        **status_order**  list of statuses  default ``False``
            Each status will appear on top of all later statuses.  If list 
            empty or ``False``, will ignore.
            Any statuses not appearing in list will simply be below those on the
            list and will not have priority by status.
            
        **timelabel** string (default '$t$')
            the horizontal label to be used on the time series plots
            
        **pos** dict   (default None)
            overrides self.pos for this display (but does not overwrite 
            self.pos.  Use set_pos if you want to do this)
                
        **statuses_to_plot** list of statuses to plot.
            If given, then the other nodes will be left invisible when plotting
            but I think this requires networkx v2.3 or later.

        ****nx_kwargs**
            any networkx keyword arguments to go into the network plot.
                
            
        '''
#        if not self.SIR and ts_plots:
#            ts_plots = [x for x in ts_plots if x != 'R']
        
        if frame_times is None:
            frame_times = np.linspace(0,self._t_[-1], 101)
        if statuses_to_plot is None:
            statuses_to_plot = self._possible_statuses_
        if ts_plots is None:
            ts_plots = statuses_to_plot
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
            
        if status_order is False:
            status_order = []
            
        #First we draw all of the nodes with their original status, and without
        #putting particular status on top.  All nodes are in place, and their color
        #can be updated at a later time.
        #
        #Then we select the nodes whose status puts them on top initially
        #
        #For each status that goes on top, we draw it in a way that we'll be
        #able to redraw that status at a later time.
        drawn_nodes = self._display_graph_(pos, initial_status, nodelist, False, statuses_to_plot, graph_ax, **nx_kwargs)
        elevated = {status: [node for node in self.G if initial_status[node] == status] for status in status_order}

        drawn_elevated = {}
        for status in reversed(status_order):
            drawn_elevated[status]=[self._draw_specific_status_(pos, elevated[status], status, graph_ax, **nx_kwargs)] #making each a list so that I can change the entry in the list while still passing the same object
        #WARNING I'm defining a dict and while that definition is happening
        #it's drawing things
        
        
        
        
        if ts_plots:
            ts_axes, time_markers = self._display_time_series_(fig, frame_times[0], ts_plots, ts_list, timelabel)
        else:
            ts_axes, time_markers = [], []
        plt.tight_layout()
        
        fargs = (pos, nodelist, drawn_nodes, drawn_elevated, status_order, graph_ax, ts_axes, time_markers, nx_kwargs)

        ani = FuncAnimation(fig, self._update_ani_, frames = frame_times, fargs = fargs, repeat=False)

        return ani
                
        #to show, do
        #simulation.animation()
        #plt.show()
        #to save do
        #ani = simulation.animation()
        #ani.save(filename, fps=5, extra_args=['-vcodec', 'libx264'])

