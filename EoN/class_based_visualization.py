import matplotlib.pyplot as plt
import networkx as nx
import random
from EoN import *
from matplotlib.animation import FuncAnimation

#from matplotlib.lines import Line2D


#class Animation(object):
#    def __init__():

class EoN_Plot(object):
    def __init__(self, timeseries = ['S', 'I', 'R'], timelabel = r'$t$', SIR = True):
        
        self.timeseries = timeseries

        timeseries_count = len(self.timeseries)
        
        self.timeseries_axi = [None for ts in timeseries] #plural of axes = axi?

        if timeseries_count==0:
            self.fig = plt.figure()
            self.network_axes = self.fig.add_subplot(1, 1, 1)
        else:
            self.fig = plt.figure(figsize = (10,4))
            self.network_axes = self.fig.add_subplot(1, 2, 1)

            self.timeseries.reverse()
            
            ax = self.fig.add_subplot(timeseries_count, 2, 2*(timeseries_count))
            self.timeseries_axi[-1] = ax
            ax.set_xlabel(timelabel)
            for index in range(1,timeseries_count):
                    ax = self.fig.add_subplot(timeseries_count, 2, 2*(timeseries_count-index), sharex=self.timeseries_axi[-1])
                    self.timeseries_axi[-index-1] = ax
                    plt.setp(ax.get_xticklabels(), visible=False)
            self.timeseries.reverse()

    def plot_network(self, G, status, pos=None, nodelist=None, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}, **nx_kwargs):
        r'''
        Plots the network in the axes self.networkx_ax.  Nodes are colored according to their status.
        if no ordered nodelist is provided, then the nodes are plotted so that 'I' nodes appear on top, while the order of the 'S' and 'R' nodes is random.  When the network is very dense this highlights 'I' nodes and allows the final state to more accurately represent the proportion of nodes having each status.
        '''
        colorlist = []
        if pos is None:
            pos = nx.spring_layout(G)
        if nodelist is None:
            nodelist = list(G.nodes())
            I_nodes = [node for node in nodelist if status[node] == 'I']
            other_nodes = [node for node in nodelist if status[node]!='I']
            random.shuffle(other_nodes)
            nodelist = other_nodes + I_nodes
            edgelist = list(G.edges())
        else:
            nodeset = set(nodelist)
            edgelist = [edge for edge in G.edges() if edge[0] in nodeset and edge[1] in nodeset]
        for node in nodelist:
            colorlist.append(colordict[status[node]])
        nx.draw_networkx_edges(G, pos, edgelist=edgelist, ax = self.network_axes, **nx_kwargs)
        nx.draw_networkx_nodes(G, pos, nodelist = nodelist, node_color=colorlist, ax=self.network_axes, **nx_kwargs)
        self.network_axes.set_xticks([])
        self.network_axes.set_yticks([])

    def plot_timeseries(self, t, S=None, I=None, R=None, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}, **kwargs):
        for ax, plot_type in zip(self.timeseries_axi, self.timeseries):
            for label, data in [('S', S), ('I', I), ('R', R)]:
                if label in plot_type and data is not None:
                    ax.plot(t, data, colordict[label], **kwargs)

    def plot_timeseries_from_analytic_model(self, G, tau, gamma, initial_status, model = EBCM_from_graph, tmin = 0, tmax = 10, tcount = 1001, SIR = True, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}, **kwargs):
        r''' Uses one of the analytic models to predict the curve.
        The analytic model needs to be one of the *_from_graph models.

        Arguments:
            G (networkx Graph)
                The graph that we will use to initialize the parameters for the model.

            initial status (dict)
                States the status of each node at t=tmin

            model (function)
                A function like the *_from_graph models

            
        '''
        
        initial_infecteds = [node for node in G.nodes() if initial_status[node]=='I']

        if SIR:
            initial_recovereds = [node for node in G.nodes() if initial_status[node]=='R']
            
            t, S, I, R = model(G, tau, gamma, initial_infecteds=initial_infecteds, 
                    initial_recovereds = initial_recovereds, tmin = tmin, tmax=tmax,
                    tcount=tcount, return_full_data=False)
            
        else:
            t, S, I = model(G, tau, gamma, initial_infecteds=initial_infecteds, 
                    tmin = tmin, tmax=tmax,
                    tcount=tcount, return_full_data=False)
            R=None

        self.plot_timeseries(t, S=S, I=I, R=R, colordict=colordict, **kwargs)
        
    def highlight_time(self, t):
        for ax in self.timeseries_axi:
            ax.axvline(x=t, linestyle='--', color='k')

    def legend(self, timeseries_plot_index = 0, loc='best'):
        #timeseries plot can be:
        #integer from 0 to len(timeseries)-1
        # or string 'All'
        if not self.timeseries_axi:
            raise EoN.EoNError("no time series axes defined")
        elif timeseries_plot_index == 'All':
            for ax in self.timeseries_axi:
                ax.legend(loc=loc)
        else:
            self.timeseries_axi[timeseries_plot_index].legend(loc=loc)

    def savefig(self, filename):
        self.fig.tight_layout()
        self.fig.savefig(filename)

class EoN_Animation(object):
    def __init__(self, G, node_history, frame_times, pos=None, nodelist=None, t=None, S=None, I=None, R=None, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'},  timeseries = ['S', 'I', 'R'], timelabel = r'$t$', SIR = True, **kwargs):
        self.G = G
        self.node_history = node_history
        self.frame_times = frame_times

        self.kwargs = kwargs
        
        if pos is None:
            self.pos = nx.spring_layout(G)
        else:
            self.pos = pos

        if nodelist is None:
            nodelist = list(G.nodes())
            random.shuffle(nodelist)
        self.nodelist=nodelist
        
        self.colordict = colordict

        self.timeseries = timeseries
        timeseries_count = len(self.timeseries)
        
        self.timeseries_axi = [None for ts in self.timeseries] #plural of axes = axi?

        if timeseries_count==0:
            self.fig = plt.figure()
            self.network_axes = self.fig.add_subplot(1, 1, 1)
        else:
            self.fig = plt.figure(figsize = (10,4))
            self.network_axes = self.fig.add_subplot(1, 2, 1)

            self.timeseries.reverse()
            
            ax = self.fig.add_subplot(timeseries_count, 2, 2*(timeseries_count))
            self.timeseries_axi[-1] = ax
            ax.set_xlabel(timelabel)
            for index in range(1,timeseries_count):
                    ax = self.fig.add_subplot(timeseries_count, 2, 2*(timeseries_count-index), sharex=self.timeseries_axi[-1])
                    self.timeseries_axi[-index-1] = ax
                    plt.setp(ax.get_xticklabels(), visible=False)
            self.timeseries.reverse()
            if t is not None:
                self.plot_timeseries(t, S, I, R, label='simulation')

    def plot_timeseries(self, t, S=None, I=None, R=None, colordict={'S':'#009a80','I':'#ff2020', 'R':'gray'}, **kwargs):
        for ax, plot_type in zip(self.timeseries_axi, self.timeseries):
            for label, data in [('S', S), ('I', I), ('R', R)]:
                if label in plot_type and data is not None:
                    ax.plot(t, data, colordict[label], **kwargs)

    def legend(self, timeseries_plot_index = 0, loc='best'):
        #timeseries plot can be:
        #integer from 0 to len(timeseries)-1
        # or string 'All'
        if not self.timeseries_axi:
            raise EoN.EoNError("no time series axes defined")
        elif timeseries_plot_index == 'All':
            for ax in self.timeseries_axi:
                ax.legend(loc=loc)
        else:
            self.timeseries_axi[timeseries_plot_index].legend(loc=loc)

    def _initialize(self):
        initial_status = EoN.get_statuses(self.G, self.node_history, self.frame_times[0])
        colorlist = [self.colordict[initial_status[node]] for node in self.nodelist]

        nodeset = {node for node in self.nodelist}
        edgelist = [edge for edge in self.G.edges() if edge[0] in nodeset and edge[1] in nodeset]
        nx.draw_networkx_edges(self.G, pos=self.pos, edgelist=edgelist, ax = self.network_axes)

        drawn_nodes = nx.draw_networkx_nodes(self.G, pos=self.pos, ax = self.network_axes, nodelist = self.nodelist, color=colorlist, **self.kwargs)
        Inodelist = [node for node in self.nodelist if initial_status[node] == 'I']
        drawn_I = [nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=Inodelist, color = self.colordict['I'], ax = self.network_axes, **self.kwargs)]

        self.network_axes.set_xticks([])
        self.network_axes.set_yticks([])
        
        time_markers = [None for ax in self.timeseries_axi]
        self._highlight_time(self.frame_times[0], time_markers)
        return drawn_nodes, drawn_I, time_markers
        
    
    def _animation_update(self, t, drawn_nodes, drawn_I, time_markers):
        status = EoN.get_statuses(self.G, self.node_history, t)
        Inodelist = [node for node in self.nodelist if status[node] == 'I']
        drawn_nodes.set_color([self.colordict[status[node]] for node in self.nodelist])
        drawn_I[0].remove()
        drawn_I[0] = nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=Inodelist, color = self.colordict['I'], ax = self.network_axes, **self.kwargs)
        for tm in time_markers:
            tm.remove()
        self._highlight_time(t, time_markers)

    def _animation_update_saveframe(self, t, drawn_nodes, drawn_I, time_markers, filename_base, filetype):
        self._animation_update(t, drawn_nodes, drawn_I, time_markers)
        plt.savefig(filename_base + str(t).replace('.', 'p') + filetype)
        

    def save_animation_frames(self, filename_base = 'tmp', filetype = '.png'):
        drawn_nodes, drawn_I,  time_markers = self._initialize()
        fargs = (drawn_nodes, drawn_I, time_markers, filename_base, filetype)
        ani = FuncAnimation(self.fig, self._animation_update_saveframe, frames = self.frame_times, fargs = fargs, repeat=False)
        plt.show()
        
    def save_animation(self, filename_base = 'tmp', filetype = '.mp4', fps = 5, extra_args = ['-vcodec', 'libx264']):
        drawn_nodes, drawn_I, time_markers = self._initialize()
        fargs = (drawn_nodes, drawn_I, time_markers)
        ani = FuncAnimation(self.fig, self._animation_update, frames = self.frame_times, fargs = fargs, repeat=False)
        ani.save(filename_base+filetype, fps=fps, extra_args=extra_args)
        plt.show()
        
    def _highlight_time(self, t, time_markers):
        for index, ax in enumerate(self.timeseries_axi):
            time_markers[index] = ax.axvline(x=t, linestyle='--', color='k')
    
