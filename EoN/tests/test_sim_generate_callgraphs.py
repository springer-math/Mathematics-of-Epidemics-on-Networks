from EoN import *
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

import numpy as np
import networkx as nx


G = nx.grid_2d_graph(50, 50)  # each node is (u,v) where 0<=u,v<=99
initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]

class TestSimSweepParameters:
    graphviz = GraphvizOutput()

    @classmethod
    def setup_class(cls):
        print("setup_class() before any methods in this class")

    @classmethod
    def teardown_class(cls):
        print("teardown_class() after any methods in this class")

    def test_Gillespie_SIS_graph_gen(self):
        TestSimSweepParameters.graphviz.output_file = 'Gillespie_SIS.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I = EoN.Gillespie_SIS(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)

    def test_Gillespie_SIR__graph_gen(self):
        TestSimSweepParameters.graphviz.output_file = 'Gillespie_SIR.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I, R = EoN.Gillespie_SIR(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)

    def test_fast_SIS_sweep_gamma(self):
        TestSimSweepParameters.graphviz.output_file = 'Fast_SIS.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I = EoN.fast_SIS(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)

    def test_fast_SIR_sweep_gamma(self):
        TestSimSweepParameters.graphviz.output_file = 'Fast_SIR.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I, R = EoN.fast_SIR(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)

    def test_basic_discrete_SIS_sweep_p(self):
        TestSimSweepParameters.graphviz.output_file = 'Basic_discrete_SIS.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I = EoN.basic_discrete_SIS(G, 1.0, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True)

    def test_basic_discrete_SIR_sweep_p(self):
        TestSimSweepParameters.graphviz.output_file = 'Basic_discrete_SIR.png'
        with PyCallGraph(output=TestSimSweepParameters.graphviz):
            t, S, I, R = EoN.basic_discrete_SIR(G, 1.0, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True)

