from EoN import *

import numpy as np
import networkx as nx


G = nx.grid_2d_graph(50, 50)  # each node is (u,v) where 0<=u,v<=99
initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]

class TestSimSweepParameters:

    @classmethod
    def setup_class(cls):
        print("setup_class() before any methods in this class")

    @classmethod
    def teardown_class(cls):
        print("teardown_class() after any methods in this class")

    # example of how to run data dependent tests in
    def test_evens(self):
        for i in range(0, 10, 2):
            yield check_even, i, i*3

    def check_even(n, nn):
        assert n % 2 == 0 or nn % 2 == 0

    def test_Gillespie_SIS_type(self):
        t, S, I = EoN.Gillespie_SIS(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=False, tmax=10)

        assert t.__len__() == 3545
        assert S.__len__() == 3545
        assert I.__len__() == 3545

    def test_Gillespie_SIS_sweep_gamma(self):
        [EoN.Gillespie_SIS(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_Gillespie_SIR_sweep_gamma(self):
        [EoN.Gillespie_SIR(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_fast_SIS_sweep_gamma(self):
        [EoN.fast_SIS(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_fast_SIR_sweep_gamma(self):
        [EoN.fast_SIR(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_Gillespie_SIS_sweep_tau(self):
        [EoN.Gillespie_SIS(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_Gillespie_SIR_sweep_tau(self):
        [EoN.Gillespie_SIR(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_fast_SIS_sweep_tau(self):
        [EoN.fast_SIS(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_fast_SIR_sweep_tau(self):
        [EoN.fast_SIR(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]

    def test_basic_discrete_SIS_sweep_p(self):
        [EoN.basic_discrete_SIS(G, i, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True) for i in np.arange(0.0, 1.0, 0.1)]

    def test_basic_discrete_SIR_sweep_p(self):
        [EoN.basic_discrete_SIR(G, i, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True) for i in np.arange(0.0, 1.0, 0.1)]

