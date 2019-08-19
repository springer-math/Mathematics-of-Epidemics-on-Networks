import networkx as nx
from EoN import *
import matplotlib.pyplot as plt
import numpy as np
import scipy
import random
from collections import defaultdict

colors = ['#5AB3E6', '#FF2000', '#009A80', '#E69A00', '#CD9AB3', '#0073B3', '#F0E442']

class TestSample:

    @classmethod
    def setup_class(cls):
        print("setup_class() before any methods in this class")

    @classmethod
    def teardown_class(cls):
        print("teardown_class() after any methods in this class")

    def test_discrete_SIR(self):
        print('testing discrete_SIR')
        passed = True
        G = nx.fast_gnp_random_graph(1000, 0.004)

        def test_trans_fxn(u, v):
            '''
            Transmissions occur if one odd and one even.  So it is basically bipartite.
            '''
            if (u + v) % 2 == 0:
                return False
            else:
                return True

        sim = EoN.discrete_SIR(G, test_transmission=test_trans_fxn, args=(), initial_infecteds=[0, 2, 4, 6, 8, 10],
                               return_full_data=True)
        # by initial condition and transmission rule, infection generations alternate parity.
        for node in G:
            if 'I' in sim.node_history(node)[1]:
                idx = sim.node_history(node)[1].index('I')

                if (node + sim.node_history(node)[0][idx]) % 2 == 1:  # should always be False
                    print('Error', node, sim.node_history(node))
                    passed = False
        print('number infected', sim.R()[-1])
        if not passed:
            print('failed')
        else:
            print('passed')
        assert passed


    def test_basic_discrete_SIR(self):
        print('testing basic_discrete_SIR, percolation_based_discrete_SIR, and EBCM_discrete_from_graph')
        plt.clf()
        N = 1000000
        initial_count = 5000
        G = nx.fast_gnp_random_graph(N, 4. / N)
        p = 0.4
        sim = EoN.basic_discrete_SIR(G, p, initial_infecteds=range(initial_count), return_full_data=True)
        t, S, I, R = sim.summary()
        sim2 = EoN.percolation_based_discrete_SIR(G, p, initial_infecteds=range(initial_count), return_full_data=True)
        t2, S2, I2, R2 = sim2.summary()
        t3, S3, I3, R3 = EoN.EBCM_discrete_from_graph(G, p, rho=float(initial_count) / N)
        t4, S4, I4, R4 = EoN.EBCM_discrete_from_graph(G, p, initial_infecteds=range(initial_count))
        print(t)
        print(S)
        print(I)
        print(R)
        print(t[0:4], S[0:4], I[0:4], R[0:4])
        print(t2[0:4], S2[0:4], I2[0:4], R2[0:4])
        print(t3[0:4], S3[0:4], I3[0:4], R3[0:4])
        print(t4[0:4], S4[0:4], I4[0:4], R4[0:4])
        plt.plot(t, S, label='basic sim', alpha=0.3)
        plt.plot(t, I, alpha=0.3)
        plt.plot(t, R, alpha=0.3)
        plt.plot(t2, S2, '--', label='percolation based', alpha=0.3)
        plt.plot(t2, I2, '--', alpha=0.3)
        plt.plot(t2, R2, '--', alpha=0.3)
        plt.plot(t3, S3, '-.', label='Discrete EBCM', alpha=0.3)
        plt.plot(t3, I3, '-.', alpha=0.3)
        plt.plot(t3, R3, '-.', alpha=0.3)
        plt.plot(t4, S4, ':', label='Discrete EBCM 2', alpha=0.3)
        plt.plot(t4, I4, ':', alpha=0.3)
        plt.plot(t4, R4, ':', alpha=0.3)
        plt.legend(loc='upper right')
        filename = 'basic_discrete_SIR_test'
        plt.savefig(filename)
        print("check {} for good match".format(filename))


    def test_estimate_SIR_prob_size(self):
        print('testing estimate_SIR_prob_size')
        N = 1000000
        G = nx.fast_gnp_random_graph(N, 5. / N)
        for p in scipy.linspace(0.1, 0.5, 5):
            P, A = EoN.estimate_SIR_prob_size(G, p)
            gamma = 1.
            tau = p * gamma / (1 - p)
            P2, A2 = EoN.estimate_directed_SIR_prob_size(G, tau, 1.0)
            t, S, I, R = EoN.EBCM_discrete_from_graph(G, p)
            print("should all be approximately the same: ", R[-1] / G.order(), A, A2)


    def test_SIR_dynamics(self):
        print("test_SIR_dynamics")
        plt.clf()
        reduced_report = scipy.linspace(0, 15, 31)
        G = nx.configuration_model([1, 5, 10] * 100000)
        N = G.order()
        initial_size = 10000
        gamma = 1.
        tau = 0.3
        t, S, I, R = EoN.fast_SIR(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, label='fast_SIR', alpha=0.3)
        plt.plot(t, I, alpha=0.3)
        plt.plot(t, R, alpha=0.3)

        t, S, I, R = EoN.Gillespie_SIR(G, tau, gamma, initial_infecteds=range(initial_size))
        plt.plot(t, S, '--', label='Gillespie_SIR', alpha=0.3)
        plt.plot(t, I, '--', alpha=0.3)
        plt.plot(t, R, '--', alpha=0.3)

        t, S, I, R = EoN.EBCM_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, ':', label='EBCM', alpha=0.3)
        plt.plot(t, I, ':', alpha=0.3)
        plt.plot(t, R, ':', alpha=0.3)
        t, S, I, R = EoN.EBCM_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, 'x', label='EBCM', alpha=0.3)
        plt.plot(t, I, 'x', alpha=0.3)
        plt.plot(t, R, 'x', alpha=0.3)

        t, S, I, R = EoN.SIR_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '-.', label='compact pairwise', alpha=0.3)
        plt.plot(t, I, '-.', alpha=0.3)
        plt.plot(t, R, '-.', alpha=0.3)
        t, S, I, R = EoN.SIR_compact_pairwise_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, 's', label='compact pairwise', alpha=0.3)
        plt.plot(t, I, 's', alpha=0.3)
        plt.plot(t, R, 's', alpha=0.3)

        t, S, I, R = EoN.SIR_super_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '.', label='super compact pairwise', alpha=0.3)
        plt.plot(t, I, '.', alpha=0.3)
        plt.plot(t, R, '.', alpha=0.3)
        t, S, I, R = EoN.SIR_super_compact_pairwise_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, 'd', label='super compact pairwise', alpha=0.3)
        plt.plot(t, I, 'd', alpha=0.3)
        plt.plot(t, R, 'd', alpha=0.3)

        t, S, I, R = EoN.SIR_effective_degree_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, 'x', label='effective degree', alpha=0.3)
        plt.plot(t, I, 'x', alpha=0.3)
        plt.plot(t, R, 'x', alpha=0.3)
        t, S, I, R = EoN.SIR_effective_degree_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, '^', label='effective degree', alpha=0.3)
        plt.plot(t, I, '^', alpha=0.3)
        plt.plot(t, R, '^', alpha=0.3)

        t, S, I, R = EoN.SIR_compact_effective_degree_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '+', label='compact effective degree', alpha=0.3)
        plt.plot(t, I, '+', alpha=0.3)
        plt.plot(t, R, '+', alpha=0.3)
        t, S, I, R = EoN.SIR_compact_effective_degree_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, 'v', label='compact effective degree', alpha=0.3)
        plt.plot(t, I, 'v', alpha=0.3)
        plt.plot(t, R, 'v', alpha=0.3)

        t, S, I, R = EoN.SIR_heterogeneous_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, 'o', label='heterogeneous pairwise', alpha=0.3)
        plt.plot(t, I, 'o', alpha=0.3)
        plt.plot(t, R, 'o', alpha=0.3)
        t, S, I, R = EoN.SIR_heterogeneous_pairwise_from_graph(G, tau, gamma, initial_infecteds=range(initial_size))
        S, I, R = EoN.subsample(reduced_report, t, S, I, R)
        t = reduced_report
        plt.plot(t, S, '>', label='heterogeneous pairwise', alpha=0.3)
        plt.plot(t, I, '>', alpha=0.3)
        plt.plot(t, R, '>', alpha=0.3)

        # t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = float(initial_size)/N)
        # plt.plot(t, S, '-.', label = 'heterogeneous meanfield')
        # plt.plot(t, I, '-.')
        # plt.plot(t, R, '-.')
        plt.axis(xmin=-20, xmax=15)
        plt.legend(loc='center left')
        plt.savefig('SIR_dynamics')


    def test_SIS_dynamics(self):
        print("test_SIS_dynamics")
        plt.clf()
        reduced_report = scipy.linspace(0, 15, 31)

        G = nx.configuration_model([1, 5, 10] * 100000)
        G = nx.Graph(G)
        G.remove_edges_from(G.selfloop_edges())
        N = G.order()
        initial_size = 5000
        gamma = 1.
        tau = 0.3

        print('\tfast_SIS')
        t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds=range(initial_size), tmax=15)
        plt.plot(t, S, label='fast_SIS', alpha=0.3)
        plt.plot(t, I, alpha=0.3)
        print('\tfast_nonMarkov_SIS')

        def trans_time_fxn(source, target, rec_delay, tau):
            r = []
            d = random.expovariate(tau)
            while d < rec_delay:
                r.append(d)
                d += random.expovariate(tau)
            return r

        def rec_time_fxn(u, gamma):
            return random.expovariate(gamma)

        t, S, I = EoN.fast_nonMarkov_SIS(G, trans_time_fxn, rec_time_fxn, trans_time_args=(tau,), rec_time_args=(gamma,),
                                         initial_infecteds=range(initial_size), tmax=15)

        plt.plot(t, S, ':', label='fast_nonMarkov_SIS', alpha=0.3)
        plt.plot(t, I, ':', alpha=0.3)

        print('\tGillespie_SIS')
        t, S, I = EoN.Gillespie_SIS(G, tau, gamma, initial_infecteds=range(initial_size), tmax=15)
        plt.plot(t, S, '--', label='Gillespie_SIS', alpha=0.3)
        plt.plot(t, I, '--', alpha=0.3)

        print('\tSIS_compact_pairwise')
        t, S, I = EoN.SIS_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '-.', label='compact pairwise', alpha=0.3)
        plt.plot(t, I, '-.', alpha=0.3)
        t, S, I = EoN.SIS_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        S, I = EoN.subsample(reduced_report, t, S, I)
        plt.plot(reduced_report, S, 's', label='compact pairwise', alpha=0.3)
        plt.plot(reduced_report, I, 's', alpha=0.3)

        print('\tSIS_super_compact_pairwise')
        t, S, I = EoN.SIS_super_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, 'o', label='super compact pairwise', alpha=0.3)
        plt.plot(t, I, 'o', alpha=0.3)
        t, S, I = EoN.SIS_super_compact_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        S, I = EoN.subsample(reduced_report, t, S, I)
        plt.plot(reduced_report, S, 'd', label='super compact pairwise', alpha=0.3)
        plt.plot(reduced_report, I, 'd', alpha=0.3)

        print('\tSIS_effective_degree')
        t, S, I = EoN.SIS_effective_degree_from_graph(G, tau, gamma, tmax=15, rho=float(initial_size) / N)
        plt.plot(t, S, 'x', label='effective degree', alpha=0.3)
        plt.plot(t, I, 'x', alpha=0.3)
        t, S, I = EoN.SIS_effective_degree_from_graph(G, tau, gamma, tmax=15, rho=float(initial_size) / N)
        S, I = EoN.subsample(reduced_report, t, S, I)
        plt.plot(reduced_report, S, '^', label='effective degree', alpha=0.3)
        plt.plot(reduced_report, I, '^', alpha=0.3)

        print('\tSIS_compact_effective_degree')
        t, S, I = EoN.SIS_compact_effective_degree_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '+', label='compact effective degree', alpha=0.3)
        plt.plot(t, I, '+', alpha=0.3)
        t, S, I = EoN.SIS_compact_effective_degree_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        S, I = EoN.subsample(reduced_report, t, S, I)
        plt.plot(reduced_report, S, '>', label='compact effective degree', alpha=0.3)
        plt.plot(reduced_report, I, '>', alpha=0.3)

        print('\tSIS_heterogeneous_pairwise')
        t, S, I = EoN.SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        plt.plot(t, S, '.', label='heterogeneous pairwise', alpha=0.3)
        plt.plot(t, I, '.', alpha=0.3)
        t, S, I = EoN.SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho=float(initial_size) / N)
        S, I = EoN.subsample(reduced_report, t, S, I)
        plt.plot(reduced_report, S, 'v', label='heterogeneous pairwise', alpha=0.3)
        plt.plot(reduced_report, I, 'v', alpha=0.3)

        # t, S, I = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = float(initial_size)/N)
        # plt.plot(t, S, '-.', label = 'heterogeneous meanfield')
        # plt.plot(t, I, '-.')
        plt.axis(ymin=0, xmax=15)
        plt.legend()
        print('saving figure')
        plt.savefig('SIS_dynamics.pdf')


    def test_SIS_simulations(self):
        print("test_SIS_simulations")
        plt.clf()
        tau = 0.1
        gamma = 0.3
        G = nx.configuration_model([1, 5, 10] * 100000)
        G = nx.Graph(G)
        G.remove_edges_from(G.selfloop_edges())
        N = G.order()
        initial_size = 5000
        for counter in range(10):
            print('fast_SIS')
            t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds=range(initial_size), tmax=20)
            plt.plot(t, S, '-.', color='b', alpha=0.3)
            plt.plot(t, I, '-.', color='b', alpha=0.3)
            print('Gillespie_SIS')
            t, S, I = EoN.Gillespie_SIS(G, tau, gamma, initial_infecteds=range(initial_size), tmax=20)
            plt.plot(t, S, '--', color='r', alpha=0.3)
            plt.plot(t, I, '--', color='r', alpha=0.3)
            plt.title('curves should overlie to show event driven and gillespie agree')
            plt.savefig('SIS_sims')


    def test_SIR_final_sizes():
        print("test_SIR_final_sizes")
        plt.clf()
        G = nx.configuration_model([3, 6, 3, 6, 20] * 10000)
        N = G.order()
        tau = 0.2
        gamma = 1
        t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds=range(5000))
        plt.plot(t, S, color=colors[0], label='simulation', alpha=0.3)
        plt.plot(t, I, color=colors[0], alpha=0.3)
        plt.plot(t, R, color=colors[0], alpha=0.3)
        infected_nodes = EoN.get_infected_nodes(G, tau, gamma, initial_infecteds=range(5000))
        A = len(infected_nodes)
        print(A)
        plt.plot([0, 10], [A, A], label=r'percolation based', alpha=0.3)
        A = EoN.Attack_rate_cts_time_from_graph(G, tau, gamma, rho=0.1)
        plt.plot([0, 10], [A * N, A * N], '-.', label='analytic', alpha=0.3)
        plt.legend(loc='upper right')
        plt.savefig('test_SIR_final_sizes')


    # ======================================================================
    # ERROR: EoN.tests.test_from_joel.test_SIR_individual_based
    # ----------------------------------------------------------------------
    # Traceback (most recent call last):
    #   File "c:\users\tting\appdata\local\programs\python\python36\lib\site-packages\nose\case.py", line 198, in runTest
    #     self.test(*self.arg)
    #   File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_from_joel.py", line 347, in test_SIR_individual_based
    #     t, S, I, R = EoN.SIR_individual_based(G, nodelist, X, Y, tau, gamma=gamma, tmax=20)
    # TypeError: SIR_individual_based() got multiple values for argument 'gamma'
    # test_SIR_individual_based
    # def test_SIR_individual_based():
    #     print("test_SIR_individual_based")
    #     plt.clf()
    #
    #     G = nx.configuration_model([3, 10] * 10000)
    #     tau = 0.3
    #     gamma = 1
    #     N = G.order()
    #
    #     initial_infecteds = scipy.random.choice(G.nodes(), size=100, replace=False)
    #     rho = len(initial_infecteds) * 1. / N
    #
    #     t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds=initial_infecteds)
    #     plt.plot(t, S, '--', label='simulation', alpha=0.3)
    #     plt.plot(t, I, '--', alpha=0.3)
    #     plt.plot(t, R, '--', alpha=0.3)
    #     nodelist = G.nodes()
    #     X = (1 - rho) * scipy.ones(N)
    #     Y = rho * scipy.ones(N)
    #
    #     t, S, I, R = EoN.SIR_individual_based(G, nodelist, X, Y, tau, gamma=gamma, tmax=20)
    #
    #     plt.plot(t, S, label='individual-based equations', alpha=0.3)
    #     plt.plot(t, I, alpha=0.3)
    #     plt.plot(t, R, alpha=0.3)
    #     plt.legend(loc='upper right', alpha=0.3)
    #     plt.savefig('SIR_individual_based')

    # ======================================================================
    # ERROR: EoN.tests.test_from_joel.test_SIS_individual_based
    # ----------------------------------------------------------------------
    # Traceback (most recent call last):
    #   File "c:\users\tting\appdata\local\programs\python\python36\lib\site-packages\nose\case.py", line 198, in runTest
    #     self.test(*self.arg)
    #   File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_from_joel.py", line 363, in test_SIS_individual_based
    #     t, S, I = EoN.SIS_individual_based(G, nodelist, Y, 0.3, gamma=1, tmax=20)
    # TypeError: SIS_individual_based() got multiple values for argument 'gamma'
    # def test_SIS_individual_based():
    #     print('test_SIS_individual_based')
    #     G = nx.configuration_model([3, 10] * 1000)
    #     nodelist = G.nodes()
    #     N = G.order()
    #     rho = 1. / N
    #     Y = rho * scipy.ones(N)
    #     t, S, I = EoN.SIS_individual_based(G, nodelist, Y, 0.3, gamma=1, tmax=20)
    #     plt.clf()
    #     plt.plot(t, S)
    #     plt.plot(t, I)
    #     plt.savefig('SIS_individual_based')


    def test_pair_based(self):
        print("test_pair_based")
        G = nx.fast_gnp_random_graph(1000, 0.004)
        nodelist = G.nodes()
        Y0 = scipy.array([1 if node < 10 else 0 for node in nodelist])

        print('testing SIS_pair_based')
        t, S, I = EoN.SIS_pair_based(G, 2, 0.5, rho=0.01, tmax=5, tcount=101)
        print('creating SIS fig')
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.savefig('SIS_pair_based')

        print('testing SIR_pair_based')
        t, S, I, R = EoN.SIR_pair_based(G, 2, 0.5, tmax=5, tcount=101)
        print('creating SIR fig')
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.savefig('SIR_pair_based')

    # ======================================================================
    # ERROR: EoN.tests.test_from_joel.test_SIR_pair_based2
    # ----------------------------------------------------------------------
    # Traceback (most recent call last):
    #   File "c:\users\tting\appdata\local\programs\python\python36\lib\site-packages\nose\case.py", line 198, in runTest
    #     self.test(*self.arg)
    #   File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_from_joel.py", line 401, in test_SIR_pair_based2
    #     t, S, I, R = EoN.SIR_pair_based2(G, tau, gamma=gamma, nodelist=nodelist, Y0=Y0, tmax=5, tcount=101)
    # AttributeError: module 'EoN' has no attribute 'SIR_pair_based2'
    # def test_SIR_pair_based2():
    #     print("test_SIR_pair_based2")
    #     G = nx.fast_gnp_random_graph(1000, 0.004)
    #     nodelist = G.nodes()
    #     Y0 = scipy.array([1 if node < 10 else 0 for node in nodelist])
    #     tau = 2
    #     gamma = 0.5
    #     t, S, I, R = EoN.SIR_pair_based2(G, tau, gamma=gamma, nodelist=nodelist, Y0=Y0, tmax=5, tcount=101)
    #     print('creating SIR fig2')
    #     plt.clf()
    #     plt.plot(t, S)
    #     plt.plot(t, I)
    #     plt.plot(t, R)
    #     plt.savefig('SIR_pair_based2')
    #     print('done with pair_based2')


    def test_SIS_homogeneous_meanfield(self):
        print("testing SIS_homogeneous_meanfield")
        S0 = 99
        I0 = 1
        n = 1
        tau = 3
        gamma = 2
        t, S, I = EoN.SIS_homogeneous_meanfield(S0, I0, n, tau, gamma, tmax=10)

        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.savefig('SIS_homogeneous_meanfield')


    def test_SIR_homogeneous_meanfield(self):
        print("testing SIR_homogeneous_meanfield")
        S0 = 97
        I0 = 2
        R0 = 1
        n = 1
        tau = 3
        gamma = 2
        t, S, I, R = EoN.SIR_homogeneous_meanfield(S0, I0, R0, n, tau, gamma, tmax=10)
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.savefig('SIR_homogeneous_meanfield')


    def test_SIS_homogeneous_pairwise(self):
        print("test_SIS_homogeneous_pairwise")
        S0 = 99
        I0 = 1
        SI0 = 10
        SS0 = 980
        n = 10
        tau = 1
        gamma = 5
        t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmax=5)
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.savefig('SIS_homogeneous_pairwise')


    def test_SIR_homogeneous_pairwise(self):
        print("test_SIR_homogeneous_pairwise")
        S0 = 97
        I0 = 2
        R0 = 1
        SI0 = 18
        SS0 = 941
        n = 10
        tau = 0.4
        gamma = 1
        t, S, I, R = EoN.SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmax=10)
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.savefig('SIR_homogeneous_pairwise')


    def test_SIS_heterogeneous_meanfield(self):
        print("testing SIS_heterogeneous_meanfield")
        Sk0 = scipy.arange(100) * 100
        Ik0 = scipy.arange(100)

        t, S, I = EoN.SIS_heterogeneous_meanfield(Sk0, Ik0, 1, 10, tmax=1)
        print("plotting SIS_heterogeneous_meanfield")
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.savefig('SIS_heterogeneous_meanfield')


    def test_SIR_heterogeneous_meanfield(self):
        print("testing SIR_heterogeneous_meanfield")
        Sk0 = scipy.arange(100) * 100
        Ik0 = scipy.arange(100)
        Rk0 = 0 * scipy.arange(100)

        t, S, I, R = EoN.SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, 0.1, 5, tmax=5)
        print("plotting SIR_heterogeneous_meanfield")
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.savefig('SIR_heterogeneous_meanfield')


    def test_SIS_heterogeneous_pairwise(self):
        print("test_SIS_heterogeneous_pairwise")

        # graph will be 2 stars: both with 3 "leaves".  one of them has central node infected
        SkSl0 = scipy.array([[0, 0, 0, 0], [0, 0, 0, 3], [0, 0, 0, 0], [0, 3, 0, 0]])
        SkIl0 = scipy.array([[0, 0, 0, 0], [0, 0, 0, 3], [0, 0, 0, 0], [0, 0, 0, 0]])
        IkIl0 = scipy.zeros((4, 4))

        print((SkSl0 + SkIl0).T / (scipy.array([1, 0, 0, 0]) + scipy.arange(4.)))
        Sk0 = sum((SkSl0 + SkIl0).T / (scipy.array([1, 0, 0, 0]) + scipy.arange(4.)))
        Ik0 = sum((SkIl0.T + IkIl0).T / (scipy.array([1, 0, 0, 0]) + scipy.arange(4.)))

        Sk0[0] = 1
        print('Sk0', Sk0)
        print('Ik0', Ik0)
        tau = 3
        gamma = 1
        #    print(SkIl0, SkSl0
        #    print(Sk0
        #    print(Ik0
        t, S, I = EoN.SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, gamma, tmax=10)
        plt.clf()
        plt.plot(t, S, label='pure IC')
        plt.plot(t, I)

        G = nx.Graph()
        G.add_edges_from([(1, 2), (1, 3), (1, 4), (5, 6), (5, 7), (5, 8)])
        G.add_node(0)
        t, S, I = EoN.SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho=1. / 9, tmax=10)
        plt.plot(t, S, '-.', label='uniform')
        plt.plot(t, I, '-.')
        plt.legend(loc='upper right')
        plt.title('starting from different IC')
        plt.savefig('SIS_heterogeneous_pairwise')


    def test_SIR_heterogeneous_pairwise(self):
        print("SIR_heterogeneous_pairwise not yet tested")


    def test_SIS_compact_pairwise(self):
        print("testing SIS_compact_pairwise")
        EoN.EoNError('changing order of arguments')
        Sk0 = scipy.arange(100) * 100
        Ik0 = scipy.arange(100)
        SI0 = Ik0.dot(scipy.arange(100))
        SS0 = Sk0.dot(scipy.arange(100)) - SI0
        II0 = 0
        tau = 0.1
        gamma = 0.3
        t, S, I = EoN.SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmax=5)
        print("plotting SIS_compact_pairwise")
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.savefig('SIS_compact_pairwise')


    def test_SIR_compact_pairwise(self):
        EoN.EoNError('changing order of arguments')
        print("testing SIR_compact_pairwise")
        Sk0 = scipy.arange(100) * 100
        I0 = sum(scipy.arange(100))
        R0 = 0
        SI0 = 1000
        SS0 = Sk0.dot(scipy.arange(100)) - SI0
        tau = 0.1
        gamma = 0.3
        t, S, I, R = EoN.SIR_compact_pairwise(Sk0, I0, R0, SI0, SS0, tau, gamma, tmax=5)
        print("plotting SIR_compact_pairwise")
        plt.clf()
        plt.plot(t, S)
        plt.plot(t, I)
        plt.plot(t, R)
        plt.savefig('SIR_compact_pairwise')


    def test_SIS_super_compact_pairwise(self):
        print("SIS_super_compact_pairwise not yet tested")


    def test_SIR_super_compact_pairwise(self):
        print("SIR_super_compact_pairwise not yet tested")


    def test_SIS_effective_degree(self):
        print("SIS_effective_degree not yet tested")


    def test_SIR_effective_degree(self):
        print("SIR_effective_degree not yet tested")


    def test_SIS_compact_effective_degree(self):
        print("SIS_compact_effective_degree not yet tested")


    def test_SIR_compact_effective_degree(self):
        print("SIR_compact_effective_degree not yet tested")


    def test_ErdősRényi_million_Fast_Gillespie_SIR(self):
        print("testing ErdősRényi_million_Fast_Gillespie_SIR")
        N = 10 ** 6  # number of individuals
        kave = 5  # expected number of partners
        G = nx.fast_gnp_random_graph(N, kave / (N - 1))  # Erdős-Rényi graph

        rho = 0.005  # initial fraction infected
        tau = 0.3  # transmission rate
        gamma = 1.0  # recovery rate
        t1, S1, I1, R1 = EoN.fast_SIR(G, tau, gamma, rho=rho)
        t2, S2, I2, R2 = EoN.Gillespie_SIR(G, tau, gamma, rho=rho)

        plt.plot(t1, I1, label='fast_SIR')
        plt.plot(t2, I2, label='Gillespie_SIR')
        plt.legend()
        plt.savefig('test_ErdősRényi_million_Fast_Gillespie_SIR')


    def test_ErdősRényi_million_Fast_Gillespie_SIS(self):
        N = 10 ** 5  # number of individuals
        kave = 5  # expected number of partners
        G = nx.fast_gnp_random_graph(N, kave / (N - 1))  # Erdős-Rényi graph

        rho = 0.005  # initial fraction infected
        tau = 0.3  # transmission rate
        gamma = 1.0  # recovery rate
        t1, S1, I1 = EoN.fast_SIS(G, tau, gamma, rho=rho, tmax=30)
        t2, S2, I2 = EoN.Gillespie_SIS(G, tau, gamma, rho=rho, tmax=30)

        plt.plot(t1, I1, label='fast_SIS')
        plt.plot(t2, I2, label='Gillespie_SIS')
        plt.legend()
        plt.savefig('test_ErdősRényi_million_Fast_Gillespie_SIS')


    def test_fast_nonMarkov_SIR(self):
        def rec_time_fxn_gamma(u):
            # gamma(shape, scale = 1.0)
            return np.random.gamma(3, 0.5)

        def trans_time_fxn(u, v, tau):
            if tau > 0:
                return np.random.exponential(1. / tau)
            else:
                return float('Inf')

        N = 10 ** 6  # number of individuals
        kave = 5  # expected number of partners
        G = nx.fast_gnp_random_graph(N, kave / (N - 1))  # Erdős-Rényi graph
        tau = 0.3

        for cntr in range(10):
            t, S, I, R = EoN.fast_nonMarkov_SIR(G, trans_time_fxn=trans_time_fxn,
                                                rec_time_fxn=rec_time_fxn_gamma, trans_time_args=(tau,))
            plt.plot(t, R)

        plt.savefig('test_fast_nonMarkov_SIR')


    def test_SIS_Pairwise_Model(self):
        N = 10000
        gamma = 1
        rho = 0.05
        kave = 20
        tau = 2 * gamma / kave
        S0 = (1 - rho) * N
        I0 = rho * N
        SI0 = (1 - rho) * kave * rho * N
        SS0 = (1 - rho) * kave * (1 - rho) * N
        t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, kave, tau, gamma,
                                               tmax=10)
        plt.plot(t, S, label='S')
        plt.plot(t, I, label='I')
        plt.legend()
        plt.savefig('test_SIS_Pairwise_Model')


    def test_SIR_EBCM(self):
        gamma = 1
        tau = 1.5
        kave = 3
        rho = 0.01
        phiS0 = 1 - rho

        def psi(x):
            return (1 - rho) * np.exp(-kave * (1 - x))

        def psiPrime(x):
            return (1 - rho) * kave * np.exp(-kave * (1 - x))

        N = 1

        t, S, I, R = EoN.EBCM(N, psi, psiPrime, tau, gamma, phiS0, tmax=10)

        plt.plot(t, S, label='S')
        plt.plot(t, I, label='I')
        plt.plot(t, R, label='R')
        plt.legend()
        plt.savefig('test_SIR_EBCM')


    def test_Gillespie_simple_contagion(self):
        N = 100000
        G = nx.fast_gnp_random_graph(N, 5. / (N - 1))

        # they will vary in the rate of leaving exposed class.
        # and edges will vary in transition rate.
        # there is no variation in recovery rate.

        node_attribute_dict = {node: 0.5 + random.random() for node in G.nodes()}
        edge_attribute_dict = {edge: 0.5 + random.random() for edge in G.edges()}

        nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
        nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')
        #
        # These individual and partnership attributes will be used to scale
        # the transition rates.  When we define `H` and `J`, we provide the name
        # of these attributes.

        # We show how node and edge attributes in the contact network 'G' can be used
        # to scale the transmission rates.  More advanced techniques are shown in
        # the documentation

        H = nx.DiGraph()
        H.add_node('S')  # This line is actually unnecessary.
        H.add_edge('E', 'I', rate=0.6, weight_label='expose2infect_weight')
        H.add_edge('I', 'R', rate=0.1)

        J = nx.DiGraph()
        J.add_edge(('I', 'S'), ('I', 'E'), rate=0.1, weight_label='transmission_weight')
        IC = defaultdict(lambda: 'S')
        for node in range(200):
            IC[node] = 'I'

        return_statuses = ('S', 'E', 'I', 'R')

        t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses,
                                                       tmax=float('Inf'))

        plt.plot(t, S, label='Susceptible')
        plt.plot(t, E, label='Exposed')
        plt.plot(t, I, label='Infected')
        plt.plot(t, R, label='Recovered')
        plt.legend()
        plt.savefig('test_Gillespie_simple_contagion')


    def test_Two_Cooperative_SIR_Diseases_oscillatory(self):
        # N = 1000000
        N = 100000
        G = nx.fast_gnp_random_graph(N, 5. / (N - 1))

        print('got G')

        # In the below:
        # 'SS' means an individual susceptible to both diseases
        # 'SI' means susceptible to disease 1 and infected with disease 2
        # 'RS' means recovered from disease 1 and susceptible to disease 2.
        # etc.

        H = nx.DiGraph()  # DiGraph showing possible transitions that don't require an interaction
        H.add_node('SS')  # we actually don't need to include the 'SS' node in H.
        H.add_edge('SI', 'SR', rate=1)
        H.add_edge('IS', 'RS', rate=1)
        H.add_edge('II', 'IR', rate=1)
        H.add_edge('II', 'RI', rate=1)
        H.add_edge('IR', 'RR', rate=0.5)
        H.add_edge('RI', 'RR', rate=0.5)

        # In the below the edge (('SI', 'SS'), ('SI', 'SI')) means an
        # 'SI' individual connected to an 'SS' individual can lead to a transition in which
        # the 'SS' individual becomes 'SI'.  The rate of this transition is 0.2.
        #
        # Note that `IR` and `RI` individuals are more infectious than other individuals.
        #
        J = nx.DiGraph()  # DiGraph showing transitiona that do require an interaction.
        J.add_edge(('SI', 'SS'), ('SI', 'SI'), rate=0.2)
        J.add_edge(('SI', 'IS'), ('SI', 'II'), rate=0.2)
        J.add_edge(('SI', 'RS'), ('SI', 'RI'), rate=0.2)
        J.add_edge(('II', 'SS'), ('II', 'SI'), rate=0.2)
        J.add_edge(('II', 'IS'), ('II', 'II'), rate=0.2)
        J.add_edge(('II', 'RS'), ('II', 'RI'), rate=0.2)
        J.add_edge(('RI', 'SS'), ('RI', 'SI'), rate=1)
        J.add_edge(('RI', 'IS'), ('RI', 'II'), rate=1)
        J.add_edge(('RI', 'RS'), ('RI', 'RI'), rate=1)
        J.add_edge(('IS', 'SS'), ('IS', 'IS'), rate=0.2)
        J.add_edge(('IS', 'SI'), ('IS', 'II'), rate=0.2)
        J.add_edge(('IS', 'SR'), ('IS', 'IR'), rate=0.2)
        J.add_edge(('II', 'SS'), ('II', 'IS'), rate=0.2)
        J.add_edge(('II', 'SI'), ('II', 'II'), rate=0.2)
        J.add_edge(('II', 'SR'), ('II', 'IR'), rate=0.2)
        J.add_edge(('IR', 'SS'), ('IR', 'IS'), rate=1)
        J.add_edge(('IR', 'SI'), ('IR', 'II'), rate=1)
        J.add_edge(('IR', 'SR'), ('IR', 'IR'), rate=1)

        return_statuses = ('SS', 'SI', 'SR', 'IS', 'II', 'IR', 'RS', 'RI', 'RR')

        # initial_size = 650
        initial_size = 65
        IC = defaultdict(lambda: 'SS')
        for individual in range(initial_size):
            IC[individual] = 'II'

        print('got IC')

        t, SS, SI, SR, IS, II, IR, RS, RI, RR = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses,
                                                                               tmax=float('Inf'))

        plt.semilogy(t, IS + II + IR, '-.', label='Infected with disease 1')
        plt.semilogy(t, SI + II + RI, '-.', label='Infected with disease 2')

        plt.legend()
        plt.savefig('test_Two_Cooperative_SIR_Diseases_oscillatory')


    def test_Gillespie_complex_contagion(self):
        def transition_rate(G, node, status, parameters):
            # this function needs to return the rate at which ``node`` changes status
            #
            r = parameters[0]
            if status[node] == 'S' and len([nbr for nbr in G.neighbors(node) if status[nbr] == 'I']) > 1:
                return 1
            else:  # status[node] might be 0 or length might be 0 or 1.
                return 0

        def transition_choice(G, node, status, parameters):
            # this function needs to return the new status of node.  We assume going
            # in that we have already calculated it is changing status.
            #
            # this function could be more elaborate if there were different
            # possible transitions that could happen.  However, for this model,
            # the 'I' nodes aren't changing status, and the 'S' ones are changing to 'I'
            # So if we're in this function, the node must be 'S' and becoming 'I'
            #
            return 'I'

        def get_influence_set(G, node, status, parameters):
            # this function needs to return any node whose rates might change
            # because ``node`` has just changed status.  That is, which nodes
            # might ``node`` influence?
            #
            # For our models the only nodes a node might affect are the susceptible neighbors.

            return {nbr for nbr in G.neighbors(node) if status[nbr] == 'S'}

        parameters = (2,)  # this is the threshold.  Note the comma.  It is needed
        # for python to realize this is a 1-tuple, not just a number.
        # ``parameters`` is sent as a tuple so we need the comma.

        N = 60000
        deg_dist = [2, 4, 6] * int(N / 3)
        G = nx.configuration_model(deg_dist)

        for rho in np.linspace(3. / 80, 7. / 80, 8):  # 8 values from 3/80 to 7/80.
            print(rho)
            IC = defaultdict(lambda: 'S')
            for node in G.nodes():
                if np.random.random() < rho:  # there are faster ways to do this random selection
                    IC[node] = 'I'

            t, S, I = EoN.Gillespie_complex_contagion(G, transition_rate, transition_choice,
                                                      get_influence_set, IC, return_statuses=('S', 'I'),
                                                      parameters=parameters)

            plt.plot(t, I)

        plt.savefig('test_Gillespie_complex_contagion')


    def test_Snapshot_Dynamics_And_TransmissionTree(self):
        G = nx.karate_club_graph()

        nx_kwargs = {"with_labels": True}
        sim = EoN.Gillespie_SIR(G, 1, 1, return_full_data=True)
        sim.display(time=1, **nx_kwargs)
        plt.savefig('test_Snapshot_Dynamics_And_TransmissionTree_1')

        T = sim.transmission_tree()
        Tpos = EoN.hierarchy_pos(T)

        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        nx.draw(T, Tpos, ax=ax, node_size=200, with_labels=True)
        plt.savefig('test_Snapshot_Dynamics_And_TransmissionTree_2')


    def test_Animation_Dynamics_SIR_With_Vaccination_In_Lattice(self):
        G = nx.grid_2d_graph(100, 100)  # each node is (u,v) where 0<=u,v<=99
        # we'll initially infect those near the middle
        initial_infections = [(u, v) for (u, v) in G if 45 < u < 55 and 45 < v < 55]

        H = nx.DiGraph()  # the spontaneous transitions
        H.add_edge('Sus', 'Vac', rate=0.01)
        H.add_edge('Inf', 'Rec', rate=1.0)

        J = nx.DiGraph()  # the induced transitions
        J.add_edge(('Inf', 'Sus'), ('Inf', 'Inf'), rate=2.0)

        IC = defaultdict(lambda: 'Sus')  # initial condition
        for node in initial_infections:
            IC[node] = 'Inf'

        return_statuses = ['Sus', 'Inf', 'Rec', 'Vac']

        color_dict = {'Sus': '#009a80', 'Inf': '#ff2000', 'Rec': 'gray', 'Vac': '#5AB3E6'}
        pos = {node: node for node in G}
        tex = False
        sim_kwargs = {'color_dict': color_dict, 'pos': pos, 'tex': tex}

        sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax=30,
                                             return_full_data=True, sim_kwargs=sim_kwargs)

        times, D = sim.summary()
        #
        # imes is a numpy array of times.  D is a dict, whose keys are the entries in
        # return_statuses.  The values are numpy arrays giving the number in that
        # status at the corresponding time.

        newD = {'Sus+Vac': D['Sus'] + D['Vac'], 'Inf+Rec': D['Inf'] + D['Rec']}
        #
        # newD is a new dict giving number not yet infected or the number ever infected
        # Let's add this timeseries to the simulation.
        #
        new_timeseries = (times, newD)
        sim.add_timeseries(new_timeseries, label='Simulation',
                           color_dict={'Sus+Vac': '#E69A00', 'Inf+Rec': '#CD9AB3'})

        sim.display(time=6, node_size=4, ts_plots=[['Inf'], ['Sus+Vac', 'Inf+Rec']])
        plt.savefig('test_Animation_Dynamics_SIR_With_Vaccination_In_Lattice')

        ani = sim.animate(ts_plots=[['Inf'], ['Sus+Vac', 'Inf+Rec']], node_size=4)
        ani.save('test_Animation_Dynamics_SIR_With_Vaccination_In_Lattice.mp4', fps=5, extra_args=['-vcodec', 'libx264'])



