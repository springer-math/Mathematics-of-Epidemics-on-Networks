import networkx as nx
from EoN import *
import matplotlib.pyplot as plt
import scipy
import random

colors = ['#5AB3E6', '#FF2000', '#009A80', '#E69A00', '#CD9AB3', '#0073B3', '#F0E442']


def test_discrete_SIR():
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


def test_basic_discrete_SIR():
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
    filename = 'basic_discrete_SIR_test.pdf'
    plt.savefig(filename)
    print("check {} for good match".format(filename))


def test_estimate_SIR_prob_size():
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


def test_SIR_dynamics():
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
    plt.savefig('SIR_dynamics.pdf')


def test_SIS_dynamics():
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


def test_SIS_simulations():
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
        plt.savefig('SIS_sims.pdf')


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
    plt.savefig('test_SIR_final_sizes.pdf')

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
#     plt.savefig('SIR_individual_based.pdf')

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
#     plt.savefig('SIS_individual_based.pdf')


def test_pair_based():
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
    plt.savefig('SIS_pair_based.pdf')

    print('testing SIR_pair_based')
    t, S, I, R = EoN.SIR_pair_based(G, 2, 0.5, tmax=5, tcount=101)
    print('creating SIR fig')
    plt.clf()
    plt.plot(t, S)
    plt.plot(t, I)
    plt.plot(t, R)
    plt.savefig('SIR_pair_based.pdf')

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
#     plt.savefig('SIR_pair_based2.pdf')
#     print('done with pair_based2.pdf')


def test_SIS_homogeneous_meanfield():
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
    plt.savefig('SIS_homogeneous_meanfield.pdf')


def test_SIR_homogeneous_meanfield():
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
    plt.savefig('SIR_homogeneous_meanfield.pdf')


def test_SIS_homogeneous_pairwise():
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
    plt.savefig('SIS_homogeneous_pairwise.pdf')


def test_SIR_homogeneous_pairwise():
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
    plt.savefig('SIR_homogeneous_pairwise.pdf')


def test_SIS_heterogeneous_meanfield():
    print("testing SIS_heterogeneous_meanfield")
    Sk0 = scipy.arange(100) * 100
    Ik0 = scipy.arange(100)

    t, S, I = EoN.SIS_heterogeneous_meanfield(Sk0, Ik0, 1, 10, tmax=1)
    print("plotting SIS_heterogeneous_meanfield")
    plt.clf()
    plt.plot(t, S)
    plt.plot(t, I)
    plt.savefig('SIS_heterogeneous_meanfield.pdf')


def test_SIR_heterogeneous_meanfield():
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
    plt.savefig('SIR_heterogeneous_meanfield.pdf')


def test_SIS_heterogeneous_pairwise():
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
    plt.savefig('SIS_heterogeneous_pairwise.pdf')


def test_SIR_heterogeneous_pairwise():
    print("SIR_heterogeneous_pairwise not yet tested")


def test_SIS_compact_pairwise():
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
    plt.savefig('SIS_compact_pairwise.pdf')


def test_SIR_compact_pairwise():
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
    plt.savefig('SIR_compact_pairwise.pdf')


def test_SIS_super_compact_pairwise():
    print("SIS_super_compact_pairwise not yet tested")


def test_SIR_super_compact_pairwise():
    print("SIR_super_compact_pairwise not yet tested")


def test_SIS_effective_degree():
    print("SIS_effective_degree not yet tested")


def test_SIR_effective_degree():
    print("SIR_effective_degree not yet tested")


def test_SIS_compact_effective_degree():
    print("SIS_compact_effective_degree not yet tested")


def test_SIR_compact_effective_degree():
    print("SIR_compact_effective_degree not yet tested")


# test_basic_discrete_SIR()
#
# '''tests of ODE models -----'''
# # test_my_odeint()
# # test_SIR_individual_based()
# # test_SIS_individual_based()
# test_pair_based()
# test_SIS_homogeneous_meanfield()
# test_SIR_homogeneous_meanfield()
# test_SIS_homogeneous_pairwise()
# test_SIR_homogeneous_pairwise()
# test_SIS_heterogeneous_meanfield()
# test_SIR_heterogeneous_meanfield()
# test_SIS_heterogeneous_pairwise()
# test_SIR_heterogeneous_pairwise()
# test_SIS_compact_pairwise()
# test_SIR_compact_pairwise()
# test_SIS_super_compact_pairwise()
# test_SIR_super_compact_pairwise()
# test_SIS_effective_degree()
# test_SIR_effective_degree()
# test_SIS_compact_effective_degree()
# test_SIR_compact_effective_degree()
#
# '''tests of simulation models -----'''
# test_discrete_SIR()
# test_basic_discrete_SIR()
# test_estimate_SIR_prob_size()
# test_SIR_dynamics()
# test_SIS_dynamics()
# test_SIS_simulations()
# test_SIR_final_sizes()

# G= nx.fast_gnp_random_graph(10000,0.001)

# t, S, I, R = EoN.basic_discrete_SIR(G, 0.2)
# print(t, S, I, R)

# t,S,I,R, infection_times, recovery_times = EoN.basic_discrete_SIR(G,0.2,return_full_data=True)
# print(t, S, I, R, infection_times, recovery_times)

# t, S, I, R, infection_times, recovery_times = EoN.basic_discrete_SIR(G,0.2,initial_infecteds=range(4), return_full_data=True)
# print(t, S, I, R, infection_times, recovery_times)



# t, S, I, R, infection_times, recovery_times = EoN.fast_SIR(G, 1, 5, return_full_data = True)
# print(t, S, I, R, infection_times, recovery_times)



# initial_infecteds = range(15)
# t, S, I, Itimes, Rtimes = EoN.fast_SIS(G, 1, 5, initial_infecteds = initial_infecteds, tmax = 10, return_full_data = True)
# print(t, S, I)
# py.plot(t,I)
# py.savefig('fast_SIS_test.pdf')

# for node in initial_infecteds:
#    print(node)
#    print(Itimes[node])
#    print(Rtimes[node])
#    print('\n')

# first_Recs = [Rtimes[node][0] for node in initial_infecteds]
# print(first_Recs, sum(first_Recs)*1./len(initial_infecteds))


# print('\n\n\nStill need to fix:')
# print('SIS_pair_based')
# print('what is up with SIR_pair_based2?')
# print('alternate_dSIR_pair_based2?')
# print('get_rate_functions called in SIR_pair_based2 - can we get rid of all of this?')
#
# print('\n')
#
# print('homogeneous_mean_field from graphs...')
#
# print('attack_rate_non_Markovian')
