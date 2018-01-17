import cProfile

import EoN

import networkx as nx

N=1000000
kave = 5.
tau = 1.
gamma = 1.
rho = 0.001

G = nx.fast_gnp_random_graph(N, kave/(N-1.))

def foo():
    t, S, I = EoN.fast_SIS(G, tau, gamma, rho=rho, tmax=5)


print('about to profile')
cProfile.run('foo()', 'stats.dat')

import pstats
p = pstats.Stats('stats.dat')
p.strip_dirs().sort_stats(-1).print_stats()
