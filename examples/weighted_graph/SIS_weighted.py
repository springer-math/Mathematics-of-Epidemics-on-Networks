import networkx as nx
import EoN
import matplotlib.pyplot as plt


N= 100000
rho = 0.001
gamma = 1

G = nx.configuration_model([2,6]*int(N/2)) #N nodes, half have degree 6 and half degree 2

#assign edge weights to be product of degree.  Also give another weight to be inverse of product of degrees
weight_sum = 0
inv_weight_sum = 0

for edge in G.edges():
    G.edge[edge[0]][edge[1]]['weight'] = G.degree(edge[0])*G.degree(edge[1])
    G.edge[edge[0]][edge[1]]['inv_weight'] = 1./(G.degree(edge[0])*G.degree(edge[1]))
    weight_sum += G.degree(edge[0])*G.degree(edge[1])
    inv_weight_sum += 1./(G.degree(edge[0])*G.degree(edge[1]))

#first do it with weight, scaled so that average weight is 1.
t, S, I = EoN.fast_SIS(G, G.number_of_edges()/weight_sum, gamma, rho = rho, transmission_weight= 'weight', tmax = 5)
plt.plot(t, I, label = 'weight')


t, S, I = EoN.fast_SIS(G, G.number_of_edges()/inv_weight_sum, gamma, rho = rho, transmission_weight= 'inv_weight', tmax = 5)
plt.plot(t, I, label = 'inv_weight')

plt.legend(loc = 'lower right')
plt.savefig('SIS_weighted.pdf')