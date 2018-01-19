import EoN
import networkx as nx
import matplotlib.pyplot as plt
import scipy



def sim_and_plot(G, tau, gamma, rho, tmax, tcount, ax):
    t, S, I = EoN.fast_SIS(G, tau, gamma, rho=rho, tmax = tmax)
    report_times = scipy.linspace(0, tmax, tcount)
    I = EoN.subsample(report_times, t, I)
    ax.plot(report_times, I/N, color='grey', linewidth=5, alpha=0.3)
    
    t, S, I, = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, 
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N, '--')    
    t, S, I = EoN.SIS_compact_pairwise_from_graph(G, tau, gamma, rho=rho,
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N)
 
    t, S, I = EoN.SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, 
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N, '-.')


tcount = 1001
tmax = 10.
gamma = 1.
N=10000
rho = 0.05

fig = plt.figure(1)
main = plt.axes()


target_kave = 50.
G = nx.fast_gnp_random_graph(N, target_kave/(N-1.))
kave = sum(G.degree(node) for node in G.nodes())
ksqave =sum(G.degree(node)**2 for node in G.nodes())
tau_c = gamma*kave/ksqave
tau = 2*tau_c

sim_and_plot(G, tau, gamma, rho, tmax, tcount, main)


inset = plt.axes([0.45,0.175,0.45,0.45])
target_kave = 10.

G = nx.fast_gnp_random_graph(N, target_kave/(N-1.))
kave = sum(G.degree(node) for node in G.nodes())
ksqave =sum(G.degree(node)**2 for node in G.nodes())
tau_c = gamma*kave/ksqave
tau = 2*tau_c

sim_and_plot(G, tau, gamma, rho, tmax, tcount, inset)

main.set_xlabel('$t$')
main.set_ylabel('Prevalence')
plt.savefig('fig5p3.png')