import EoN
import networkx as nx
import matplotlib.pyplot as plt
import scipy
import random

        

def sim_and_plot(G, tau, gamma, rho, tmax, tcount, ax):
    t, S, I, R= EoN.fast_SIR(G, tau, gamma, rho = rho, tmax = tmax)
    report_times = scipy.linspace(0, tmax, tcount)
    I = EoN.subsample(report_times, t, I)
    ax.plot(report_times, I/N, color='grey', linewidth=5, alpha=0.3)
    
    t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, 
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N, '--')    
    t, S, I, R = EoN.SIR_compact_pairwise_from_graph(G, tau, gamma, rho=rho,
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N)
 
    t, S, I, R = EoN.SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, 
                                                    tmax=tmax, tcount=tcount)
    ax.plot(t, I/N, '-.')

N=50000
gamma = 1
rho = 0.05
tmax = 15
tcount = 1001

deg_seq = [30]*int(N/2) + [70]*int(N/2)
G = nx.configuration_model(deg_seq)
kave = sum(deg_seq)/N
ksqave = sum(k*k for k in deg_seq)/N

tau= 2*gamma*kave/ksqave

fig = plt.figure(1)
main = plt.axes()
sim_and_plot(G, tau, gamma, rho, tmax, tcount, main)
plt.ylabel('Prevalence')


deg_seq = [5]*int(N/2) + [15]*int(N/2)
G = nx.configuration_model(deg_seq)

kave = (sum(deg_seq)/N)
ksqave = sum(k*k for k in deg_seq)/N

tau= 2*gamma*kave/ksqave

fig = plt.figure(1)
inset = plt.axes([0.5,0.45,0.4,0.4])
sim_and_plot(G, tau, gamma, rho, tmax, tcount, inset)


plt.savefig('fig5p5.pdf')