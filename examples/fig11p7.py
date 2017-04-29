import EoN
import networkx as nx
import matplotlib.pyplot as plt
import scipy

from collections import Counter

#
#
# This code generates a Watts-Newman graph and then finds its degree
# distribution.  Using this degree distribution it predicts the epidemic
# dynamics.  
#

N=10**4  #population size
k=10
tau = 0.5 #transmission rate
gamma = 1. #recovery rate
threshold=100 #outbreaks larger than 100 are considered an epidemic.

def get_deg_probs(G):
    r'''Gets the proportion of nodes in the graph with each observed
    degree'''
    degree_count = Counter(G.degree().values())
    degree_prob = {x:degree_count[x]/float(G.order()) for x in degree_count}
    return degree_prob

def plot_deg_dist(fig,degree_prob):
    deg, prob = zip(*sorted(degree_prob.items()))
    fig.plot(deg, prob, 'x', label = 'p={}'.format(p))

def get_xmax(t,L):
    for index in range(len(L)):
        if L[-index]>1:
            break
    print L[-index]
    #print L[-index:]
    return t[-index]
 
def SIR_process(G, degree_prob, tau, gamma, tmax = 10):
    N = G.order()
    plt.figure(2)
    plt.clf()
    plt.figure(3)
    plt.clf()
    plt.figure(5)
    plt.clf()
    for index, starting_node in enumerate([x*N/100. for x in range(100)]):
        plt.figure(2)
        t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds = [starting_node])
        subt = scipy.linspace(0, t[-1], 101)
        subI, subR = EoN.subsample(subt, t, I, R)
        plt.plot(subt,subI)
        if R[-1]>500:
            plt.figure(3)
            shift = EoN.get_time_shift(t, R, threshold)
            plt.plot(subt-shift, subI)
            plt.figure(5)
            plt.plot(subt-shift, subR)
    #t, S, I, R = EoN.EBCM(degree_prob, tau, gamma, rho)
    rho = 1./N
    def psi(x):
        return sum(degree_prob[k]*x**k for k in degree_prob)
    def psiPrime(x):
        return sum(k*degree_prob[k]*x**(k-1) for k in degree_prob)
    t, S, I, R = EoN.EBCM_uniform_introduction(N, psi, psiPrime, tau, gamma, rho, tmax=tmax)
    shift = EoN.get_time_shift(t, R, threshold)

    plt.figure(2)
    #plt.savefig('sw_SIR_epi_N{}_p{}_k{}_tau{}.pdf'.format(N,p,k,tau))
    plt.figure(3)
    plt.plot(t-shift, I, '--')
    plt.xlabel('$t$', fontsize = 18)
    plt.ylabel('$I$', fontsize = 18)
    #plt.set_xtick_labels(fontsize = 15)
    xmax = get_xmax(t-shift,I)
    plt.axis(xmax=xmax)
    plt.savefig('sw_SIR_epi_N{}_p{}_k{}_tau{}_shifted.pdf'.format(N, p, k, tau))
    plt.figure(5)
    plt.plot(t-shift, R, '--')
    #plt.savefig('sw_SIR_epi_N{}_p{}_k{}_tau{}_shifted_R.pdf'.format(N, p, k, tau))
    
def SIS_process(G, degree_prob, tmax, tau, gamma):
    N=G.order() 
    plt.figure(5)
    plt.clf()
    plt.figure(6)
    plt.clf()
    for index, starting_node in enumerate([x*N/10. for x in range(10)]):
        plt.figure(5)
        t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds = [starting_node], tmax = tmax)
        print I[-1]
        subt = scipy.linspace(0, tmax, 501)
        subI = EoN.subsample(subt, t, I)
        plt.plot(subt,subI)
        if I[-1]>100:
            plt.figure(6)
            shift = EoN.get_time_shift(t, I, 1000)
            plt.plot(subt-shift, subI)
    plt.figure(5)
    plt.savefig('sw_SIS_epi_N{}_p{}_k{}_tau{}.pdf'.format(N, p, k, tau))
    plt.figure(6)
    plt.savefig('sw_SIS_epi_N{}_p{}_k{}_tau{}_shifted.pdf'.format(N, p, k, tau))
    

plt.figure(1)    
deg_dist_plot = plt.gca()
SIR_tmaxes=[20]#, 15, 15, 15]#[40, 35, 25, 25]#[10000,100,100,10,10,10,10,10,10,10,10,10]
ps = [0.01]#, 0.1, 0.5, 1]#[0.2, 0.4, 1, 2]#[0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
for SIR_tmax, p in zip(SIR_tmaxes, ps):
    print p
    if p>1:
        raise MyError('p>1')
    G = nx.watts_strogatz_graph(N, k, p)
    degree_prob = get_deg_probs(G)
    plot_deg_dist(deg_dist_plot, degree_prob)
    print 'SIR'
    SIR_process(G, degree_prob, tau, gamma, tmax = SIR_tmax)
    #print 'SIS'
    #SIS_process(G, degree_prob, SIS_tmax, tau, gamma)
plt.figure(1)
plt.legend()
plt.savefig('sw_deg_dist_N{}_p{}_k{}_tau{}.pdf'.format(N, p, k, tau))

plt.figure(4)

