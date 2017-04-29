import EoN
import networkx as nx
import matplotlib.pyplot as plt
import scipy
import random

print r"Warning, book says \tau=2\gamma/<K>, but it's really 1.5\gamma/<K>"
print r"Warning - for the power law graph the text says k_{max}=110, but I believe it is 118."

N=1000
gamma = 1.
iterations = 200
rho = 0.05
tmax = 15
tcount = 101

kave = 20


tau = 1.5*gamma/kave


def simulate_process(graph_function, iterations, tmax, tcount, rho, kave, tau, gamma, symbol):
    Isum = scipy.zeros(tcount)
    report_times = scipy.linspace(0,tmax,tcount)
    for counter in range(iterations):
        G = graph_function()
        t, S, I = EoN.fast_SIS(G, tau, gamma, rho=rho, tmax=tmax)
        I = EoN.subsample(report_times, t, I)
        Isum += I
    plt.plot(report_times, Isum*1./(N*iterations), symbol)

#regular
symbol = 'o'
graph_function = lambda : nx.configuration_model(N*[kave])
simulate_process(graph_function, iterations, tmax, tcount, rho, kave, tau, gamma, symbol)

#bimodal
symbol='x'
graph_function = lambda: nx.configuration_model((N/2)*[5,35])
simulate_process(graph_function, iterations, tmax, tcount, rho, kave, tau, gamma, symbol)

#erdos-renyi
symbol = 's'
graph_function = lambda : nx.fast_gnp_random_graph(N, kave/(N-1.))
simulate_process(graph_function, iterations, tmax, tcount, rho, kave, tau, gamma, symbol)

symbol = 'd'
pl_kmax = 118
pl_kmin = 7
pl_alpha = 2.
Pk={}
for k in range(pl_kmin, pl_kmax+1):
    Pk[k] = k**(-pl_alpha)
valsum = sum(Pk.values())
for k in Pk.keys():
    Pk[k] /= valsum
    
#print sum(k*Pk[k] for k in Pk.keys())
def generate_sequence(Pk, N):
    while True:
        sequence = []
        for counter in range(N):
            r = random.random()
            for k in Pk.keys():
                if r< Pk[k]:
                    break
                else:
                    r-=Pk[k]
            sequence.append(k)
        if sum(sequence)%2==0:
            break        
    return sequence

graph_function = lambda : nx.configuration_model(generate_sequence(Pk,N))
simulate_process(graph_function, iterations, tmax, tcount, rho, kave, tau, gamma, symbol)

symbol = '--'
S0 = (1-rho)*N
I0 = rho*N
t, S, I = EoN.SIS_homogeneous_meanfield(S0, I0, kave, tau, gamma, tmax=tmax, tcount=tcount)
plt.plot(t, I/N, symbol)

symbol = '-'
S0 = (1-rho)*N
I0 = rho*N
SI0 = (1-rho)*N*kave*rho
SS0 = (1-rho)*N*kave*(1-rho)
t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, kave, tau, gamma, tmax=tmax, tcount=tcount)
plt.plot(t, I/N, symbol)

plt.savefig('fig4p11.pdf')