import EoN
import networkx as nx
import matplotlib.pyplot as plt
import scipy


N=1000
gamma = 1.
iterations = 200
rho = 0.05
tmax = 10
tcount = 1001

report_times = scipy.linspace(0,tmax,tcount)

deg_dist1 = [18,22]*(N/2)
deg_dist2 = [5,35]*(N/2)
ax1 = plt.gca()#axes([0.1,0.1,0.9,0.9])
ax2 = plt.axes([0.44,0.2,0.4,0.4])
for deg_dist, ax in zip([deg_dist1, deg_dist2], [ax1, ax2]):
    kave = sum(deg_dist1)*1./N
    tau = 2*gamma/kave
    Isum = scipy.zeros(tcount)

    for counter in range(iterations):
        G = nx.configuration_model(deg_dist)
        t, S, I = EoN.fast_SIS(G, tau, gamma, tmax=tmax, rho=rho)
        I = I*1./N
        I = EoN.subsample(report_times, t, I)
        Isum += I
    ax.plot(report_times, Isum/iterations, color='grey', linewidth=5, alpha=0.3)
    
    S0 = (1-rho)*N
    I0 = rho*N
    
    t, S, I = EoN.SIS_homogeneous_meanfield(S0, I0, kave, tau, gamma, tmin=0, tmax=tmax, 
                                tcount=tcount)
    ax.plot(t, I/N, '--')

    SI0 = (1-rho)*N*kave*rho
    SS0 = (1-rho)*N*kave*(1-rho)
    t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, kave, tau, gamma, tmin = 0, 
                                tmax=tmax, tcount=tcount)
    ax.plot(t, I/N)

plt.savefig('fig4p10.pdf')