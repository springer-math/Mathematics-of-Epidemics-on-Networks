import networkx as nx
import EoN
import matplotlib.pyplot as plt


r'''This code demonstrates a change in parameters at time t1.  To do this, it
will run to time t1 and then stop, and restart with the new parameters.

To demonstrate that stopping and restarting doesn't do anything weird, we first
run to time t0, stop, and then restart with the original parameters.  The resulting
simulation should just look like nothing changes.

Another way to do this may be through writing a custom transmission rate
function and using the nonMarkov version of this function.
'''

def get_affected_nodes_at_end(infection_time, recovery_time):
    recovered = set(recovery_time.keys())
    infected = set(node for node in infection_time if node not in recovered)
    return infected, recovered
t0= 2
t1 = 4
tmax = 8

gamma = 1
tau0 = 1
tau1 = 0.5

N= 100000
kave = 4
rho = 0.001

G = nx.fast_gnp_random_graph(N, kave/(N-1.))

times0, S0, I0, R0, infection_time, recovery_time = EoN.fast_SIR(G, tau0, gamma, rho = rho, tmax = t0, return_full_data=True)

infected, recovered = get_affected_nodes_at_end(infection_time, recovery_time)


times1, S1, I1, R1, infection_time, recovery_time = EoN.fast_SIR(G, tau0, gamma, initial_infecteds = infected, initial_recovereds = recovered, tmin = t0, tmax = t1, return_full_data=True)

infected, recovered = get_affected_nodes_at_end(infection_time, recovery_time)


times2, S2, I2, R2 = EoN.fast_SIR(G, tau1, gamma, initial_infecteds = infected, initial_recovereds=recovered, tmin=t1, tmax=tmax)


plt.plot(times0, I0)
plt.plot(times1, I1)#the first two have the same parameters, so the transition should be as if it were a single simulation
plt.plot(times2, I2)#the infectiousness reduced, so a sharp change should be visible
plt.savefig('SIR_change_tau.pdf')