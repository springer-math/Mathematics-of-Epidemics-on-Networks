import networkx as nx
import EoN
import matplotlib.pyplot as plt


r'''This code demonstrates a change in parameters at time t1.  To do this, it
will run to time t1 and then stop, and restart with the new parameters.

To demonstrate that stopping and restarting doesn't do anything weird, we first
run to time t0, stop, and then restart with the original parameters.  The resulting
simulation should just look like nothing changes.

Once the nonMarkovian SIS function exists, another way to do this may be
through writing a custom transmission rate
function 
'''

    
    


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


times0, S0, I0, infection_time, recovery_time = EoN.fast_SIS(G, tau0, gamma, rho = rho, tmax = t0, return_full_data=True)


#the infected nodes are those that either were infected and never recovered 
#(which won't be in recovery_time)
#or have been infected multiple times and their last infection was after
#their last recovery

infected = set(node for node in infection_time if node not in recovery_time or infection_time[node][-1]> recovery_time[node][-1])
times1, S1, I1, infection_time, recovery_time = EoN.fast_SIS(G, tau0, gamma, initial_infecteds= infected, tmin=t0, tmax=t1, return_full_data=True)

infected = set(node for node in infection_time if node not in recovery_time or infection_time[node][-1]> recovery_time[node][-1])
times2, S2, I2, infection_time, recovery_time = EoN.fast_SIS(G, tau1, gamma, initial_infecteds= infected, tmin=t1, tmax=tmax, return_full_data=True)

plt.plot(times0, I0)
plt.plot(times1, I1)
plt.plot(times2, I2)

plt.savefig('SIS_change_tau.pdf')
