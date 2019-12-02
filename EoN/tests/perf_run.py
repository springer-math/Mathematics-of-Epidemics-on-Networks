#!/usr/bin/env python3
import sys
import os
sys.path.append( os.path.dirname(os.path.realpath(__file__)) + "/../..")

import pyperf

# G = nx.grid_2d_graph(50, 50)  # each node is (u,v) where 0<=u,v<=99
# initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]

runner = pyperf.Runner()

runner.timeit(name="Gillespie_SIS", stmt="EoN.Gillespie_SIS(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=False, tmax=10)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")

runner.timeit(name="Gillespie_SIR", stmt="EoN.Gillespie_SIR(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=False, tmax=10)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")

runner.timeit(name="Fast_SIS", stmt="EoN.fast_SIS(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")

runner.timeit(name="Fast_SIR", stmt="EoN.fast_SIR(G, 1.0, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")

runner.timeit(name="basic_discrete_SIS", stmt="EoN.basic_discrete_SIS(G, 1.0, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")

runner.timeit(name="basic_discrete_SIR", stmt="EoN.basic_discrete_SIR(G, 1.0, initial_infecteds=initial_infections, rho = None, tmin = 0, tmax = 10, return_full_data = True)",
                       setup="import EoN; import numpy as np; import networkx as nx; G = nx.grid_2d_graph(50, 50); initial_infections = [(u, v) for (u, v) in G if 23 < u < 27 and 23 < v < 27]")






