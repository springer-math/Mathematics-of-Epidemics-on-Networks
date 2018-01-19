import networkx as nx
import EoN
import matplotlib.pyplot as plt

r'''Each time the infected population hits 40000, we'll introduce an intervention.
To do this, we run it for short intervals until we overshoot, but then backtrack to the time at which
40000 is hit, and then run again with intervention in place.'''


