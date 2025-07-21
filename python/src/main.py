import fmi_loader
import random
import networkx as nx

graph = fmi_loader.load("../../data/theta_spanner.fmi")

for i in range(5000):
    source = random.choice(range(graph.number_of_nodes()-1))
    target = random.choice(range(graph.number_of_nodes()-1))
    if source == target:
        continue
    # l = nx.shortest_path_length(graph, source, target)
    # if l == float('inf'):
        # print(f"Shortest path from {source} to {target}: {l}")

#  get number of connected components
n = nx.number_connected_components(graph.to_undirected())
print(f"Number of connected components: {n}")