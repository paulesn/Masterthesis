import fmi_loader
import random
import networkx as nx

graph = fmi_loader.load("../../data/test.fmi")
graph.to_undirected()

for i in range(1):
    source = 3
    target = 1
    if source == target:
        continue
    l = nx.shortest_path_length(graph, source, target, weight='weight')
    print(l)
    # if l == float('inf'):
        # print(f"Shortest path from {source} to {target}: {l}")

#  get number of connected components
n = nx.number_connected_components(graph.to_undirected())
print(f"Number of connected components: {n}")