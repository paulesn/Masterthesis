'''

This class is meant to load the fmi graph files
'''

import networkx as nx
from tqdm import tqdm
import random

def load(path):
    '''
    Load the graph from the given path
    '''
    with open(path, 'r') as f:
        node_count = -1
        edge_count = -1
        count = 0
        fail_count = 0
        G = nx.DiGraph()
        with tqdm(f) as pbar:
            for line in pbar:
                if str(line).startswith('#'):
                    continue # skip comments
                elif line == '\n':
                    continue
                elif ' ' not in line:
                    if node_count == -1:
                        node_count = int(line)
                        print(f"Node count: {node_count}")
                    elif edge_count == -1:
                        edge_count = int(line)
                        pbar.total = edge_count+node_count+10
                        print(f"Edge count: {edge_count}")
                    else:
                        print("Error: too many lines")
                        return
                else:
                    if count < node_count:
                        count += 1
                        node = line.split(' ')
                        G.add_node(int(node[0]), X=float(node[2]), Y=float(node[3]), index=int(node[1]))
                    else:
                        edge = line.split(' ')
                        try:
                            u = int(edge[0])
                        except:
                            print(f"Error: failed to parse edge ({edge[0]}, {edge[1]})")
                            fail_count += 1
                            continue
                        v = int(edge[1])
                        if u not in G or v not in G:
                            # print(f"Error: edge ({edge[0]}, {edge[1]}) not found")
                            fail_count += 1
                            continue
                        G.add_edge(int(u), int(v), weight=float(edge[2]))
    print(f"Failed to add {fail_count} edges")
    return G


def generate_graph_cutout(source_node, size, graph, random_state=42):
    '''
    Generate a cutout of the graph by using a source node and randomly adding edges
    '''
    if graph.number_of_nodes() < size:
        raise("Graph is toosmall for the requested cutout size")
    upper_limit = 2*size
    cutout = set([source_node])
    frontier = set([source_node])
    for n in graph.neighbors(source_node):
        frontier.add(n)
    random.seed(random_state)
    with tqdm(total=size) as pbar:
        while len(cutout) < size:
            x = random.choice(list(frontier))
            frontier.remove(x)
            for n in graph.neighbors(x):
                frontier.add(n)
            frontier = frontier - cutout
            cutout.add(x)
            pbar.update(1)
    return graph.subgraph(list(cutout))


def generate_fmi_data_from_graph(graph):
    '''
    returns a string in the fmi format
    :param graph:
    :return:
    '''
    fmi_data = ""
    fmi_data += f"{graph.number_of_nodes()}\n"
    fmi_data += f"{graph.number_of_edges()}\n"
    for node in graph.nodes():
        fmi_data += f"{node} {graph.nodes[node]['index']} {graph.nodes[node]['X']} {graph.nodes[node]['Y']}\n"
    for edge in graph.edges():
        fmi_data += f"{edge[0]} {edge[1]} {graph.edges[edge]['weight']}\n"
    return fmi_data