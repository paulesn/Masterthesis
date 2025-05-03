## k-(hop)-spanner

In graph theory, a **hop spanner** (also known simply as a **spanner**) is a subgraph that approximates distances in the original graph. Specifically, a **hop spanner** is defined in terms of the number of edges (hops) used to connect nodes.
### **Formal Definition:**

Given a graph $G = (V, E)$ and a positive integer parameter $k$, a $k$-hop spanner (or simply $k$-spanner) is a subgraph $H = (V, E{\prime})$, with $E{\prime} \subseteq E$, such that for every pair of vertices $u, v \in V$:
$$
\text{dist}_H(u,v) \leq k \cdot \text{dist}_G(u,v)
$$
Here:
- $\text{dist}_G(u,v)$ denotes the shortest path length (in hops or edges) between nodes $u$ and $v$ in graph $G$.
- $\text{dist}_H(u,v)$ denotes the shortest path length between the same nodes $u$ and $v$ in the subgraph $H$.
### **Key Points:**

- **Stretch Factor** $k$: A hop spanner is characterized by its **stretch factor**, $k$, representing how much longer paths can become in the subgraph compared to the original graph. A lower $k$ indicates a closer approximation to the original graph’s distances.
- **Purpose**: Spanners reduce complexity while maintaining approximate distances, thus enabling faster and simpler computation of shortest paths, efficient routing, and network design.
- **Applications**: Network routing, distributed systems, computational geometry, approximation algorithms, and optimizing large-scale graph operations.

## Threshold graph

Given a **complete graph** $G=(V,E)$ where each edge $e\in E$ has an associated **weight** $w(e)$, you construct a hop spanner based on a threshold $j$ by:

- **Step 1**: Start with the original complete graph $G$.
- **Step 2**: **Remove all edges** from the graph whose weight is greater than $j$. Formally, the resulting graph $G_j = (V, E_j)$ is defined by:
  $$
    E_j = \{e \in E \mid w(e) \leq j\}
    $$
The resulting subgraph $G_j$ is called a **threshold graph** or **threshold-based subgraph**.