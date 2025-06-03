
## Problem statement
Given a set of Points in $\mathbb R^2$ we have a complete Graph $G = (P, P \times P)$ with edge weights
$$
\omega(p,q) = |pq|^\delta  C_p
$$
 - $\delta$ is a constant
 - $|pq|$ is the euclidian distance
 - $C_p$ is a node dependent offset.

A path $\pi$ from any node $p$ to any other node $q$ is considered the $t$-approximate shortest path, if $\pi$ is at most $t$ times longer than any other path $\pi´$ from $p$ to $q$.



> [!note] 
> The paper approaches the problem from another direction. While I want to minimize $k$ (the number of hops) with a fixed maximum edge weight, the paper has a fixed $k$ and tries to identify a route with a minimal maximum edge weight

