>[!info] 
>Given a start $s$ and a target $t$, we draw an edge from $s$ to $t$.
>(1) We choose the longest edge $(u,v)$in our path, identify the point $p$ closest to its midpoint and replace the edge with the new edges $(u,m)$ and $(m,v)$. 
>Repeat (1) until all edges are shorter than the jump range

 - [ ] Stop condition iff there is no path in jump range
 - [ ] Speedup after path is found by checking if there are shortcuts on the path
Does not work