
> [!tip] 
> Given a Set of points on one dimension and a jump distance of $j$ and the start end end points at a distance $d$, the worst case needs $2\frac{d}{j}$ jumps.

Proof
We space $\frac{d}{j}$ points between $s$ and $t$ so that the pairwise distance between all points is greater than $d$. We call these points $p_i$ In this scenario, there is no valid path from $s$ to $t$. 
To add an valid path, we place another points $m_i$ at a fixed distance $d_2<d$ next to each point $p_i$.
Then we have a path $s, p_1,m_1, p_2, m_2, ..., p_{\lfloor\frac{d}{j}\rfloor},m_{\lfloor\frac{d}{j}\rfloor}$
Adding any point will not increase the length of the path and removing any point will break the path. 