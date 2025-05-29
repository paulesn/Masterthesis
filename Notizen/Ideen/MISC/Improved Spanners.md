> [!note] 
> The base idea of this approach is to create a spanner (maybe based on WSPD) and replace all edges longer than the jump distance with the actual path that only uses edges within the jump distance.
> 
> - The idea assumes a fixed jump range
> - Maybe the idea can use jump amount as distance mentric
> - finally design X spanner for X distances and always use the one with the largest jump distance smaller than the actual jump distance

- [ ] #Idee: könnte man das nicht gestacked machen? also dass man mehrere Edge list hat die man über einander legt. das würde storage reduzieren

- A WSPD-based t-spanner has $\mathcal O(|V|)$ edges with $s = 4\frac{t+1}{t-1}$.
- [ ] How many of those edges are longer than the jump range $j$?


> [!note] Hop metric
> The WSPD based spanner work for any metric:
> A Metric is:
> A metric space is a pair $(S, \delta)$, where $S$ is a (finite or infinite) set, whose elements are called points, and $\delta: S \times S \to R$ is a function that assigns a distance $\delta(p, q)$ to any two points $p$ and $q$ in $S$, and that satisfies the following three conditions: 
> 1. For all points $p$ and $q$ in $S, \delta(p, q) \geq 0$. 
> 2. For all points $p$ and $q$ in $S, \delta(p, q) = 0$ if and only if $p = q$. 
> 3. For all points $p$ and $q$ in $S, \delta(p, q) = \delta(q, p)$. 
> 4. For all points $p, q$, and $r$ in $S, \delta(p, q) \leq \delta(p, r) + \delta(r, q)$.
> 
> We want to show that the smallest number of hops is a metric

1. For all points $p$ and $q$ in $S, \delta(p, q) \geq 0$
   There are no negative numbers of hops.
2. For all points $p$ and $q$ in $S, \delta(p, q) = 0$ if and only if $p = q$
   not hopping ($p=q$) results in $0$ hops.
   everything else needs at most $1$ hop.
3. For all points $p$ and $q$ in $S, \delta(p, q) = \delta(q, p)$
   As long as we ignore neutron stars, this applies
4. For all points $p, q$, and $r$ in $S, \delta(p, q) \leq \delta(p, r) + \delta(r, q)$
   Assuming there is a way to jump from $p$ to $q$ using $r$ which is shorter than the shortest way $p \to q$, then it is not the shortest way and the way over $r$ becomes the shortest way and the condition holds.


