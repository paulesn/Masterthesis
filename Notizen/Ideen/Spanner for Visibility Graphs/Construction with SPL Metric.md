Given a visiblity graph $G = (V,E)$, and for each vertice $v \in V$ we know its position in $\mathbb R^2: pos(v)$ 
We construct a Well separated pairs decomposition on V with an $s>4$.
Each time we identify a well separated pair $A,B$, we choose any $r_A \in A, r_B\in B$.
we start a single source ($r_A$) to multiple targets ($\forall v \in A$) shortest path query on $G$. If there is any shortest path that is more than $t$ times as long as the direct euklidian distance it represents, the pair is handled as if it was not well separated.  

and calculate its shortest path $\pi$ in $G$. All Edges in $\pi$ are then added to our WSPD-spanner $S$.


[^1]: This is in $\mathcal O(n^2)$ and therefore potentially to expensive


This could be archieved using: The inner distance of set in a well seperated pair is not the radius of circle encompassing all points it is the diameter of the subgraph of the visibility graph of that point set. This can be calculated in $\mathcal O(|V|+|E|)$

---

## Shortest Path Length is a metric 

> [!tip] Lemma 
> The shortest Path distance between to nodes in a Graph, which edge weights are the euklidian distance between the nodes is a metric.

Proof:

1. For all points $p$ and $q$ in $S, δ(p, q) ≥ 0$. 
   as the edge weights are all positive and added together no shortest path can be negative
   
2. For all points $p$ and $q$ in $S, δ(p, q) = 0$ if and only if $p = q$. 
   The shortest path between a node and itself is always $0$.
   
3. For all points $p$ and $q$ in $S, δ(p, q) = δ(q, p)$. 
   The edges are undirected therefore the shortest paths are identical in both directions
   
4. For all points $p, q$, and $r$ in $S, δ(p, q) ≤ δ(p, r) + δ(r, q)$.
   if there would be a $r$ so that $δ(p, q) > δ(p, r) + δ(r, q)$ then $(p\to r \to q)$ would be a shorter path than the shortest path which is a. Thats a contradiction.

Therefore shortest path length in our visibility graph is a metric and can be used as the basis for a spanner.

**q.e.d.**

We denote $|pq|_{SPD}$ as the shortest path distance between the nodes $p$ and $q$.
## Invalidating Well Separated Pairs 

> [!tip] Lemma 
> During the top-down Quad Tree construction of a well separated pairs decomposition, if we find a well separated pair $A,B$ with more than one node inside $A$ or $B$. The pair can be treated as not well separated and the resulting graph will be a valid spanner with at most the same stretch as if all pairs are treated normally.

Proof:
Let $A$ and $B$ be well separated. Then all child combinations in the Quadtree: $C = \{X,Y| X \text{ child of or equal to } A, Y \text{ child of or equal to } B\} \backslash \{A,B\}$ are well separated, because their radius is at most the same as the radius in the pair $A,B$ and the distance is at most the distance in $A,B$. 
Therefore, we can choose a pair of nodes for each combination in $C$. each of them is a valid replacement for the pair in the original pairing as we could have chosen them in the original pairing.

**q.e.d.**


## Construction 

- Build a WSPD using SPD as a metric.
	- for each pair $A,B$ found, start a multiple sources multiple targets shortest path search with all elements of $A$ as sources and all elements of $B$ as targets. The first path found is added to the spanner. 
		- Maybe this could use A* with the distance to the center of $B$ as a heuristic.  

## Size

> [!tip] Lemma
> A spanner build with the construction detailed above has $\mathcal O(n\mathcal \ell)$ edges.
> With $\ell$ being the average unweighted shortest path length 

Proof:


using the corollary of the book:
> [!PDF|] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=186&selection=263,1,346,7|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.168]]
> > **Corollary 9.4.7 (WSPD-Spanner)**. Let $S$ be a set of $n$ points in $\mathbb R^d$ and let $t > 1$ be a real number. In $\mathcal O(\frac{n \log n + n}{(t − 1)^d} )$ time, we can construct a $t$-spanner for $S$ having $\mathcal O(\frac n{(t − 1)^d} )$ edges. Thus, the WSPD-spanner has $\mathcal O(n)$ edges and can be constructed in $\mathcal O(n \log n)$ time.
> 

Assuming we connect each well separated pair just with any direct edge, then we have, according to corollary 9.4.7, $\mathcal O(n)$ edges. 
In our construction, each edge is replaced by the shortest path connecting the two sets of a pair.
Therefore it would be multiplied with the average unweigthed shortest path length.



### Example Graphs and their size

beispiel finden für wspd graphen die große clusterdurch die visibility co,ponente verhindern dadurch den spanner wachsen lassen
Ziel: zeigen dass es visiility graphen hat bei deinen jede wspd quadratisch ist 

#### Full Graph 
$\ell = 1 \implies \mathcal O (n)$


#### Path Graph 
We have $1$ path of length $n$
we have $2$ paths of length $n-1$
we have $3$ paths of length $n-2$
we have $p$ paths of length $n-(p-1)$

$$
\begin{align}
l &= n-(p-1) \\
l &= n-p-1 \\
l+1 &= n-p \\
l+1-n &= -p \\
-(l+1-n) &= p \\
p &= -l-1+n\\
p &= n-l-1\\
\end{align}
$$

$$
\begin{align}
\ell &= \frac{\sum_{l=2}^n (l(n+1-l))}{\sum_{l=2}^n (n+1-l)}\\
&= \frac{\sum_{l=2}^n (l\cdot n+l-l^2))}{\sum_{l=2}^n (n+1-l)}\\
&= \frac{\sum_{l=2}^n (l\cdot n+l-l^2))}{\sum_{l=2}^n (n+1-l)}\\
&= \frac{n+4}{3} & \text{for } n > 1 \text{ by Wolfram Alpha}\\
\end{align}
$$

$\implies \mathcal O (n^2)$


## Construction Time

According to corollary 9.4.7 it takes $\mathcal O(|V| \log |V|)$ time to construct a WSPD spanner. As we need $\mathcal O(|V|+|E|)$ to identify the diameter of a subgraph and $\mathcal O((|E|+|V|) \log |V|)$ to find the shortest path between two sets

$$
\begin{align}
(|V|\cdot \log |V|)\cdot (|V|+|E| + ((|E|+|V|) \log |V|)) \\
|V|\cdot \log |V|\cdot (|V|+|E| + ((|E|+|V|) \log |V|)) \\
|V|\cdot (\log |V|\cdot |V|+ \log |V|\cdot|E| + \log |V|\cdot((|E|+|V|) \log |V|)) \\
|V|\cdot (\log |V|\cdot |V|+ \log |V|\cdot|E| + \log |V|\cdot(|E|+|V|) \log 
|V|) \\
|V|\cdot (\log |V|\cdot |V|+ \log |V|\cdot|E| + (\log |V|)^2\cdot(|E|+|V|)) \\
|V|\cdot\log |V|\cdot |V|+ |V|\cdot\log |V|\cdot|E| + |V|\cdot(\log |V|)^2\cdot(|E|+|V|)) \\
|V|^2\cdot\log |V|+ |V|\cdot\log |V|\cdot|E| + (\log |V|)^2\cdot|V|\cdot(|E|+|V|)) \\
|V|^2\cdot\log |V|+ |V|\cdot\log |V|\cdot|E| + (\log |V|)^2\cdot(|V|\cdot|E|+|V|^2)) \\
\end{align}
$$