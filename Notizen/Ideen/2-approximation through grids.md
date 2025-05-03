
> [!important] Lemma
> Given a set of points in 2D and a jump distance $j$, we can embed the points in an uniform grid with a cell size of $\frac{j}{2\sqrt 2}$. 
> Each cell with at least one point inside of it becomes a vertice and two vertices $A$ and $B$ are connected by an edge if
> - the represented cells are adjacent.
> - the cell between them is empty and there is at least one pair of points $a\in A$ and $b \in B$ with $||ab||\leq j$ 
> 
> A shortest path found in this graph has at most 3 times the amount of hops than the optimal hop path

Proof:
given an optimal path of jumps $\pi = (p_0, p_1, ..., p_n)$.
For each $p_i$ the next point $p_{i+1}$ has to be closer than $j$.
For this proof we define the point $p_i$ is in $P$, the cell adjacent to $P$ is $A$, the cell adjacent to $A$ is $B$ and the cell adjacent to $B$ is $C$:
![[Pasted image 20250430104352.png]]
Therefore, w.A.d.A. we have four cases:
1. $p_{i+1} \in P$
   Then the optimal path needs more jumps than the path found by the algorithm
2. $p_{i+1} \in A$
   As $A$ is an adjacent cell to $P$ we add the edge $(P,A,1)$ to the search graph 
3. $p_{i+1} \in B$
   This can have two subcases
	1. $A$ is empty, then because $p_{i+1} \in B$, there is a point in $b \in B$ and $p \in P$ so that $||pb||\leq j$ and we add the edge $(P,B,1)$ to the search graph.
	2. A is non-empty. Then $(P,A,1)$ is part of the search graph and we repeat case 2 from $A$. In total, this gives us $2$ hops for the $1$ hop jump  $p_{i}\to p_{i+1}$.
4. $p_{i+1}\in C$
   This can have three subcases
	1. Neither $A$ nor $B$ are empty. Then we have the edges $(P,A,1)$ and $(P,B,1)$ in the search graph and we repeat case 2 from $B$
	2. $A$ is non-empty and $B$ is empty. Then there is an edge $(P,A,1)$ and we repeat case 3 from $A$, resulting in up to $3$ hops for the $1$ hop jump  $p_{i}\to p_{i+1}$.
	3. $A$ is empty, then because $p_{i+1} \in C$, there is a point in $c \in C$ and $p \in P$ so that $||pc||\leq j$ and we add the edge $(P,C,1)$ to the search graph.

In the wort case, each hop in the optimal path is case 4, resulting in 3 hops found by the algorithm for each optimal hop.
q.e.d.

With a cell size of $c = \frac{j}{2\sqrt 2}$ we can use the pythagoras to identify this distance:
$$
\begin{align}
(0.5d)^2 &= \left(\frac{j}{2\sqrt 2}\right)^2+\left(\frac{j}{2\sqrt 2}\right)^2\\
(0.5d)^2 &= 2\left(\frac{j}{2\sqrt 2}\right)^2\\
0.5d &= \sqrt2\left(\frac{j}{2\sqrt 2}\right)\\
0.5d &= \frac{\sqrt2j}{2\sqrt 2}\\
0.5d &= \frac{j}{2}\\
d &= j
\end{align}

$$


or in the two-neighborhood around the adjacent cell.
If $p_{i+1}$ is in the two neighborhood, the adjacent cell can be non-empty, from which the adjacent cell still is used as the next hop in the graph and then the 

## Algorithm
Add all points to an uniform grid with cell size $\frac{j}{2\sqrt 2}$.
The cell containing $s$ becomes our start node. 
 - For each adjacent cell $a$, if $a$ contains nodes, add $a$ to the graph with distance $1$. 
 - otherwise check all nodes in cells in the two neigborhood around a but not in the one neigborhood around s if they are in jump distance to any node in a. if yes, add that cell to the graph with edge cost 1.
	- stop this search if any connection is found as multiple are not needed
	- [ ] find speedup approaches for this

this way we build the graph our dijkstra is running on while running the dijkstra.

![[Pasted image 20250430093455.png]]

#### Optimization
 - Given a search graph solution, the shortest path is searched by running dijkstra on only the nodes in the cells that are pare of the shortest path in the search graph.
 - If any adjecent cell is empty, the check if a node in another cell is in jump range can be done in $O((p+n) \log n )$ with $p$ being the number of points in the other cell and $n$ the number of points in the source cell, by creating a kd tree on the source cell and then executing a closest neighboor query on each point in the target cell until one is close enough. The kd tree can be reused for other target cells.
 - the check if a cell is empty can be done like this:
   store all x coords in a tree. each node contains a sorted list of all y coords if the points with this x,y, combination is closer than the max jump range of all ships to the x coord line. for performance calculation we assume max = $\infty$. 
   per binary search find the first x coord within the rectangle ($\log n$) then with binary search, search its y values. return TRUE if a point is found, go to the next higher x value if not. return FALSE if the rectangle is left. this should be possible in $O(\log^2 N)$

### Time Complexity


If we have a cell size $\frac{j}{2\sqrt 2}^2$ and $N$ points and all cells are edge alinged with a rectangle of size $R$ that contains all points and these points are uniformly distributed in the rectangle
The probability for any given cell with arriving at least 1 node:
$$
P(\exists p \in C) = 1-P(p\notin C) = 1 - \left(1- \frac{R^2}{\frac{j}{2\sqrt 2}^2}\right)^N
$$

Assuming each case shown in the proof above has the same probability, 

| case | subcase        | probability                                                                                                                                               | hops per optimal hop |
| ---- | -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------- |
| 1    | -              | $\frac14$                                                                                                                                                 | 1                    |
| 2    | -              | $\frac14$                                                                                                                                                 | 1                    |
| 3    | A empty        | $\frac14\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)$                                                                         | 1                    |
| 3    | A non-empty    | $\frac14\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N$                                                                                          | 2                    |
| 4    | A, B empty     | $\frac14\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)$ | 1                    |
| 4    | B empty        | $\frac14\cdot\left(\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)$     | 2                    |
| 4    | A, B non-empty | $\frac14\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N$                                   | 3                    |
|      |                |                                                                                                                                                           |                      |

Therefore the expected number of hops for a path is:
$$
\begin{align}
E(X) &= [0.25\cdot 1] + [0.25\cdot 1] + [\frac14\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right) \cdot 1] + [\frac14\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N \cdot 2] + [\frac14\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right) \cdot 1] + [\frac14\cdot\left(\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)\cdot\left(1 - \left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\right)\cdot 2] + [\frac14\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\cdot\left(1- \frac{\frac{j}{2\sqrt 2}^2}{R^2}\right)^N\cdot 3]
\\
&= 1 + \frac34 (1 - \frac{j^2}{8 R^2})^N

\end{align}
$$
As $j << R$, we now that $0 < \frac{j^2}{8 R^2} < 1 \implies 0 < 1 - \frac{j^2}{8 R^2}< 1 \implies 0 < (1 - \frac{j^2}{8 R^2})^N< 1$
Therefore we know
$$
1 < E(X) < 1.75
$$

For each point in the path. We check for each adjacent rectangle if it is empty which is $O(8\cdot \log^2 d)$ with $d$ being the average amount of points in a cell (density).
If a point is empty, which has the probability $P(cell empty) = e^{-d}$, up to 10 cells, are checked if they are empty. Each check takes $O(d \log d)$

This means each step in the exploration takes
$$
8\cdot (\log^2 d + e^{-d}\cdot (d \log d))
$$
In the dijkstra this is the amount of times needed to explore a node. therefore the total time for the dijkstra is:
$$
V(8\cdot (\log^2 d + e^{-d}\cdot (d \log d))) + (V+E)\log V
$$
But $V$ is not the number of points but the number of grid cells, from which we have $\left\lceil\frac {|points|}d\right\rceil$. Also from the cells we know that we can have at most 45 edges per node:
$$
\left\lceil\frac {|points|}d\right\rceil(8\cdot (\log^2 d + e^{-d}\cdot (d \log d))) + \left(\left\lceil\frac {|points|}d\right\rceil+45\left\lceil\frac {|points|}d\right\rceil\right)\log \left\lceil\frac {|points|}d\right\rceil
$$
$$
\boxed{\mathcal O\left(\left\lceil\frac {|points|}d\right\rceil(\log^2 d + e^{-d}\cdot (d \log d)) + \left\lceil\frac {|points|}d\right\rceil\log \left\lceil\frac {|points|}d\right\rceil\right)}
$$




```tikz 
\begin{document} 
\begin{tikzpicture}[domain=0:3, scale=1.8] 

    % Draw grid

    \draw[step=1cm,gray,very thin] (-5,-5) grid (5,5);

    % Draw the point at (0,0)

    \filldraw[red] (0,0) circle (2pt);
    \draw[red,thick] (0,0) circle (2.828);

	\filldraw[blue] (1,1) circle (2pt);
    \draw[blue,thick] (1,1) circle (2.828);

	%\filldraw[green] (0,0) circle (2pt);
    %\draw[green,thick] (0,0) circle (2.828);
\end{tikzpicture} 
\end{document} 
```

