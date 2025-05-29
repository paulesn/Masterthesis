When using the Well seperated pairs decomposiotion to build the subgraph, we have multiple methods of coosing the representant of the well seperated pairs. 
Here are 3 approaches and their proofs listed

## WSPD Construction Center Representant

> [!tip] Lemma
> Given a Graph $G = (V,E)$, a $t$-Spanner $S = (V,E_S)$ for $G$ constructed with the WSPD construction method, $s=4+t$ and a threshold  $c$. The subgraph $S' = (V, \{e \in E_S | d(e)\leq 2c\})$ contains a path from $a \in V$ to $b \in V$ iff the  Graph $G' = (V_S, \{e \in E | d(e)\leq c\})$ contains a path from $a$ to $b$
> 
> This path is at most twice as long as the same path in $G'$

^4d33e0
This proof assumes the following wspd-spanner construction:
> From top down buid a quad tree. Each cell with each cell is tested if it is well separated.
> If they are in the cell the center point of the enclosing circle is choosen as a representant. the representants are connected. for each cell all members are connected with their representant. 
> if the cell pair is not well seperated test it with each pairing of the cell and the children of the other cell as well as all children pairings recusrsivly.

$\exists (a \to b) \in G' \implies \exists (a \to b) \in S'$:
We have a path $\pi$ from $a$ to $b$ in $G'$. With $\pi = p_0, p_1, ..., p_n$ with $p_0 = a$ and $p_n = b$.
For any $0<i \leq n$ there is an edge $(p_{i-1},p_i)$ with length shorter than $c$. Otherwise the path would not be in $G'$.
By Definition, there exists a well separated pair of points $P_i$ and $P_{i-1}$ with $p_i$ and $p_{i-1}$ as respective members. Let $r_i$ and $r_{i-1}$ the center points of $P_i$ and $P_{i-1}$.
Then, at most $p_i$ and $p_{i-1}$ can be on a line between $r_i$ and $r_{i-1}$ equally far away from their represenant and both the furthest point away from their representant. Then, the circle enclosing  $r_i$ and $r_{i-1}$ has the radius $r= |p_ir_i|$ and $|r_ir_{i-1}| = 2r+|p_ip_{i-1}|$.

$$
\begin{align}
2r+|p_ip_{i-1}| &= s\cdot r &| -2r\\
|p_ip_{i-1}| &= s\cdot r -2r\\
|p_ip_{i-1}| &= (s-2)r &|s=4+t\\
|p_ip_{i-1}| &= (4+t-2)r \\
|p_ip_{i-1}| &= (2+t)r \\
\end{align}
$$
that means that the edge from $p_i$ to its representative is still in $S'$. And given an Edge $e$ in $G'$ its representative edge is still in $S'$.
Therefore, there is always a path $p_{i-1}, r_{i-1},r_i, p_i$ in $S'$, that connects any two nodes of a path still in $G'$. Per induction, this shows that for each edge in $\pi$ there exists a path in $S'$.

We know that the proxy route for a jump of length $d$ is at most $2r+d$. 
With $r= \frac{d}{2+t}$. therefore it is 
$$
\begin{align}
|p_ip_{i-1}|_{S'} &= \frac{2d}{2+t}+d \\
&= \frac{2d+d(2+t)}{2+t} \\
&= \frac{d(2(2+t))}{2+t} \\
&= \frac{d(4+2t)}{2+t}\\
&= d\frac{2(2+t)}{2+t}\\
&= 2d\\
\end{align}

$$

**q.e.d.**

```tikz
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (0) at (0,0) {};
\node[below of = 0] (0t){$r_{i-1}$};
\draw (0,0) circle (3);

\node[fill, circle] (1) at (3,0) {};
\node[below of = 1] (1t) {$p_{i-1}$};

\node[fill, circle] (2) at (6,0) {};
\node[below of = 2] (2t){$p_{i}$};

\node[fill, circle] (3) at (9,0) {};
\node[below of = 3] (3t){$r_{i}$};
\draw (9,0) circle (3);

\draw[-] (1) to node[below] {$|p_{i-1}p_{i}|$} (2);
\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```


---

## WSPD Construction by Book

> [!tip] Lemma
> Given a Graph $G = (V,E)$, a $t$-Spanner $S = (V,E_S)$ for $G$ constructed with the WSPD construction method, $s=4+t$ and a threshold  $c$. The subgraph $S' = (V, \{e \in E_S | d(e)\leq 1.5c\})$ contains a path from $a \in V$ to $b \in V$ if the  Graph $G' = (V_S, \{e \in E | d(e)\leq c\})$ contains a path from $a$ to $b$


this time we use the wspd spanner construction from the Spanner Book [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf]]:
> The WSPD is constructed top down. For each well separated pair, a random representant is choosen and an edge is placed between them.
> The distance between the two circles in the WSPD is defined as the distance between the circles and not the midpoints of the circles

```tikz
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (0) at (-0.2,0.3) {};
\node[below of = 0] (0t){$r_{i-1}$};
\draw (0,0) circle (3);

\node[] (1) at (3,0) {};
%\node[below of = 1] (1t) {$p_{i-1}$};

\node[] (2) at (6,0) {};
%\node[below of = 2] (2t){$p_{i}$};

\node[fill, circle] (3) at (9.5,-0.1) {};
\node[below of = 3] (3t){$r_{i}$};
\draw (9,0) circle (3);

\draw[<->] (1) to node[below] {$s\cdot r$} (2);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```
Proof: ^bd63d2

$\exists (a \to b) \in G' \implies \exists (a \to b) \in S'$:
for each edge $(p_1,p_2)$ in the path $(a\to b)$ in $G'$ There is a well seperated pair with (w.l.o.G.) $p_1 \in A$ and $p_2 \in B$.
This is the worst case:

```tikz
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (0) at (-3,0) {};
\node[below of = 0] (0t){$r_1$};
\draw (0,0) circle (3);
\node at (0,0) {A};

\node[fill, circle] (1) at (3,0) {};
\node[below of = 1] (1t) {$p_1$};

\node[fill, circle] (2) at (6,0) {};
\node[below of = 2] (2t){$p_2$};

\node[fill, circle] (3) at (12,0) {};
\node[below of = 3] (3t){$r_2$};
\draw (9,0) circle (3);
\node at (9,0) {B};

\draw[<->] (1) to node[below] {$c = s\cdot r$} (2);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```

The two nodes in our path $p_1$ and $p_2$ are exactly the jump range $j=s\cdot r$ apart from each other and they are both at the opposite ends of their respective circles to the representants of their set. Therefore they are both $2r$ from their representants apart. and $|r_1r_2| = (2+s)r$.

Because $s$ has to be at least $>4$ (lets say it is $4+o\in \mathbb{R_+}$):
$$
\begin{align}
|r_1,r_2| = (4+s)r &= (4+4+o)r\\
(4+s)r &= (8+o)r\\
\\
j = s\cdot r &= (4+o)r\\
\\
c'(4+o)r &> (6+o)r\\
1.5j &> |r_1,r_2|
\end{align}
$$

Now we have to know if there is an edge between $p_1$ and $r_1$ as well as  $p_2$ and $r_2$. 
Without loss of generality we only look at $p_1$ and $r_1$.
These points are at most $\frac j4$ apart from each other. Therefore there is an Edge  $(p_1, r_1) \in G'$.
Therefore this argumentation can be repeated recursively until $A$ and $B$ each contain only one node. Then $p_1 = r_1$.



---

## WSPD Construction Closest Pair Representant

> [!tip] Lemma
> Given a full Graph $G = (V,E)$, a $t$-Spanner $S = (V,E_S)$ for $G$ constructed with the WSPD construction method, $s>4$ and a threshold  $c$. The subgraph $S' = (V, \{e \in E_S | d(e)\leq c\})$ contains a path from $a \in V$ to $b \in V$ iff the  Graph $G' = (V_S, \{e \in E | d(e)\leq c\})$ contains a path from $a$ to $b$
> 
> $S'$ is an $2-$Spanner of $G'$


this time we use the wspd spanner construction from the Spanner Book [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf]]:
> The WSPD is constructed top down. For each well separated pair, a random representant is choosen and an edge is placed between them.
> The distance between the two circles in the WSPD is defined as the distance between the circles and not the midpoints of the circles
> **BUT** we modify it that the representant choosen for the edge is always the pair closest between $A$ and $B$.

Proof:

$\exists (a \to b) \in S' \implies \exists (a \to b) \in G'$:
Because $G$ is a full graph, every graph using only $V$ as nodes is a subgraph of $G$. Therefore $G'$ is a subgraph of $G$ that contains all edges with a length shorter than $c$ and $S'$ is a subgraph of $G$ that contains some edges with a length shorter than $c$ but none that are longer than $c$. Therefore $S'$ is also an subgraph of $G'$.
Because $S'$ is a subgraph of $G'$, if a path exists in $S'$ it also exists in $G'$.

$\exists (a \to b) \in G' \implies \exists (a \to b) \in S'$:

For each edge $(p,q)$ in the path there is an well seperated pair $A$ and $B$ with (w. l . o. G.) $p\in A$ and $p\in B$. 
Per Definition the edge $(r_A,r_B)$ in $S$ which is connecting $A$ and $B$ is at most as long as any other edge between two points $p' \in A$ and $q' \in B$. Therefore $A$ and $B$ are connected in $S'$. if there is any connection in $A$ and $B$  in $G'$.
By definition, $(r_A,r_B)$ is at least s times longer than any edge between two points in $A$ therefore, if $(r_A,r_B)$ is part of $S'$, all Edges $\{(v_1,v_2)| v_1,v_2 \in A\}$ are in $G'$.
Then the existence of a path $(p\to r_A)$ is argumented recursively until $A$ only contains a single element. At this point $p = r_A \implies \exists (p \to r_A)$. 
That means in if $p$ and $q$ are connected in $G'$ then there is a path $(p,r_A, r_B, q)$ in $S'$.

For each edge $(p,q)$ in the shortest path in $G'$, the spanner needs up to three components: $(p\to r_A)(r_A, r_B)(r_B\to q)$. Based on the spanner property, the edge $(p,r_A)$ and $(r_B, q)$ are at most $\frac{|(r_A, r_B)|}s$ long. That means in our spanner $S'$ they are at most $\frac{t'|(r_A, r_B)|}s$ The worst case for $|(r_A, r_B)|$ is $|(p,q)|$ with $(r_A, r_B) \neq (p, q)$. (Because if $(r_A, r_B) = (p, q)$ the other components would have a length of $0$)
Then the total length of the worst path that represents $(p,q)$ in $S'$:

$$
\begin{align}
d_{S'}(p,q) &\leq 2t'\frac{d_{G'}(p,q)}s + d_{G'}(p,q)\\
d_{S'}(p,q) &\leq 2t'\frac{d_{G'}(p,q)}s + \frac{s\cdot d_{G'}(p,q)}s\\
d_{S'}(p,q) &\leq \frac{2t'd_{G'}(p,q) + s\cdot d_{G'}(p,q)}s\\
d_{S'}(p,q) &\leq \frac{d_{G'}(p,q)(2t'+s)}s\\
t'd_{G'}(p,q) &\leq \left(\frac{2t'}s+1\right) d_{G'}(p,q)\\
t'&\leq \left(\frac{2t'}s+1\right)\\
t'&\leq \left(t'\frac{2}s+1\right) &s>4\implies \frac2s<\frac12\\
t'&\leq \left(t'\frac{2}s+1\right) < \left(t'\frac{1}2+1\right)\\
t'&\leq  \left(\frac{t'}2+1\right) &|\cdot 2\\
2t'&\leq  t'+2 &|-t'\\
t'&\leq  2 \\
\end{align}
$$

**q.e.d.**

- [ ] S' mit c+ epsilon filtern damit man bei der konstruktion keine vergleiche machen muss
	- [ ] visiblity graphen nachdenken
- [ ] 

---






- c' sollte von c abhängig sein -> Beweis nochmal ordentlich durchrechnen
- Kreis mit barrieren zwischen jedem knoten dass man nur gegenüber sieht 