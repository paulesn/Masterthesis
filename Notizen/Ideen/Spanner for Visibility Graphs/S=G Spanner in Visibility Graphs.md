> [!tip] Lemma 
> Given a number $n \in \mathbb N$ one can place $n$ nodes together with visibility blockers in $\mathbb R^2$ so that the resulting visibility graph $G$ has only one connected component and the WSPD-spanner $S$ (using a separation factor $s>2$) based on $G$ are the same graph.
> For large $n: |E| \in \mathcal{O}(n^2)$

$G = S$ iff all well separated pairs are pairs with only one node in both sets. 

- $n=1$
  $G = S$ because neither has edges
- $n=2$
  $G = S$ because both only have one edge
- $n=3$
  place all points equidistant from each other (a triangle). Then each pair between 2 and 1 node is as far apart as the two nodes are in the set.
- $n=4$
  place the points in a square. any two points that can be in a set together are either adjacent or diagonally to each other
	- **adjacent**: for each other node in $G$ there is one node in the set that is equally far apart than the nodes in the set are far apart from each other
	- **diagonal**: the nodes in the set are further apart from each other than from any other node in the graph
	$\implies G=S$
	Alternatively, one could place the fourth point in the center of the triangle created by the other points. see $n=5$.
- $n = 5$
  Space four nodes in a square and one in the center of the square:
-
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{4}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	%\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}
\node[draw, circle] (0) at (0,0) {0};

\end{tikzpicture}
\end{document}
```
Now if a set of two nodes is used to create a well separated pair, the nodes can either be adjacent on the outside of the square $((1,4),(1,2),(2,3)(3,4))$, $0$ is part of the set or the nodes in the set are on opposite side of the square.
	- **adjacent**: for each other node in $G$ there is one node in the set that is equally far apart than the nodes in the set are far apart from each other or the last node is $0$, which is closer.
	- **0 is part**: all nodes are equidistant from $0$. therefore no node can be well separated from a set including $0$.
	- **opposite** the nodes in the set as far apart as any two nodes in the graph can be.
	$\implies G=S$

- $n=6$


```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.65) to (\b+\i*\a:1.9);
}

\draw[-] (1*\a:2) to node[right] {a} (3*\a:2);
\draw[-, red] (1*\a:2) to (2*\a:2);
\end{tikzpicture}
\end{document}
```
  
  Given the graph $G$ shown above. The lines represent visibility blocker on the line between each adjacent node. Therefore they can only block the visibility between adjacent nodes.
  For no well separated pair with more than one node in a set to form, we need the distance between two nodes that are not adjacent $a$ has to be chosen so that $s\cdot a$ is larger than the diameter $2r$ of the circle.
  Then, for any pair of nodes $a,b$ (here $(0),(1)$, because neighbors neighbor is the closest node) and for any other node $c$ ( here $(4)$ or $(5)$, $2r$ being the distance to the furthest node away)  $|ab| > min(s|ac|,s|bc|)$ and therefore no two nodes can be in a set that is well separated to any node. therefore only single node cluster are well separated.
  
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,2,3,4,5,6}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.65) to (\b+\i*\a:1.9);
}
\draw[-] (1*\a:2) to node[inner sep=0pt,text width=0pt] (a) {} (3*\a:2);
\node[fill, circle, inner sep=0pt,text width=3pt,] at (0,0) {};
\draw[-] (0,0) to node[right] {r} (1*\a:2);
\draw[-] (2*\a:2) to (a);
\draw[-] (1*\a:2) to node[left] {?} (a);


\draw[-, red] (1*\a:2) to (2*\a:2);
\end{tikzpicture}
\end{document}
```
  
  With $6$ nodes we calculate the length of $|(1)(3)|$ with:
  $$
  \begin{align}
	|(1)(3)| &= 2\left(\underbrace{2\sin\left(\frac{180}6\right)}_\text{length of a side of the polygon}\cdot \sin\left(\underbrace{\frac{(6-2)\cdot180}{6}}_\text{inside angle}\div 2\right)\right)\\
	 &= 2\left(2\sin(30)\cdot \sin\left(\frac{4\cdot180}{6}\div 2\right)\right)\\
	 &= 2\left(2\sin(30)\cdot \sin\left(\frac{720}{6}\div 2\right)\right)\\
	 &= 2\left(2\sin(30)\cdot \sin(120\div 2)\right)\\
	 &= 2(2\sin(30)\cdot \sin(60))\\
	 &= \sqrt{3}
\end{align}
  
  $$
  $s|(1)(3)|> 2r$ so that no well separated pair  with more than one node in a set can be created. With $r=1$ we get 
  $$
s>\frac{2}{|(1)(3)|}
   $$
## General case

The argument for $n=6$ can be generalized:
$$
  \begin{align}
	s>\frac{2}{|(a)(b)|} &= \frac{2}{2\left(\underbrace{2\sin\left(\frac{180}n\right)}_\text{length of a side of the polygon}\cdot \sin\left(\underbrace{\frac{(n-2)\cdot180}{n}}_\text{inside angle}\div 2\right)\right)}\\
\end{align}
  $$


With this we can create the following table:

|      | Length $d$ of a side of the polygon | Inside Angle $\alpha$ of the polygon | $2d\sin(\alpha) = \|(1)(3)\|$ | $2/\|(a)(b)\| = s$ |
| ---- | ----------------------------------- | ------------------------------------ | ----------------------------- | ------------------ |
| $6$  | $1$                                 | $120$                                | $\approx1.732$                | $\approx1.155$     |
| $7$  | $\approx 0.868$                     | $\approx128.571$                     | $\approx1.564$                | $\approx1.279$     |
| $8$  | $\approx 0.765$                     | $135$                                | $\approx1.414$                | $\approx1.414$     |
| $9$  | $\approx0.684$                      | $140$                                | $\approx1.286$                | $\approx1.556$     |
| $10$ | $\approx0.618$                      | $144$                                | $\approx1.176$                | $\approx1.701$     |
| $11$ | $\approx0.563$                      | $\approx147.273$                     | $\approx1.081$                | $\approx1.850$     |
| $12$ | $\approx0.518$                      | $150$                                | $\approx 1$                   | $\approx2$         |
| $13$ | $\approx0.479$                      | $\approx152.308$                     | $\approx0.929$                | $\approx2.152$     |
| $14$ | $\approx0.445$                      | $\approx154.286$                     | $\approx0.8687$               | $\approx2.305$     |
| $15$ | $\approx0.416$                      | $156$                                | $\approx0.813$                | $\approx2.459$     |
| $16$ | $\approx0.390$                      | $157.5$                              | $\approx0.765$                | $\approx2.613$     |
| $17$ | $\approx0.368$                      | $\approx158.824$                     | $\approx0.722$                | $\approx2.768$     |
| $18$ | $\approx0.347$                      | $160$                                | $\approx0.684$                | $\approx2.924$     |
| $19$ | $\approx0.329$                      | $\approx161.053$                     | $\approx0.649$                | $\approx3.080$     |
> [!attention] 
> The value for $s$ in the table is the maximum $s$ with which the WSPD still contains sets with more than one node.
## Addition of valid n's

Given $n_1$ and $n_2$ with size $6$ or larger, we can combine them by placing them on top of each other with shifted degree. Additionally we place visibility blocker between all nodes that are adjacent and from different $n$'s:


```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[fill, circle, inner sep=0pt,text width=3pt,] (\i) at (\i*\a:2) {};
	\node[blue, fill, circle, inner sep=0pt,text width=3pt,] (b\i) at (\b+\i*\a:2) {};
	\draw[-] (\b+\i*\a:1.7) to (\b+\i*\a:1.8);
	\draw[-] (\i*\a:1.7) to (+\i*\a:1.8);
	\draw[-] (\c+\i*\a:1.95) to (\c+\i*\a:2.05);
	\draw[-] (\b+\c+\i*\a:1.95) to (\b+\c+\i*\a:2.05);
}
\end{tikzpicture}
\end{document}
```

As the visibility blocker can be placed freely on the line that connects the two nodes they are suppost to block, they can be moved so that they are all in line with another visibility blocker and the center of the circle:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[fill, circle, inner sep=0pt,text width=3pt,] (\i) at (\i*\a:2) {};
	\node[blue, fill, circle, inner sep=0pt,text width=3pt,] (b\i) at (\b+\i*\a:2) {};
	\draw[-] (\c+\i*\a:1.7) to (\c+\i*\a:2.05);
	\draw[-] (\b+\c+\i*\a:1.7) to (\b+\c+\i*\a:2.05);
}
\end{tikzpicture}
\end{document}
```

Now, the furthest point from each point is the same as the one before the two $n$'s were combined. and the same is true for the closest point for each point. Therefore no multi node well separated pair can be found.

This neither needs a equidistant shift:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{20}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[fill, circle, inner sep=0pt,text width=3pt,] (\i) at (\i*\a:2) {};
	\node[blue, fill, circle, inner sep=0pt,text width=3pt,] (b\i) at (\b+\i*\a:2) {};
	\draw[-] (\c+\i*\a:1.7) to (\c+\i*\a:2.05);
	\draw[-] (\b+\c+\i*\a:1.7) to (\b+\c+\i*\a:2.05);
}
\end{tikzpicture}
\end{document}
```

Nor do $n_1$ and $n_2$ need to be of the same size:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[fill, circle, inner sep=0pt,text width=3pt,] (\i) at (\i*\a:2) {};
	%\draw[-] (\i*\a:1.7) to (\c+\i*\a:2.05);
	\draw[-] (\b+\i*\a:1.7) to (\b+\i*\a:1.81);
}


\renewcommand{\max}{8}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{40}

\foreach \i in {1,...,\max}{
	\node[blue, fill, circle, inner sep=0pt,text width=3pt,] (b\i) at (\c+\i*\a:2) {};
	%\draw[-] (\c+\i*\a:1.7) to (\c+\i*\a:2.05);
	\draw[-, blue] (\c+\b+\i*\a:1.8) to (\c+\b+\i*\a:1.9);
}

\foreach \i in {-2.5,30, 50, 75, 90, 115,125, 145,177.5, 210, 230, 255, 270, 295,305, 325}{
	\draw[-, red] (\i:1.95) to (\i:2.05);
}

%\draw[-, red] (\c+0*\a:2) to (\c+1*\a:2);
\end{tikzpicture}
\end{document}
```
Because there are infinitly many points on a line between to points but only finite many points where two lines from a finite amount of points can cross, there can be always a placement for a visibility blocker, that is only blocking the intended point. 
The radius of the circle can be increased so that the visibility blocker can have an area without unintended blocking.

## Generalization to any n

Because we have all $n<12$, we can build the following table, which shows that for any $n$ we can construct a graph using addition and the graphs we already have: 

|       | $1$          | $2$          | $3$          | $4$          | $5$          | $6$          | $7$          | $8$          | $9$          |
| ----- | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |
| $10$  | $11$         | $6+6$        | $6+7$        | $7+7$        | $9+6$        | $10+6$       | $10+7$       | $10+8$       | $10+9$       |
| $10k$ | $10(k-1)+11$ | $10(k-1)+12$ | $10(k-1)+13$ | $10(k-1)+14$ | $10(k-1)+15$ | $10(k-1)+16$ | $10(k-1)+17$ | $10(k-1)+18$ | $10(k-1)+19$ |
A graph $G$ constructed this way will have a spanner $S=G$ as long as the separation factor $s>2$.


## Number of Edges 

If we assume the nodes are evenly spaced on the circle , at most $2k+1$ nodes are skipped by the shortest edge, if $k$ graphs are combined according to the table above.
For a graph with $n = 10k+a$ nodes each node still has $10k+a-2k+1 = 8k+a+1$ nodes. 
In total that are 
$$
\begin{align}
|E| &= n(8k+a+1)\\
&= (10k+a)(8k+a+1)\\
&= 10k(8k+a+1)+a(8k+a+1)\\
&= 80k^2+10ka+10k+8ka+a^2+a)\\
&= 80k^2+18ka+10k+a^2+a)\\
\end{align}
$$
Because $k\approx \frac{n}{10}$ for large $k$: 
$$
\begin{align}
&= 80\left(\frac{n}{10}\right)^2+18\left(\frac{n}{10}\right)a+10\left(\frac{n}{10}\right)+a^2+a)\\
\end{align}
$$
Because $a$ is a small constant, $\left(\frac{n}{10}\right)^2$ is dominating the term:
$$
80\left(\frac{n}{10}\right)^2 = 0.8n^2 \implies \mathcal O(n^2)
$$

nochmal zeigen: für ein fixes $s$ (z.B. 4) kann man einen Graphen bauen mit quadratischen spanner edges. (zwei linien gegenüber)
- erst wenn blocker nicht teil des graphen sind
- Dann wenn dei blocker teil des graphen sind