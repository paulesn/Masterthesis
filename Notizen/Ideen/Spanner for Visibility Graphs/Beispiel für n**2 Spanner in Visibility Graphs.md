## n = 1-5
- 1 is trivial
- 2 is trivial
- 3, if equally spaced, is trivial
- 4, if placed in a square, any two points that are in a set together are as close to each other as any potential node for the other side of the well separated pair
- 5
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{5}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}
\foreach \i in {1,...,\max}{
	
}

\draw[-] (1*\a:2) to (3*\a:2);
\end{tikzpicture}
\end{document}
```

all nodes a node can see are equally spaced

### n =6-12
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\a}{60}
\renewcommand{\b}{30}

\begin{document}
\begin{tikzpicture}

\node[draw, circle] (0) at (0*\a:1) {0};
\node[draw, circle] (1) at (1*\a:1) {1};
\node[draw, circle] (2) at (2*\a:1) {2};
\node[draw, circle] (3) at (3*\a:1) {3};
\node[draw, circle] (4) at (4*\a:1) {4};
\node[draw, circle] (5) at (5*\a:1) {5};

\draw[-] (\b+0*\a:.75) to (\b+0*\a:1.2);
\draw[-] (\b+1*\a:.75) to (\b+1*\a:1.2);
\draw[-] (\b+2*\a:.75) to (\b+2*\a:1.2);
\draw[-] (\b+3*\a:.75) to (\b+3*\a:1.2);
\draw[-] (\b+4*\a:.75) to (\b+4*\a:1.2);
\draw[-] (\b+5*\a:.75) to (\b+5*\a:1.2);

\draw[-] (0) to node [below] {d} (2);

\end{tikzpicture}
\end{document}
```

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{10}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}
\foreach \i in {1,...,\max}{
	
}

\draw[-] (1*\a:2) to (3*\a:2);
\end{tikzpicture}
\end{document}
```

Given the graph $G$ shown above. The lines represent visibility blocker exactly at the midpoint between each adjacent node. Therefore they can only block the visibility between adjacent nodes
The distance between two nodes that are not adjacent $a$ has to be chosen so that $s\cdot a$ is larger than the radius $r$ of the circle.
Then, for any pair of nodes $a,b$ (here $0,1$, because neighbors neighbor is the closest node) and for any other node $c$ ( here $4$, $2r$ being the distance to the furthest node away)  $|ab| > min(s|ac|,s|bc|)$ and therefore no two nodes can be in a set that is well separated to any node. therefore only single node cluster are well separated.

With $n$ nodes the distance between 2 adjacent nodes is $2r\sin(\frac{\pi}n)$. And as the relative distances in a circle are always (relative to each other) the same we can choose $r=1$

With $n$ nodes, the angle $\alpha$ between the edge $(0,1)$ and $(1,2)$ is $\frac{360}{n}$.
By the law of cosines: $a^2 = b^2 + c^2 - 2bc \cdot \cos(\alpha)$

With:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\a}{60}
\renewcommand{\b}{30}

\begin{document}
\begin{tikzpicture}

\node[draw, circle] (0) at (0*\a:1) {0};
\node[draw, circle] (1) at (1*\a:1) {1};
\node[draw, circle] (2) at (2*\a:1) {2};

\draw[-] (0) to node [below] {a} (2);
\draw[-] (0) to node [right] {b} (1);
\draw[-] (1) to node [above] {c} (2);

\end{tikzpicture}
\end{document}
```
$$
\begin{align}
	a^2 &= b^2 + c^2 - 2bc \cdot \cos(\alpha)\\
	a &= \sqrt{b^2 + c^2 - 2bc \cdot \cos(\alpha)}\\
	\text{with } b = c \\
	a &= \sqrt{b^2 + b^2 - 2b^2 \cdot \cos(\alpha)}\\
	a &= \sqrt{b^2 (1 + 1 - 2 \cdot \cos(\alpha))}\\
	a &= \sqrt{b^2 (2 - 2 \cdot \cos(\alpha))}\\
	a &= \sqrt{2b^2(1-\cos(\alpha))}\\
	a &= b\sqrt{2(1-\cos(\alpha))}\\
	\text{with }  b = 2\sin(\frac{\pi}n)\\
	\text{with }  \alpha = \frac{360}{n}\\
	\text{with }  d = a\\
	d &= 2\sin(\frac{\pi}n)\sqrt{2(1-\cos(\frac{360}{n}))}\\
\end{align}
$$

Now we want 
$$
\begin{align}
s*d &> 2r\\
s2\sin(\frac{\pi}n)\sqrt{2(1-\cos(\frac{360}{n}))} &>2r
\end{align}
$$
Which is true for:$n = \{6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 21, 22\}$


## n >12

Given a graph with $n>12$ we can still construct a circle in which $G=S$.
We construct a circle with $12 <n<24$ and $n$ divisable by 3:

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{18}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}

\draw[-] (1*\a:2) to (4*\a:2);
\end{tikzpicture}
\end{document}
```


and place a visual blocker line from $m$ the midpoint between adjacent nodes $a$ and $b$ to the closest point on the line from $a$ to the other adjacent node of $b$. This way all edges which are not skipping two neighboors are visually blocked.

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{18}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}

\draw[-, red] (0*\a:2) to (3*\a:2) to (6*\a:2)to (9*\a:2)to (12*\a:2) to (15*\a:2) to (18*\a:2);
\draw[-, blue] (1*\a:2) to (4*\a:2) to (7*\a:2)to (10*\a:2)to (13*\a:2) to (16*\a:2) to (1*\a:2);
\draw[-, green] (2*\a:2) to (5*\a:2) to (8*\a:2)to (11*\a:2)to (14*\a:2) to (17*\a:2) to (2*\a:2);

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at ($(5,0)+(\i*\a:2)$) {\i};
	%\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}
\draw[-] ($(5,0)+(0*\a:2)$) to ($(5,0)+(1*\a:2)$) to($(5,0)+(2*\a:2)$) to($(5,0)+(3*\a:2)$) to($(5,0)+(4*\a:2)$) to ($(5,0)+(5*\a:2)$) to ($(5,0)+(6*\a:2)$);

\end{tikzpicture}
\end{document}
```

The visusal outline for each node is similar to the one of the graph with a third of the nodes. Therefore the closest node is at the same distance as in the graph with $n/3$ nodes and the furthest is still $2r$ away. Therefore still no well separated pair with more than one node can be constructed.

this argument can  be repeated recursivly so that we can construct a graph for any $n$ in 
$$
\{k3n | n\in\{1,2,3,4,5,6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 21, 22\}, k\in\mathbb N]\}
$$

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{15}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}
\begin{tikzpicture}[]

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at (\i*\a:2) {\i};
	\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}

\draw[-, red] (0*\a:2) to (3*\a:2) to (6*\a:2)to (9*\a:2)to (12*\a:2) to (15*\a:2) to (18*\a:2);
\draw[-, blue] (1*\a:2) to (4*\a:2) to (7*\a:2)to (10*\a:2)to (13*\a:2) to (16*\a:2) to (1*\a:2);
\draw[-, green] (2*\a:2) to (5*\a:2) to (8*\a:2)to (11*\a:2)to (14*\a:2) to (17*\a:2) to (2*\a:2);

\renewcommand{\max}{5}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\foreach \i in {1,...,\max}{
	\node[draw, circle] (\i) at ($(5,0)+(\i*\a:2)$) {\i};
	%\draw[-] (\b+\i*\a:1.75) to (\b+\i*\a:2.2);
}
\draw[-] ($(5,0)+(0*\a:2)$) to ($(5,0)+(1*\a:2)$) to($(5,0)+(2*\a:2)$) to($(5,0)+(3*\a:2)$) to($(5,0)+(4*\a:2)$) to ($(5,0)+(5*\a:2)$) to ($(5,0)+(6*\a:2)$);

\end{tikzpicture}
\end{document}
```


## Everything 






```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}
\usepackage{pgfplots}
\pgfplotsset{width=10cm,compat=1.9}

\renewcommand{\max}{10}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}

\begin{document}

\def\s{4}    % example value for s
\def\r{10}    % example value for r

\begin{tikzpicture}[]
\begin{axis}[
    domain=1:10, % reasonable n values, start from 3 (triangle)
    samples=400,
    xlabel={$n$},
    ylabel={$f(n)$},
    grid=major,
    thick
]

\addplot[
    blue
]

{ \s * (2*pi*\r / x) * sqrt(2*(1 - cos(2*pi/x))) };
\end{axis}
\end{tikzpicture}
\end{document}



```
