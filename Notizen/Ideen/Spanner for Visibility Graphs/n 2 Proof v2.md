## Visibility Blocker are not part of of the Graph 

> [!tip] 
>  Given a $n \in \mathbb N$, we can place $n$ Points and Visibilityblocker in $\mathbb R^2$ so that the spanner $S$ constructed using a WSPD with a shortest path metric on the visibility graph on the Points and Visiblityblocker has $n^2$ edges.


```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]

\foreach \i in {1,...,\max}{
	\node[base] (a) at (\i,0) {};
	\node[base] (1) at (\i,5) {};
	\draw[red] (\i+0.5, 0.15) -- (\i+0.5,-0.15);
	\draw[red] (\i+0.5, 5.15) -- (\i+0.5,4.75);
}

\draw[blue] (1,0) -- (\max, 5);
\draw[green] (1,0) -- (1, 5);
	 
    
  \end{tikzpicture}
\end{document}
```

We construct a graph as seen in the figure above. If $n$ is odd, it is added to the upper line of nodes.
Then the green line has the length $d$ and the nodes adjecent to each other have a distance $a$. The red lines are visibility blockers.

The green line represents the shortest possible edge in the graph and the blue line represents the longest possible edge in the graph.
The ends of the edges, that are not in the same node, are 
$$
(\left\lfloor \frac n2 \right\rfloor-1) a
$$
apart from each other.
By pytagoras the blue line has a length of:
$$
\sqrt{((\left\lfloor \frac n2 \right\rfloor-1) a)^2+d^2}
$$

Now we can choose $s$ such that 
$$
ds > \sqrt{((\left\lfloor \frac n2 \right\rfloor-1) a)^2+d^2}
$$

We can place $d=1$ and then for a given $n$ and $s$ we receive:
$$
\begin{align}
s &> \sqrt{((\left\lfloor \frac n2 \right\rfloor-1) a)^2+1}\\
s^2 &> ((\left\lfloor \frac n2 \right\rfloor-1) a)^2+1\\
s^2-1 &> ((\left\lfloor \frac n2 \right\rfloor-1) a)^2\\
\sqrt{s^2-1} &> (\left\lfloor \frac n2 \right\rfloor-1) a\\
\frac{\sqrt{s^2-1}}{\left\lfloor \frac n2 \right\rfloor-1} &> a\\
\end{align}
$$

Therefore, the wspd can not construct any well separated pair using two or more nodes in one of the sets and the resulting spanner has $n^2$ edges.


## Visibility Blocker are part of the Graph 

in this proof, the corners of the visibility blocking polygons are part of the graph:

### V1 like above

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]
\draw[red] (0.5,1)
\foreach \i in {1,...,\max}{
	-- (\i,0)
	-- (\i.5,1)
};

\draw[red] (0.5,4)
\foreach \i in {1,...,\max}{
	-- (\i,5)
	-- (\i.5,4)
};
	 
    
  \end{tikzpicture}
\end{document}
```

This Graph has 4 relevant distances:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{6}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]
\draw[red] (0.5,1)
\foreach \i in {1,...,\max}{
	-- (\i,0)
	-- (\i.5,1)
};

\draw[red] (0.5,4)
\foreach \i in {1,...,\max}{
	-- (\i,5)
	-- (\i.5,4)
};

\draw[green] (0.5,1) -- (1.5,1);
\draw[blue] (0.5,1) -- (0.5,4);
\draw[black] (0.5,1) -- (\max.5,4);
	 
    
  \end{tikzpicture}
\end{document}
```
- the redline $r$
- the green line $g$, which represents the spacing of two points in the polygon
- the blue line $b$, which represents the distance between the two polygons
- the white line $w$, which is the longest visibility line in the graph.

We want to choose $b, g$ and $r$ in so a way, that $rs > w$. (and $bs > w$ and $gs>w$)
The number of nodes is $n$. so 
$$
w = \sqrt{b^2+\frac{ng}{2}^2}
$$

### V2

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{2}
\renewcommand{\a}{4}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]

\renewcommand{\a}{1}
\foreach \i in {1,...,8}{
	\foreach \o in {1,...,8}{
		\draw[red] (\i*\a,\o*\a) -- ++(-1*\a,0) -- ++(0,-1*\a) -- ++(1*\a,0) -- ++(0,1*\a);
	};	
};

\renewcommand{\a}{2}
\foreach \i in {1,...,4}{
	\foreach \o in {1,...,4}{
		\draw[green] (\i*\a,\o*\a) -- ++(-1*\a,0) -- ++(0,-1*\a) -- ++(1*\a,0) -- ++(0,1*\a);
		\draw[fill=red] (\i*\a,\o*\a) -- ++(-1*\b,0) -- ++(0,-1*\b) -- ++(1*\b,0) -- ++(0,1*\b);
	};	
};

\renewcommand{\a}{4}
\foreach \i in {1,...,\max}{
	\foreach \o in {1,...,\max}{
		\draw[blue] (\i*\a,\o*\a) -- ++(-1*\a,0) -- ++(0,-1*\a) -- ++(1*\a,0) -- ++(0,1*\a);
		\draw[fill=green] (\i*\a,\o*\a) -- ++(-1*\b,0) -- ++(0,-1*\b) -- ++(1*\b,0) -- ++(0,1*\b);
	};	
};



    
  \end{tikzpicture}
\end{document}
```


given $n= x^3$ points, place them in the raster above (recursivly continued) so that each unfilled cell contains exactly one node. 
Then we have a cell size $c$ at the lowest level. the furthest distance a wsp can form containing a group of nodes is the 2nd lowest level of a corner cell to a lowest level corner cell. Given a depth $d$ we get a distance $\Delta$ 
$$
\Delta = \sqrt{2((\sqrt{d}-2)c)^2}
$$
in each filled cell, of side length $x$ is the following visibility blocking polygon:

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{2}
\renewcommand{\a}{4}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]

\draw[white] (-1,-1) -- ++(-12,0) -- ++(0,-12) -- ++(12,0) -- ++(0,12);
\draw[red] (0,0) -- ++(-10,0) -- ++(0,-10) -- ++(10,0) -- ++(0,10);

\draw[blue] (-1,-1) -- ++(-8,0) -- ++(0,-8) -- ++(8,0) -- ++(0,7) -- ++(-7,0) -- ++(0,-6) -- ++(6,0) -- ++(-6.1,-0.1)-- ++(0,6.2) -- ++(7.2,0) -- ++(0,-7.2) -- ++(-8.2,0) -- ++(0,8.2) -- ++(8.1,-0.1);



    
  \end{tikzpicture}
\end{document}
```

The longest shortest path in this polygon is $0.9c$ three times and then $2(0.9^dc)$ for each depth level, as the number of "circles" is increased with the depth of the quad tree, so that the distance inside the spiral is longer than $2d/s$ on the same level. therefore in the same level, the quad tree - wspd alorithm cannot find any well separated pair as each cell has a to high radius. 

on the lowest level than we connect each cell, assuming as lower bound the filled cells are well separated with each other cell.

$$
s\Delta = s\sqrt{2((\sqrt{d}-2)c)^2} > 3\cdot0.9c + \sum_{l=1}^d2(0.9^l)c
$$

Now we can calculate the largest $s$ this works for. We can scale down $s$ by adding an $x$:
$$
s\Delta = s\sqrt{2((\sqrt{d}-2)c)^2} > 3\cdot0.9c + \sum_{l=1}^{d+x}2(0.9^l)c
$$
The colored cells then have respective for their level $6+(4+2x)l$ nodes.
