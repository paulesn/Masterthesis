
> [!summary]
> This proof assumes the following visibility graph definition:
> given a set of points $V\subset R^3$ and a set of edges $P \subset E = \{(a,b)| a,b\in V\times V\}$ a visibility graph is a subgraph of $G = (V,E)$ in which all edges in $E$ exist except if they are crossing an edge in $P$.
> As this definition is stupid, I also add that each cyclic path in $P$ has an inside and an outside. The visibility graph has no edge that is at any point inside such a circle.

## n = 7

Given the following point set with the red outlined visibility blocking polygons:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{5}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]

\foreach \i in {1,...,\max}{
	\node[base, blue] (\i) at (\i*\a:2) {};
	\filldraw[gray] (\b + \i*\a:.5) -- (\b+\c+\i*\a:5) -- (\b-\c+\i*\a:5) --cycle;
}
\end{tikzpicture}
\end{document}
```


All red edges should be of the same length $d$. The green edge should be of length $\frac d2$:
```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\renewcommand{\max}{5}
\renewcommand{\a}{360/\max}
\renewcommand{\b}{\a/2}
\renewcommand{\c}{\b/2}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]

\foreach \i in {1,...,\max}{
	\node[base, blue] (\i) at (\i*\a:2) {};
	\filldraw[gray] (\b+\i*\a:1) -- (\b+\c+\i*\a:4) -- (\b-\c+\i*\a:4) --cycle;
}

\draw[green] (\b+\c+0*\a:4) -- (\b-\c+0*\a:4) -- (0*\a:2) -- (0*\a:2);
\draw[red]  (0*\a:2) --(0,0);
\end{tikzpicture}
\end{document}
```
Then we have four kinds of nodes:
- free nodes (blue)
- center circle nodes
- outside nodes 

lets define the distance of a red edge as $d$ 

## Center Node 
given a set that includes a center node, the smallest circle including two nodes would be one center node and one blue node. 
 - The 