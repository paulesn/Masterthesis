

The Proof for the upper bound of an theta spanner in [Spanner Buch] shows that an theta spanner is a $t$-Spanner with
$$
t = \frac{1}{\cos\theta-\sin\theta}
$$
With $\theta$ being the angle of the largest cone or $\frac{2\pi}{k}$. ($k$ being the number of equal sized cones).

This Worst-Case analysis assumes that the points $s,p,t$ relevant to the constructions are at extrem points in the cone, creating an angle between them equal to $\theta$:

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}

\begin{document}
\begin{tikzpicture}[
base/.style={fill, circle, inner sep=0pt,text width=3pt}
]
\node[base] at (0,0) (s) {};
\node[below of = s, node distance = .5cm] {s};
\node[above right of = s, node distance = 1cm, yshift = -.4cm] {$\theta$};
\node[base] at (0:5cm) (p) {};
\node[below of = p, node distance = .5cm] {p};
\node[base] at (40:5cm) (t) {};
\node[above of = t, node distance = .5cm] {t};

\draw[-] (s) to +(0:9cm);
\draw[-] (s) to +(40:9cm);
	 
    
  \end{tikzpicture}
\end{document}
```

But the angle is not always this extrem. The angle can be expected to follow a continuous normal distribution over the intervall $]0,\theta[$ over the whole graph. 
Then the expected angle $\theta'$ will be $\frac{\theta}{2}$.

> [!attention] 
> This proof assumes there are no visibility Blocker in the graph

