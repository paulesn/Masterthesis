
## Visibility Graphs

For the context of visiblity graphs, this could be interesting:

 > [!todo] Lemma
> Given a Graph $G = (V,E)$, a $t$-Spanner $S = (V,E_S)$ for $G$ constructed with the WSPD construction method, $s=4+t$ and a threshold  $f: x\to \{0,1\}$. If the subgraph $S' = (V, \{e \in E_S | f(e)= 1\})$ contains a path from $a \in V$ to $b \in V$  it is at most X times longer than the shortest path from $a$ to $b$ in the Graph $G' = (V_S, \{e \in E | f(e)=1\})$

With this example graph:
 ```tikz 
\begin{document} 
\begin{tikzpicture}

\node[draw,circle] (p0) at (0,0){};
\node[draw,circle] (p1) at (1,.5){};
\node[draw,circle] (p2) at (2,.5){};
\node[draw,circle] (p3) at (1.5,-.5){};
\node[draw,circle] (p4) at (3,0){};
\node[draw,circle] (p6) at (0,-1){};
\node[draw,circle] (p7) at (1,-1.5){};
\node[draw,circle] (p8) at (2,-1.5){};
\node[draw,circle] (p9) at (3,-1){};

\draw[-] (p0) to node[above] {2} (p1)
              to node[above] {2} (p2)
              to node[above] {2} (p4)
              to node[right] {2} (p9)
              to node[below] {2} (p8)
              to node[below] {2} (p7)
              to node[below] {2} (p6)
              to node[left] {2} (p0);
              
\draw[-] (p0) to node[above] {1} (p3);
\draw[-] (p1) to node[right] {1} (p3);
\draw[-] (p2) to node[right] {1} (p3);
\draw[-] (p4) to node[below] {1} (p3);
\draw[-] (p9) to node[below] {1} (p3);
\draw[-] (p8) to node[left] {1} (p3);
\draw[-] (p7) to node[left] {1} (p3);
\draw[-] (p6) to node[above] {1} (p3);

\end{tikzpicture} 
\end{document} 
```
we get the following greedy spanner:
 ```tikz 
\begin{document} 
\begin{tikzpicture}

\node[draw,circle] (p0) at (0,0){};
\node[draw,circle] (p1) at (1,.5){};
\node[draw,circle] (p2) at (2,.5){};
\node[draw,circle] (p3) at (1.5,-.5){};
\node[draw,circle] (p4) at (3,0){};
\node[draw,circle] (p6) at (0,-1){};
\node[draw,circle] (p7) at (1,-1.5){};
\node[draw,circle] (p8) at (2,-1.5){};
\node[draw,circle] (p9) at (3,-1){};
              
\draw[-] (p0) to node[above] {1} (p3);
\draw[-] (p1) to node[right] {1} (p3);
\draw[-] (p2) to node[right] {1} (p3);
\draw[-] (p4) to node[below] {1} (p3);
\draw[-] (p9) to node[below] {1} (p3);
\draw[-] (p8) to node[left] {1} (p3);
\draw[-] (p7) to node[left] {1} (p3);
\draw[-] (p6) to node[above] {1} (p3);

\end{tikzpicture} 
\end{document} 
```

But if the filter function is chosen in a way so that it removes all Edges with $w=1$ in the original graph:
```tikz 
\begin{document} 
\begin{tikzpicture}

\node[draw,circle] (p0) at (0,0){};
\node[draw,circle] (p1) at (1,.5){};
\node[draw,circle] (p2) at (2,.5){};
\node[draw,circle] (p3) at (1.5,-.5){};
\node[draw,circle] (p4) at (3,0){};
\node[draw,circle] (p6) at (0,-1){};
\node[draw,circle] (p7) at (1,-1.5){};
\node[draw,circle] (p8) at (2,-1.5){};
\node[draw,circle] (p9) at (3,-1){};

\draw[-] (p0) to node[above] {2} (p1)
              to node[above] {2} (p2)
              to node[above] {2} (p4)
              to node[right] {2} (p9)
              to node[below] {2} (p8)
              to node[below] {2} (p7)
              to node[below] {2} (p6)
              to node[left] {2} (p0);
              

\node[draw,circle] (p10) at (6,0){};
\node[draw,circle] (p11) at (7,.5){};
\node[draw,circle] (p12) at (8,.5){};
\node[draw,circle] (p13) at (7.5,-.5){};
\node[draw,circle] (p14) at (9,0){};
\node[draw,circle] (p16) at (6,-1){};
\node[draw,circle] (p17) at (7,-1.5){};
\node[draw,circle] (p18) at (8,-1.5){};
\node[draw,circle] (p19) at (9,-1){};

\node[] at (1.5,1.5){filtered Graph:};
\node[] at (7.5,1.5){filtered Spanner:};
\end{tikzpicture} 
\end{document} 
```