
1. $L = \{ w \in \{a, b\}^* \mid w = w^R \}$
```tikz
\begin{document}
\begin{tikzpicture}[scale=3]

\node[draw, circle, minimum size=1cm] (0) at (0,0) {0};
\node[draw, circle, minimum size=1cm] (1) at (1,0) {1};
\node[draw, circle, minimum size=1cm] (2) at (.5,-.75) {2};
\node[draw, circle, minimum size=1.2cm] (3) at (.5,-.75) {};

\draw[->] (-.25,.25) to (0);
\draw[->, in=65, out=105] (0) to node[above,align=center]{$a,a,aa$\\$a,b,ba$\\$b,a,ab$\\$b,b,bb$} (0);
\draw[->, in=65, out=105] (1) to node[above,align=center]{$a,b,e$\\$ b,a,e$} (1);
\draw[->] (0) to node[above,align=center]{$a,b,e$\\$ b,a,e$} (1);
\draw[->] (1) to node[left,align=center]{$e,e,e$} (3);

\end{tikzpicture}
\end{document}
```

2. $L = \{ w \in \{a, b\}^* \mid w_a = w_b \}$
3. ends with ab
4. only contains a
5. contains exactly 4 a's
6. contains never aabbaa
7. just FIUS
8. all words with infinit length
9. epsilon
10. a joker