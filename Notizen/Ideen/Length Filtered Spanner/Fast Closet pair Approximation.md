DOESNT WORK
Given two circles set of points $A$ and $B$ and a radius $r$, for each set $S \in \{A,B\}$ we have a center point $c_S$ so that : $\forall s \in S: |c_ss|\leq r$.
We want to quickly find a pair $a,b$ so that $a \in A$ and $b \in B$ for that holds:
$$
|ab|\leq t \min(\{|a'b'| \mid a' \in A, b'\in B\})
$$


We have to sets $A$ and $B$ with center points.

```tikz
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (0) at (0,0) {};
\node[below of = 0] (0t){$c_A$};
\draw (0) circle (3);
\node at (0) {A};

\node[fill, circle] (3) at (10,0) {};
\node[below of = 3] (3t){$c_B$};
\draw (3) circle (3);
\node at (3) {B};

%\draw[<->] (3,0) to node[below] {$c = s\cdot r$} (6,0);
%\draw[<->] (0) to node[below] {$t'c = r+(s\cdot r)$} (3);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```

Now we select the middlepoint $m$ between $c_A$ and $c_B$. The closest point in $A$ and $B$ are on the blue and green line respectivly.:

```tikz
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (0) at (0,0) {};
\node[below of = 0] (0t){$c_A$};
\draw (0) circle (3);
\node at (0) {A};

\node[fill, circle] (3) at (10,0) {};
\node[below of = 3] (3t){$c_B$};
\draw (3) circle (3);
\node at (3) {B};

\node[fill, circle] (m) at (5,0) {};
\node[below of = m] (3t){$m$};
%\draw[red] (m)+(0,-4) rectangle ++(-4,4);
\begin{scope}
	\clip (m)+(0,-3) rectangle ++(3,3);
	\draw[green] (m) circle (3);
\end{scope}
\begin{scope}
	\clip (m)+(0,-4) rectangle ++(-4,4);
	\draw[blue] (m) circle (4);
\end{scope}
%\draw[<->] (3,0) to node[below] {$c = s\cdot r$} (6,0);
%\draw[<->] (0) to node[below] {$t'c = r+(s\cdot r)$} (3);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```

The worst case would be the two points on opposite intersection between the distance circle and the radius circle of the sets:
```tikz
\usetikzlibrary{calc}
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (a) at (0,0) {};
\node[below of = a] (at){$c_A$};
\draw (a) circle (3);
\node at (a) {A};

\node[fill, circle] (b) at (10,0) {};
\node[below of = b] (bt){$c_B$};
\draw (b) circle (3);
\node at (b) {B};

\node[fill, circle] (m) at (5,0) {};
\node[below of = m] (mt){$m$};
%\draw[red] (m)+(0,-4) rectangle ++(-4,4);
\begin{scope}
	\clip (m)+(0,-3) rectangle ++(3,3);
	\draw[green] (m) circle (3);
\end{scope}
\begin{scope}
	\clip (m)+(0,-4) rectangle ++(-4,4);
	\draw[blue] (m) circle (4);
\end{scope}
\node[fill=blue, circle] at ($(a)+(-53:3)$) {};
\node[fill=green, circle] at ($(b)+(147:3)$) {};
%\draw[<->] (3,0) to node[below] {$c = s\cdot r$} (6,0);
%\draw[<->] (0) to node[below] {$t'c = r+(s\cdot r)$} (3);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```

Then this creates this areal in which nodes could be that are closer pairs for either:

```tikz
\usetikzlibrary{calc}
\usetikzlibrary{through}
\begin{document}
\begin{tikzpicture}

\node[fill, circle] (a) at (0,0) {};
\node[fill, circle] (b) at (10,0) {};
\node[fill, circle] (m) at (5,0) {};

\node[fill=blue, circle] (iA) at ($(a)+(-55:3)$) {};
\node[fill=green, circle] (iB) at ($(b)+(165:3)$) {};

\begin{scope}
	\clip (m)+(0,-8) rectangle ++(7,8);
	\node [draw, green] at (m) [circle through={(iB)}] {};
\end{scope}
\begin{scope}
	\clip (m)+(0,-8) rectangle ++(-7,8);
	\node [draw, blue] at (m) [circle through={(iA)}] {};
\end{scope}

\node[below of = a] (at){$c_A$};
\draw (a) circle (3);
\node at (a) {A};

\node[below of = b] (bt){$c_B$};
\draw (b) circle (3);
\node at (b) {B};

\node[below of = m] (mt){$m$};


\begin{scope}
	
	\clip (b) circle (3);
	\node [draw, blue] at (iA) [circle through={(iB)}] {};
\end{scope}

\begin{scope}
	\clip (a) circle (3);
	\node [draw, green] at (iB) [circle through={(iA)}] {};
\end{scope}


%\draw[<->] (3,0) to node[below] {$c = s\cdot r$} (6,0);
%\draw[<->] (0) to node[below] {$t'c = r+(s\cdot r)$} (3);
%\draw[-, draw=blue] (0) ++(0,0.2) to node[above] {$|r_{i-1}r_{i}|$} ++(9,0);
%\draw[-] (0)  to node[below] {$|r_{i-1}p_{i-1}|$} (1);
%\draw[-] (2)  to node[below] {$|p_{i}r_{i}|$} (3);
\end{tikzpicture}
\end{document}
```

---


- $c_A$ is the center point of $A$ 
- $c_B$ is the center point of $B$
- $m$ is the center point between $c_A$ and $c_B$ $m = \frac {c_A +c_B}2$
- $up_A$ as is subset of A, which is *above* the line $c_A,c_B$. 
- $down_A$ is the subset of $A$ that is below the same line.
- $u_A = min(\{|pm| |p \in up_A\})$
- $d_A = min(\{|pm| |p \in down_A\})$



