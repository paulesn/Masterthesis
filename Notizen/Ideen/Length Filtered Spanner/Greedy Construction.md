
> [!tip] Lemma
> Given a Graph $G = (V,E)$, a $t$-Spanner $S = (V,E_S)$ for $G$ constructed with the greedy construction method and a threshold  $c$. The subgraph $S' = (V, \{e \in E_S | d(e)\leq c\})$ contains a path from $a \in V$ to $b \in V$ iff the  Graph $G' = (V_S, \{e \in E | d(e)\leq c\})$ contains a path from $a$ to $b$

Proof:
$\exists (a \to b) \in S' \implies \exists (a \to b) \in G'$:
As $S$ is a subgraph of $G$ [^1], and $G'$ removes at least all edges that are also removed in $S'$, every edge that is in $S'$ is also in $G'$. Therefore every path in $S'$ is also in $G'$.

$\exists (a \to b) \in G' \implies \exists (a \to b) \in S'$:
We have a path from $a$ to $b$ in $G'$. That means that there is a path from $a$ to $b$ in $G$. Which means that there are edges shorter than $c$ that can be added to create a path $(a \to b)$ before any edge longer than $c$ will be added. Therefore there is such a path in $S'$.