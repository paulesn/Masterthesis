

Given a set of Points $S: ((x,y,z,b)| x,y,z \in \mathbb R, b \in (1,4))$, a start $s \in S$, a target $t \in S$ and a jump range $r \in \mathbb R$, we want to know the shortest path from $s$ to $t$.
A path of $\pi = (\pi_1, \pi_2, \dots, \pi_n)$ is an ordered list of points out of $S$, that holds for each element:
$$
\sqrt{(x_i-x_{i+1})^2+(y_i-y_{i+1})^2+(z_i-z_{i+1})^2}\leq r\cdot b_i
$$
The length of a graph is defined by the amount of points in the list.