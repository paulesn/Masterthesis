> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=27&selection=83,0,165,1&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.9]]
> > **Definition 1.2.1 (Spanner)**. Let S be a set of n points in Rd and let t ≥ 1 be a real number. A t-spanner for S is an undirected graph G with vertex set S, such that for any two points p and q of S, there is a path in G between p and q, whose length is less than or equal to t|pq|. Any path satisfying this condition is called a t-spanner path between p and q.

> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=127&selection=17,0,165,1&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.109]]
> > **Definition 6.1.1 (Gap Property)**. Let $w \geq 0$ be a real number, and let $E$ be a set of directed edges in Rd . 
> > 1. We say that E satisfies the $w$-gap property if for any two distinct edges $(p, q)$ and $(r, s)$ in $E$, we have $|pr| > w \cdot \min(|pq|, |rs|)$. 
> > 2. We say that $E$ satisfies the strong $w$-gap property if for any two distinct edges $(p, q)$ and $(r, s)$ in $E$, we have $|pr| > w \cdot \min(|pq|, |rs|)$ and $|qs| > w \cdot \min(|pq|, |rs|)$. The Gap Theorem below bounds the total length of any set of edges that satisfies the gap property. Recall that for any directed edge $(p, q)$, $p$ is called the source, and $q$ is called the sink.
> 
> That is relevant to the Leapfrog Property

> [!PDF|red] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=29&selection=21,0,71,9&color=red|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.11]]
> > **Basic spanner problem**: Let S be a set of n points in Rd , and let t > 1 be a real number. Does there exist a t-spanner for S having at most ctd n edges, where ctd is a real number that depends only on t and d? If so, how much time does it take to compute such a t-spanner?

> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=170&selection=209,0,291,1&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.152]]
> > **Lemma 9.1.2.** Let $s > 0$ be a real number, let $A$ and $B$ be two finite sets of points that are well-separated with respect to $s$, let $p$ and $p′$ be any two points in $A$, and let $q$ and $q′$ be any two points in $B$. Then 
> >  3. $|pp′ | \leq (2/s)|pq|$, and 
> >  4. $|p′q′ | \leq (1+\frac4s)|pq|$.

> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=172&selection=68,0,97,3&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.154]]
> > **Basic spanner construction**: Construct a well-separated pair decomposition with separation ratio $s > 4$, and take one (arbitrary) edge for each pair of the decomposition. This results in a t-spanner with $t = \frac{(s + 4)}{(s − 4)}$.
> 
> This is relevant for points in any $\mathbb R^d$
> 
 
> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=186&selection=263,0,346,7&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.168]]
> > **Corollary 9.4.7 (WSPD-Spanner)**. Let $S$ be a set of n points in $\mathbb R^d$ and let $t > 1$ be a real number. In $\mathcal O(n \log n + \frac n{(t − 1)d} )$ time, we can construct a $t$-spanner for $S$ having $O(\frac n{(t − 1)d })$ edges. Thus, the WSPD-spanner has $\mathcal O(n)$ edges and can be constructed in $\mathcal O(n \log n)$ time.
> 
> This also holds for any d
> 
> 1. There is no nontrivial bound on the degree of the WSPD-spanner. However, we will show in Section 10.1.1 that the technique of Section 5.5.3 can be used to transform the WSPD-spanner to a spanner of bounded degree. 
> 2. There is no nontrivial bound on the spanner diameter of the WSPD-spanner. However, we will show in Section 10.2 that for a special choice of the representatives for each well-separated pair $\{A_i , B_i \}$, the WSPD-spanner has spanner diameter $\mathcal O(\log n)$. 
> 3. By combining the Gap Theorem (Theorem 6.1.2) of Chapter 6 with the so-called dumbbell trees that will be introduced in Chapter 11, it can be shown that the weight of the WSPD-spanner is $\mathcal O(\log n)$ times the weight of a minimum spanning tree of the point set S; see Exercises 9.12 and 11.6.

> [!PDF|important] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=191&selection=6,0,170,2&color=important|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.173]]
> > A metric space is a pair (S, δ), where S is a (finite or infinite) set, whose elements are called points, and δ : S × S −→ R is a function that assigns a distance δ(p, q) to any two points p and q in S, and that satisfies the following three conditions: 1. For all points p and q in S, δ(p, q) ≥ 0. 2. For all points p and q in S, δ(p, q) = 0 if and only if p = q. 3. For all points p and q in S, δ(p, q) = δ(q, p). 4. For all points p, q, and r in S, δ(p, q) ≤ δ(p, r) + δ(r, q).
> 
> I have to see if they can proof that the WSP spanner holds for any metric but I think the amount of hops is a valid metric and could be used here.

> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=215&selection=91,0,169,45&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.197]]
> > Length-grouping property: Either ′ ≤ 2 or ′ ≥ /β, where β is a fixed real number such that 0 < β < 1/2. That is, D and D′ either have approximately the same lengths or their lengths differ by a large amount. Empty-region property: If ′ ≤ 2 , then the distance between any head of D and any head of D′ is larger than γ , where γ is a fixed positive real number. That is, any head of D and any head of D′ are disjoint and separated by a large amount.
> 
> The dumbell theory could be used for some proofs with the maximum jump range if $\beta$ and $\gamma$ are choosen carefully

> [!PDF|note] [[Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007).pdf#page=277&selection=17,0,163,13&color=note|Giri Narasimhan, Michiel Smid - Geometric Spanner Networks-Cambridge University Press (2007), p.259]]
> > Definition 14.1.1 (Leapfrog Property). Let t > 1 be a real number. A set E of undirected edges in Rd is said to satisfy the t-leapfrog property, if for every k ≥ 2, and for every sequence {p1, q1}, {p2, q2}, . . . , {pk , qk } of k pairwise distinct edges of E, t|p1q1| < k∑ i=2 |p i qi | + t ( |p1p2| + k−1∑ i=2 |qi p i+1| + |qk q1| ) . Observe that this definition requires that the inequality holds for every permutation of the k edges, and for every labeling of their endpoints. In other words, for every set of k pairwise distinct edges, the leapfrog property gives k! 2k inequalities.


