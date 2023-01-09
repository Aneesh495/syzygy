# Mathematical Theory of Syzygy

## 1. Commutative Algebra and Groebner Bases

### 1.1 The Ideal Membership Problem
Let $R = k[x_0, \dots, x_{n-1}]$ be a multivariate polynomial ring over a field $k$. Given an ideal $I = \langle f_1, \dots, f_s \rangle$ and a test polynomial $f \in R$, the ideal membership problem asks whether $f \in I$, or equivalently whether $f$ vanishes identically on the affine algebraic variety:
$$ V(I) = \{ x \in k^n \mid f_i(x) = 0 \text{ for all } i = 1, \dots, s \}. $$

By Hilbert's Basis Theorem (1890), every ideal in $R$ is finitely generated. However, the classical Euclidean division algorithm fails for multivariate polynomials because leading terms depend on the chosen monomial ordering $\succ$.

### 1.2 Monomial Orders
A monomial ordering $\succ$ on the set of monomials $\{x^\alpha \mid \alpha \in \mathbb{N}^n\}$ is a well-ordering satisfying:
1. Total order: for any $\alpha \ne \beta$, either $x^\alpha \succ x^\beta$ or $x^\beta \succ x^\alpha$.
2. Compatibility with multiplication: $x^\alpha \succ x^\beta \implies x^{\alpha + \gamma} \succ x^{\beta + \gamma}$ for all $\gamma \in \mathbb{N}^n$.
3. Well-ordering: every non-empty set of monomials contains a least element with respect to $\succ$.

Key orders implemented in Syzygy:
- **Lexicographic (Lex)**: $x^\alpha \succ_{\text{lex}} x^\beta$ iff the leftmost non-zero entry of $\alpha - \beta$ is positive. Eliminates variables from left to right.
- **Degree Reverse Lexicographic (DegRevLex / Grevlex)**: $x^\alpha \succ_{\text{grevlex}} x^\beta$ iff $|\alpha| > |\beta|$, or $|\alpha| = |\beta|$ and the rightmost non-zero entry of $\alpha - \beta$ is negative. Highly efficient for computation and degree-compatible resolutions.
- **Elimination Block Orders**: Partitions variables into disjoint subsets $X$ and $Y$ with $X \succ Y$, guaranteeing that any polynomial in $I \cap k[Y]$ has leading term depending only on $Y$.

### 1.3 Groebner Bases and Buchberger's Criterion
A finite subset $G = \{g_1, \dots, g_t\} \subset I$ is a Groebner basis for $I$ with respect to $\succ$ if:
$$ \langle \text{LT}(g_1), \dots, \text{LT}(g_t) \rangle = \langle \text{LT}(I) \rangle. $$

**Buchberger's Criterion (1965)**:
A generating set $G$ is a Groebner basis iff for every pair $g_i, g_j \in G$, the S-polynomial:
$$ S(g_i, g_j) = \frac{\text{LCM}(\text{LM}(g_i), \text{LM}(g_j))}{\text{LT}(g_i)} g_i - \frac{\text{LCM}(\text{LM}(g_i), \text{LM}(g_j))}{\text{LT}(g_j)} g_j $$
reduces to zero modulo $G$:
$$ S(g_i, g_j) \xrightarrow{G} 0. $$

### 1.4 Gebauer-Moeller Criteria
To avoid computing S-polynomials known a priori to reduce to zero:
1. **Product Criterion**: If $\text{LCM}(\text{LM}(g_i), \text{LM}(g_j)) = \text{LM}(g_i) \cdot \text{LM}(g_j)$ (coprime leading monomials), then $S(g_i, g_j) \xrightarrow{G} 0$.
2. **Chain Criterion (Triangle Criterion)**: If there exists $g_k \in G$ such that $\text{LM}(g_k)$ divides $\text{LCM}(\text{LM}(g_i), \text{LM}(g_j))$, and both $(i, k)$ and $(k, j)$ have already been processed, then the pair $(i, j)$ is redundant.

### 1.5 Faugere's F4 Algorithm
Instead of reducing critical pairs one by one as in classic Buchberger, F4 (Faugere 1999) selects a batch of critical pairs of minimal degree, performs symbolic preprocessing to collect all necessary reductor polynomials, and builds a large sparse Macaulay matrix whose rows correspond to polynomials and columns to monomials sorted by $\succ$. Simultaneous Gaussian elimination over the coefficient field $k$ reduces all pairs in parallel.

### 1.6 Faugere's Signature-Based F5 Algorithm
Faugere's F5 algorithm (2002) equips every polynomial with a signature $(i, m)$, tracking its origin in the free module $R^s = \bigoplus_{i=1}^s R e_i$. By enforcing the **F5 Syzygy Criterion** (discarding signatures divisible by leading terms of syzygies) and the **Rewritten Criterion**, F5 eliminates all reductions to zero when the input polynomials form a regular sequence.

---

## 2. Graded Free Resolutions and Syzygies

### 2.1 Hilbert's Syzygy Theorem (1890)
Let $R = k[x_0, \dots, x_{n-1}]$ be the polynomial ring with standard grading $\deg(x_i) = 1$. Let $M$ be a finitely generated graded $R$-module (such as the quotient ring $R/I$).

A graded free resolution of $M$ is an exact sequence of graded free $R$-modules:
$$ 0 \longrightarrow F_p \xrightarrow{d_p} F_{p-1} \xrightarrow{d_{p-1}} \dots \xrightarrow{d_2} F_1 \xrightarrow{d_1} F_0 \xrightarrow{d_0} M \longrightarrow 0 $$
where each $F_i = \bigoplus_{j} R(-j)^{\beta_{i, j}}$, and each differential $d_i$ is a degree-preserving homogeneous homomorphism ($d^2 = 0$).

**Hilbert's Syzygy Theorem**:
The projective dimension of $M$ is bounded by the number of variables:
$$ \text{pd}(M) \le n. $$
Every graded free resolution over $k[x_0, \dots, x_{n-1}]$ terminates in at most $n$ homological steps.

### 2.2 Schreyer's Theorem (1980)
Schreyer's Theorem provides a constructive algorithm to compute syzygies. Given a Groebner basis $G = \{g_1, \dots, g_s\}$ of a submodule $M \subseteq R^m$:
1. Define the Schreyer module term order $\succ_{\text{Schreyer}}$ on $R^s$ induced by $G$:
   $$ x^\alpha e_i \succ_{\text{Schreyer}} x^\beta e_j \iff \text{LT}(x^\alpha g_i) \succ \text{LT}(x^\beta g_j) \text{ (breaking ties with } i < j \text{)}. $$
2. For each critical pair $S(g_i, g_j) = u_{ij} g_i - u_{ji} g_j = \sum_{k} q_k g_k$, the syzygy vector:
   $$ \sigma_{ij} = u_{ij} e_i - u_{ji} e_j - \sum_k q_k e_k \in R^s $$
   satisfies $\phi(\sigma_{ij}) = 0$.
3. The set of all $\sigma_{ij}$ forms a Groebner basis for the syzygy module $\ker(\phi)$ with respect to $\succ_{\text{Schreyer}}$!
4. Iterating Schreyer's construction yields the complete free resolution.

### 2.3 Minimization and Graded Betti Numbers
A free resolution is **minimal** if the matrices representing the differentials $d_k$ contain no non-zero constant entries (i.e. $d_k(F_k) \subseteq \mathfrak{m} F_{k-1}$ where $\mathfrak{m} = \langle x_0, \dots, x_{n-1} \rangle$).

By Nakayama's Lemma, every finitely generated graded module has an essentially unique minimal free resolution. The ranks of the free modules in the minimal resolution define the **graded Betti numbers**:
$$ \beta_{i, j}(M) = \dim_k \text{Tor}_i^R(M, k)_j. $$

In standard Macaulay2 formatting, the entry in column $i$ and row $(j - i)$ represents $\beta_{i, j}$:
- Total Betti numbers: $\beta_i = \sum_j \beta_{i, j}$.
- Projective dimension: $\text{pd}(M) = \max \{ i \mid \beta_i \ne 0 \}$.
- Castelnuovo-Mumford regularity: $\text{reg}(M) = \max \{ j - i \mid \beta_{i, j} \ne 0 \}$.

---

## 3. Real Algebraic Geometry and Solvers

### 3.1 Zero-Dimensional Varieties and the Stetter-Moeller Eigenvalue Method
When an ideal $I$ is zero-dimensional ($\dim V(I) = 0$), the quotient ring:
$$ A = k[x_0, \dots, x_{n-1}] / I $$
is a finite-dimensional vector space over $k$ of dimension $D = \dim_k(A)$. A monomial basis for $A$ is given by the standard monomials:
$$ \mathcal{B} = \{ x^\alpha \mid x^\alpha \notin \langle \text{LT}(I) \rangle \}. $$

For any coordinate variable $x_i$, the multiplication operator:
$$ M_{x_i} : A \to A, \quad [f] \mapsto [x_i \cdot f] $$
can be represented as a $D \times D$ matrix with respect to $\mathcal{B}$.

**Stetter-Moeller Theorem (1993)**:
1. The multiplication operators $M_{x_0}, \dots, M_{x_{n-1}}$ mutually commute: $M_{x_i} M_{x_j} = M_{x_j} M_{x_i}$.
2. The eigenvalues of $M_{x_i}$ are precisely the $x_i$-coordinates of the points in $V(I)$.
3. The common eigenvectors of the commuting matrices $\{M_{x_i}\}$ correspond to the evaluation homomorphisms at the points of $V(I)$.

### 3.2 Sturm's Theorem and Real Root Isolation
Let $f(x) \in \mathbb{R}[x]$ be a square-free univariate polynomial. The Sturm sequence $f_0, f_1, \dots, f_m$ is defined by:
$$ f_0 = f, \quad f_1 = f', \quad f_{k+1} = -\text{rem}(f_{k-1}, f_k). $$
Let $V(c)$ denote the number of sign variations in the sequence $[f_0(c), f_1(c), \dots, f_m(c)]$, ignoring zeros.

**Sturm's Theorem (1829)**:
The number of distinct real roots of $f$ in the half-open interval $(a, b]$ is exactly:
$$ N(a, b] = V(a) - V(b). $$

By combining Sturm sequences with Cauchy root bounds $|z| \le 1 + \max_{i < n} |a_i / a_n|$, Syzygy isolates all real roots into disjoint rational isolating intervals with exact multiplicity determinations.

### 3.3 Cylindrical Algebraic Decomposition (CAD)
Introduced by George Collins (1975), CAD decomposes $\mathbb{R}^n$ into a finite collection of disjoint connected semialgebraic cells on each of which a given set of polynomials has constant sign $(-, 0, +)$.
1. **Projection Phase**: Computes the Collins projection operator $\text{Proj}(A) \subset \mathbb{R}[x_1, \dots, x_{k-1}]$ using principal subresultant coefficients and discriminants.
2. **Base Phase**: Solves the univariate polynomials in $\mathbb{R}^1$ using Sturm isolation.
3. **Lifting Phase**: Erects cylinders over each cell and isolates real roots of the polynomials evaluated at sample points to subdivide cylinders into sections and sectors.

---

## 4. Invariant Theory of Finite Groups

Let $G \subset \text{GL}(n, k)$ be a finite matrix group acting linearly on $R = k[x_0, \dots, x_{n-1}]$ by $(g \cdot f)(x) = f(g^T x)$. The invariant subring is:
$$ R^G = \{ f \in R \mid g \cdot f = f \text{ for all } g \in G \}. $$

### 4.1 The Reynolds Operator
When $\text{char}(k) = 0$ or $\text{char}(k) \nmid |G|$, the Reynolds operator:
$$ R_G : R \to R^G, \quad R_G(f) = \frac{1}{|G|} \sum_{g \in G} g \cdot f $$
is a $k$-linear projection from $R$ onto $R^G$ satisfying:
- $R_G(f) = f$ for all $f \in R^G$.
- $R_G(f \cdot h) = f \cdot R_G(h)$ for all $f \in R^G$ and $h \in R$ (Reynolds identity).

### 4.2 Molien's Theorem (1897)
The graded dimensions of the invariant subring $R^G = \bigoplus_{d=0}^\infty R^G_d$ are encoded by the Molien series (Hilbert-Poincare series):
$$ \Phi(t) = \sum_{d=0}^\infty \dim_k(R^G_d) \, t^d = \frac{1}{|G|} \sum_{g \in G} \frac{1}{\det(I - t \cdot g)}. $$

By the **Hochster-Eagon Theorem (1971)**, $R^G$ is a Cohen-Macaulay graded ring, admitting a Hironaka decomposition:
$$ R^G = \bigoplus_{j=1}^m \eta_j \cdot k[\theta_1, \dots, \theta_n] $$
into algebraically independent primary invariants $\theta_i$ and secondary invariants $\eta_j$.

### 4.3 Noether's Degree Bound (1916)
In characteristic zero, $R^G$ is generated as an algebra by homogeneous invariants of degree at most $|G|$:
$$ \deg(f_i) \le |G|. $$

---

## 5. Tropical Geometry

Tropical geometry is algebraic geometry over the min-plus semiring $\mathbb{T} = (\mathbb{R} \cup \{+\infty\}, \oplus, \odot)$:
$$ a \oplus b = \min(a, b), \quad a \odot b = a + b. $$

### 5.1 Tropical Hypersurfaces
A tropical polynomial in $n$ variables $w \in \mathbb{R}^n$ is a concave piecewise-linear function:
$$ F(w) = \min_{\alpha \in A} \left( c_\alpha + \sum_{i=0}^{n-1} \alpha_i w_i \right). $$

The tropical hypersurface (corner locus) $\mathcal{T}(F) \subset \mathbb{R}^n$ is the non-differentiability locus where the minimum is achieved by at least two distinct terms:
$$ \mathcal{T}(F) = \{ w \in \mathbb{R}^n \mid \exists \alpha \ne \beta \text{ with } c_\alpha + \langle \alpha, w \rangle = c_\beta + \langle \beta, w \rangle = F(w) \}. $$

### 5.2 Kapranov's Theorem and Newton Polytope Duality
**Kapranov's Theorem (2000)**:
For a polynomial $f(x)$ over a valued field $K$ (e.g. Puiseux series with $t$-adic valuation), the tropical hypersurface $\mathcal{T}(\text{trop}(f))$ equals the Euclidean closure of the image under coordinate valuation of the classical variety:
$$ \mathcal{T}(\text{trop}(f)) = \overline{\text{val}(V(f) \cap (K^*)^n)}. $$

The polyhedral fan of $\mathcal{T}(F)$ is dual to the regular polyhedral subdivision of the **Newton polytope** $\text{Newt}(f) = \text{Conv}(\{\alpha \mid c_\alpha \ne 0\})$.

---

## 6. Cellular Sheaf Cohomology and Hodge Laplacians

### 6.1 Cellular Sheaves over Cell Complexes
Let $X$ be a regular cell complex (Alexandrov poset $(X, \trianglelefteq)$). A cellular sheaf $\mathcal{F}$ over $X$ with values in $\text{Vect}_k$ assigns:
1. To each cell $\sigma \in X$, a stalk vector space $\mathcal{F}(\sigma) \cong k^{d_\sigma}$.
2. To each incidence $\sigma \trianglelefteq \tau$ (where $\sigma$ is a facet of $\tau$), a linear restriction map:
   $$ \rho_{\sigma \trianglelefteq \tau} : \mathcal{F}(\sigma) \to \mathcal{F}(\tau). $$
3. Commutativity: for any codimension-2 incidence $\sigma \trianglelefteq \gamma$ with intermediate facets $\tau_1, \tau_2$:
   $$ \rho_{\tau_1 \trianglelefteq \gamma} \circ \rho_{\sigma \trianglelefteq \tau_1} = \rho_{\tau_2 \trianglelefteq \gamma} \circ \rho_{\sigma \trianglelefteq \tau_2}. $$

### 6.2 Sheaf Cohomology
The $k$-cochains are $C^k(X; \mathcal{F}) = \bigoplus_{\sigma \in X_k} \mathcal{F}(\sigma)$.
The coboundary operator $\delta^k : C^k(X; \mathcal{F}) \to C^{k+1}(X; \mathcal{F})$ is defined on a cochain $x \in C^k$ evaluated at a cell $\tau \in X_{k+1}$ by:
$$ (\delta^k x)_\tau = \sum_{\sigma \trianglelefteq \tau} [\sigma : \tau] \, \rho_{\sigma \trianglelefteq \tau}(x_\sigma) $$
where $[\sigma : \tau] \in \{-1, +1\}$ is the topological incidence number.
Because $\delta^{k+1} \circ \delta^k = 0$, the $k$-th sheaf cohomology group is:
$$ H^k(X; \mathcal{F}) = \frac{\ker(\delta^k)}{\text{im}(\delta^{k-1})}. $$
The global sections $H^0(X; \mathcal{F}) = \ker(\delta^0)$ represent consistent local assignments over the entire complex.

### 6.3 Sheaf Hodge Laplacian and Hodge Decomposition
The $k$-th Sheaf Hodge Laplacian $L_k : C^k(X; \mathcal{F}) \to C^k(X; \mathcal{F})$ is the self-adjoint linear operator:
$$ L_k = \delta_{k-1} \delta_{k-1}^* + \delta_k^* \delta_k = L_k^{\text{down}} + L_k^{\text{up}}. $$

**Hodge Decomposition Theorem for Cellular Sheaves (Hansen and Ghrist 2019)**:
1. Orthogonal decomposition of cochains:
   $$ C^k(X; \mathcal{F}) = \text{im}(\delta_{k-1}) \oplus \ker(L_k) \oplus \text{im}(\delta_k^*). $$
2. The space of harmonic cochains is canonically isomorphic to sheaf cohomology:
   $$ \mathcal{H}^k(X; \mathcal{F}) = \ker(L_k) \cong H^k(X; \mathcal{F}). $$
3. In particular:
   $$ \dim \ker(L_k) = b_k(X; \mathcal{F}). $$
For $k = 0$, $L_0 = \delta_0^* \delta_0$ generalizes the graph Laplacian matrix, with $\ker(L_0)$ spanning the global sections of $\mathcal{F}$.
