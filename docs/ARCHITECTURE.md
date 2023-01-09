# Syzygy System Architecture

## 1. System Overview

Syzygy is an exact computational algebraic geometry, commutative algebra, and cellular sheaf engine implemented in Haskell with SIMD-accelerated C kernels and an interactive WebGL variety workbench.

The system is structured as a layered pipeline from arbitrary-precision exact rings up to homological resolutions, algebraic variety solvers, and cellular sheaf Hodge theory.

```
       +-------------------------------------------------------+
       |             User Interface & Interaction              |
       |  REPL (CLI.hs)  |  ASCII Raymarcher  |  WebGL Studio  |
       +-------------------------------------------------------+
                                  |
       +-------------------------------------------------------+
       |             Homological & Geometric Engines           |
       |  Free Resolutions | Betti Tables | Hilbert Series     |
       |  Elimination      | Radical      | Primary Decomp     |
       |  Stetter Solver   | CAD          | Invariant Theory   |
       +-------------------------------------------------------+
                                  |
       +-------------------------------------------------------+
       |              Groebner Basis Core Engine               |
       |  Buchberger | Faugere F4 | Signature F5 | FGLM Basis  |
       |  Gebauer-Moeller Criteria | S-Polynomial Reductions   |
       +-------------------------------------------------------+
                                  |
       +-------------------------------------------------------+
       |              Graded Free Module Algebra               |
       |  Module Monomials | TOP / POT / Schreyer Orderings    |
       |  Module Division  | Submodule S-Vectors & Syzygies    |
       +-------------------------------------------------------+
                                  |
       +-------------------------------------------------------+
       |               Exact Algebraic Substrate               |
       |  Polynomials | Multivariate Division | GCD | Parser   |
       |  Monomials   | Lex / DegLex / Grevlex / Block Orders  |
       |  Exact Rings | BigInt | Rational | FiniteField | Dual |
       |  Linear Alg  | Matrix | Gaussian RREF | Smith Form    |
       +-------------------------------------------------------+
                                  |
       +-------------------------------------------------------+
       |             C Acceleration Subsystem (SIMD)           |
       |  64-bit Packed Exponents | F4 Sparse Row Elimination  |
       |  Sheaf Hodge Laplacian SpMV Spectral Iteration        |
       +-------------------------------------------------------+
```

---

## 2. Layered Module Hierarchy

### 2.1 Layer 0: Exact Algebraic Arithmetic (`Syzygy.Core.*`)
- `BigInt`: Arbitrary-precision integer arithmetic, extended Euclidean algorithm, Miller-Rabin probabilistic primality, and modular exponentiation.
- `Rational`: Exact rational numbers with canonical gcd normalization and zero-denominator invariants.
- `FiniteField`: Galois field Z/pZ arithmetic with modular inverse via Extended Euclidean Algorithm and batch inversion.
- `GaloisField`: Prime power field extensions GF(p^k) represented modulo irreducible polynomials.
- `AlgebraicNumber`: Real algebraic numbers encoded as isolating intervals with Thom's lemma Sturm sign determinations.
- `Interval`: Rigorous interval arithmetic with interval Newton-Raphson refinement.
- `Dual`: Dual numbers and hyper-dual numbers for forward-mode algorithmic differentiation without truncation error.

### 2.2 Layer 1: Monomial Representations and Orders (`Syzygy.Monomial.*`)
- `Order`: Term orders including Lexicographical (Lex), Degree Lexicographical (DegLex), Degree Reverse Lexicographical (DegRevLex), Block Orders (Elimination), and Matrix Orders (Robbiano's Theorem).
- `Packed`: High-performance 64-bit bit-packed exponent representations (8 variables at 8 bits each) allowing single-instruction comparison and SIMD divisibility checks.
- `Sparse`: Arbitrary-variable sparse exponent vectors backed by Data.Map.Strict.
- `Operations`: Least Common Multiples (LCM), Greatest Common Divisors (GCD), and coprime monomial detection.

### 2.3 Layer 2: Multivariate Polynomial Algebra (`Syzygy.Polynomial.*`)
- `Base`: Canonical polynomial representation as ordered lists of non-zero terms paired with a monomial ordering: `Poly [Term c] MonomialOrder`.
- `Arithmetic`: Addition, subtraction, scalar multiplication, term multiplication, and binary exponentiation.
- `Division`: Multivariate polynomial division algorithm with remainder with respect to ordered lists of divisor polynomials.
- `Univariate`: Specialized univariate polynomial algebra, polynomial square-free factorization, and remainder sequences.
- `MultivariateGCD`: Multivariate GCD computation via subresultant pseudo-division and variable reduction.
- `Homogenization`: Homogenization with respect to an auxiliary variable and standard affine dehomogenization.
- `Parser` and `Pretty`: Recursive-descent tokenization and formatting (ASCII, Unicode superscripts, LaTeX math).

### 2.4 Layer 3: Exact Linear Algebra (`Syzygy.Linear.*`)
- `Vector`: Dense vectors over arbitrary fields.
- `Matrix`: Dense rectangular matrices, transposition, block slicing, matrix-vector multiplication, and trace.
- `SparseMatrix`: Compressed Sparse Row (CSR) matrix representation for large sparse Macaulay matrices.
- `Gaussian`: Gauss-Jordan elimination, Reduced Row Echelon Form (RREF), matrix rank, null space bases, and Bareiss fraction-free determinants.
- `SmithNormalForm`: Diagonalization of integer matrices U * A * V = D for finitely generated abelian groups and integral homology.

### 2.5 Layer 4: Groebner Bases (`Syzygy.Grobner.*`)
- `SPolynomial`: Critical pairs, S-polynomial formation, and leading term cancellation.
- `Criteria`: Gebauer-Moeller product and chain criteria for early detection and elimination of redundant critical pairs.
- `Buchberger`: Classic Buchberger algorithm with normal form reduction.
- `F4`: Faugere's F4 linear algebra algorithm, constructing Macaulay matrices for batches of critical pairs and performing simultaneous reduction via Gaussian elimination.
- `F5`: Faugere's signature-based F5 algorithm avoiding all reductions to zero for regular sequences.
- `Reduce`: Minimal and unique canonical reduced Groebner basis computation.
- `HilbertDriven`: Hilbert-driven Groebner basis algorithm bounding the highest degree necessary for S-polynomial exploration.
- `FGLM`: FGLM basis conversion algorithm transforming zero-dimensional Groebner bases from Grevlex to Lex in O(n * D^3) time.

### 2.6 Layer 5: Graded Free Resolutions & Homological Syzygies (`Syzygy.Resolution.*`)
- `Module`: Graded free modules R^m, module elements, and degree shift assignments.
- `ModuleOrder`: Module term orders (Term Over Position TOP, Position Over Term POT, and Schreyer term orders).
- `ModuleGrobner`: Submodule division, module S-polynomials, and submodule Groebner bases.
- `Syzygy`: Computation of syzygy modules ker(phi) using Schreyer's Theorem.
- `FreeResolution`: Construction of finite graded chain complexes 0 -> F_p -> ... -> F_0 -> R/I -> 0.
- `Minimization`: Minimization of free resolutions via row/column elimination over the residue field k.
- `Betti`: Computation of invariant graded Betti tables beta_{i, j} = dim Tor_i(M, k)_j and Macaulay2 formatting.
- `Regularity`: Castelnuovo-Mumford regularity, projective dimension, depth, and Cohen-Macaulay checks.
- `HilbertSeries`: Hilbert-Poincare series rational functions P(t) / (1 - t)^d and variety multiplicity (degree).

### 2.7 Layer 6: Algebraic Geometry & Real Algebraic Geometry (`Syzygy.Geometry.*` & `Syzygy.Solver.*`)
- `Ideal`: Ideal intersection, sum, product, and equality via elimination.
- `Elimination`: Variable elimination via block orders and projection of algebraic varieties.
- `Implicitization`: Polynomial and rational parametric hypersurface implicitization using the Rabinowitsch trick.
- `Colon`: Colon ideals (I : J) and ideal saturation (I : J^infinity).
- `Radical`: Radical ideal membership testing via Rabinowitsch trick and Seidenberg's univariate projection algorithm.
- `Dimension`: Affine Krull dimension via independent variable subsets modulo leading ideals.
- `PrimaryDecomposition`: Zero-dimensional primary decomposition using the Gianni-Trager-Zacharias (GTZ) algorithm.
- `Singularities`: Jacobian matrix, Jacobian criterion, singular locus ideal, and algebraic variety smoothness checks.
- `Quotient`: Quotient ring structure k[x_0, ..., x_{n-1}] / I and monomial vector space bases.
- `Multiplication`: Coordinate multiplication operators M_{x_i} and Faddeev-LeVerrier characteristic polynomials.
- `Stetter`: Stetter-Moeller eigenvalue polynomial solver for zero-dimensional varieties.
- `RealRoot`: Sturm sequences, sign variation counting, Cauchy bounds, and exact real root isolation.
- `SturmHabicht`: Sturm-Habicht subresultant polynomials and Tarski queries.
- `CAD`: Cylindrical Algebraic Decomposition with Collins projection operator and cell sample evaluation.
- `SOS`: Sum of Squares polynomial decomposition and Gram matrix semi-definite representation.

### 2.8 Layer 7: Invariant Theory & Tropical Geometry (`Syzygy.Invariant.*` & `Syzygy.Tropical.*`)
- `FiniteGroup`: Representation of finite matrix groups in GL(n, k) (Cyclic, Dihedral, Symmetric, Klein 4-group).
- `Reynolds`: Reynolds projection operator R_G and polynomial group action.
- `Molien`: Molien's formula expanding the Hilbert-Poincare series 1 / |G| * sum 1 / det(I - t*g).
- `Generators`: Minimal algebra generators of invariant subrings k[x]^G and Noether's degree bound.
- `Semiring`: Min-plus and max-plus idempotent semirings for tropical geometry.
- `Hypersurface`: Tropical polynomials, corner locus hypersurfaces, and Newton polygon convex hulls.

### 2.9 Layer 8: Cellular Sheaves & Spectral Hodge Theory (`Syzygy.Sheaf.*`)
- `Complex`: Regular cell complexes, simplicial triangulations, and boundary matrices satisfying d^2 = 0.
- `CellularSheaf`: Stalk vector spaces F(sigma) and linear restriction maps rho_{sigma <| tau}.
- `Cohomology`: Coboundary operator delta^k satisfying delta^2 = 0, sheaf cohomology groups H^k(X; F), and sheaf Betti numbers.
- `Laplacian`: Sheaf Hodge Laplacians L_k = delta_{k-1} delta_{k-1}^* + delta_k^* delta_k, harmonic cochain bases, and heat diffusion.

---

## 3. Data Structures & Memory Layout

### 3.1 Polynomial Representation
Polynomials are stored in normal form as strictly sorted, strictly non-zero terms:
```haskell
data Poly c = Poly ![Term c] !MonomialOrder
data Term c = Term !c !SparseMonomial
newtype SparseMonomial = SparseMonomial { unSparseMonomial :: Map.Map Int Int }
```
- Strict spine (`!`) prevents memory leaks from unreduced thunks during iterative polynomial division.
- Terms are kept sorted in descending order according to the embedded `MonomialOrder`.
- Identity monomial 1 is represented as the empty map `Map.empty`.

### 3.2 Chain Complex Representation
```haskell
data ChainComplex c = ChainComplex
  { ccRanks         :: ![Int]         {- Ranks of free modules F_0, F_1, ..., F_p -}
  , ccShifts        :: ![[Int]]       {- Shift degrees for basis elements -}
  , ccDifferentials :: ![[[Poly c]]]  {- Differential matrices d_1, ..., d_p -}
  }
```
Each differential `d_k` is an r_{k-1} x r_k matrix of polynomials satisfying `d_{k-1} . d_k = 0`.
During minimization, non-zero constant entries are eliminated via elementary row and column operations over the residue field k.
