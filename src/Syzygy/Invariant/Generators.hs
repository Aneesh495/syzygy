{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Invariant.Generators
Description : Minimal algebra generators, Noether's bound, and subring membership
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Computes the fundamental invariants (algebra generators) of the invariant ring
k[x_0, ..., x_{n-1}]^G for a finite matrix group G in GL(n, k).

By Noether's Theorem (1916), in characteristic zero, the invariant ring is
generated as an algebra by homogeneous invariants of degree at most |G|:
  deg(f_i) <= |G|

The generators are computed degree by degree:
1. Generate all monomials of degree d.
2. Apply the Reynolds operator R_G to each monomial.
3. Extract a linearly independent basis of invariants of degree d via Gaussian elimination.
4. Filter out invariants expressible as polynomials in previously found lower-degree generators.
5. Invariants that cannot be expressed in terms of lower-degree generators are minimal generators.
-}
module Syzygy.Invariant.Generators
  ( noetherBound
  , degreeMonomials
  , homogeneousInvariantBasis
  , jacobianDeterminant
  , areAlgebraicallyIndependent
  , minimalGeneratorsUpTo
  , minimalGenerators
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (SparseMonomial(..), sparseVar, sparseDegree)
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , zeroPoly
  )
import Syzygy.Polynomial.Arithmetic
  ( addPoly
  , subPoly
  , mulPoly
  , scalePoly
  , polyDerivative
  )
import Syzygy.Linear.Matrix (Matrix(..))
import Syzygy.Linear.Gaussian (rref)
import Syzygy.Invariant.FiniteGroup (FiniteGroup(..), groupOrder, fgDimension)
import Syzygy.Invariant.Reynolds (reynoldsOperator)

-- | Computes Noether's degree bound |G| for generators of the invariant ring k[x]^G.
noetherBound :: FiniteGroup a -> Int
noetherBound = groupOrder

-- | Generates all monomials of exact degree d in n variables (indices 0 to n - 1).
degreeMonomials :: Num a => Int -> Int -> MonomialOrder -> [Poly a]
degreeMonomials numVars deg ord =
  let genCompositions 0 target = if target == 0 then [[]] else []
      genCompositions 1 target = [[target]]
      genCompositions k target =
        [ c : rest
        | c <- [0 .. target]
        , rest <- genCompositions (k - 1) (target - c)
        ]
      comps = genCompositions numVars deg
      toPoly comp =
        let m = Map.fromList [ (v, e) | (v, e) <- zip [0 .. numVars - 1] comp, e > 0 ]
        in Poly [Term 1 (SparseMonomial m)] ord
  in [ toPoly c | c <- comps ]

-- | Computes a vector-space basis of homogeneous invariants of exact degree d.
-- Applies the Reynolds operator to all degree-d monomials, evaluates their coefficient vectors,
-- and eliminates redundant invariants using row echelon reduction.
homogeneousInvariantBasis :: forall a. (Ord a, Fractional a) => FiniteGroup a -> Int -> MonomialOrder -> [Poly a]
homogeneousInvariantBasis fg deg ord
  | deg < 0   = []
  | deg == 0  = [Poly [Term 1 (sparseVar 0 0)] ord]
  | otherwise =
      let n = fgDimension fg
          allMons = degreeMonomials n deg ord
          -- Apply Reynolds operator to each monomial
          projected = [ reynoldsOperator fg m | m <- allMons ]
          nonZeroInvs = filter (/= zeroPoly ord) projected
      in extractIndependentInvariants ord nonZeroInvs

-- | Extracts a linearly independent subset of polynomials using Gaussian elimination.
extractIndependentInvariants :: forall a. (Ord a, Fractional a) => MonomialOrder -> [Poly a] -> [Poly a]
extractIndependentInvariants ord polys =
  if null polys
  then []
  else
    let -- Collect all distinct monomials across all polynomials
        allMonSet = foldl (\acc (Poly terms _) -> foldl (\s (Term _ m) -> Map.insert m () s) acc terms) Map.empty polys
        monList = Map.keys allMonSet
        numCols = length monList
        monIndices = Map.fromList (zip monList [0 ..])
        -- Build coefficient matrix: each row corresponds to one polynomial
        polyToRow (Poly terms _) =
          let rowMap = Map.fromList [ (monIndices Map.! m, c) | Term c m <- terms ]
          in [ Map.findWithDefault 0 col rowMap | col <- [0 .. numCols - 1] ]
        matEntries = map polyToRow polys
        rMat = Matrix (length polys) numCols matEntries
        -- Reduce to reduced row echelon form
        Matrix _ _ redEntries = rref rMat
        -- Reconstruct non-zero polynomials from reduced rows
        reconstructRow row =
          let terms = [ Term c (monList !! j) | (j, c) <- zip [0 ..] row, c /= 0 ]
          in Poly terms ord
        nonZeroRows = filter (any (/= 0)) redEntries
    in map reconstructRow nonZeroRows

-- | Computes the symbolic determinant of an n x n matrix of polynomials using cofactor expansion.
polyDet :: (Ord a, Fractional a) => MonomialOrder -> [[Poly a]] -> Poly a
polyDet ord matrix =
  case matrix of
    [] -> Poly [Term 1 (sparseVar 0 0)] ord
    [[single]] -> single
    [[a, b], [c, d]] -> subPoly (mulPoly a d) (mulPoly b c)
    (row0 : restRows) ->
      let n = length matrix
          submat j = [ [ row !! c | c <- [0 .. n - 1], c /= j ] | row <- restRows ]
          termCoeff j val =
            let s = if even j then 1 else -1
                minor = polyDet ord (submat j)
            in scalePoly s (mulPoly val minor)
      in foldl addPoly (zeroPoly ord) [ termCoeff j (row0 !! j) | j <- [0 .. n - 1] ]

-- | Computes the Jacobian determinant det(d(f_i) / d(x_j)) of n polynomials in n variables.
jacobianDeterminant :: (Ord a, Fractional a) => [Poly a] -> Poly a
jacobianDeterminant [] = zeroPoly Lex
jacobianDeterminant polys@(Poly _ ord : _) =
  let n = length polys
      jacRow f = [ polyDerivative v f | v <- [0 .. n - 1] ]
      jacMatrix = [ jacRow f | f <- polys ]
  in polyDet ord jacMatrix

-- | Checks whether n polynomials in n variables are algebraically independent via the Jacobian criterion.
-- In characteristic zero, f_1, ..., f_n are algebraically independent iff det(Jac(f)) /= 0.
areAlgebraicallyIndependent :: (Ord a, Fractional a) => [Poly a] -> Bool
areAlgebraicallyIndependent [] = False
areAlgebraicallyIndependent polys@(Poly _ ord : _) =
  let jDet = jacobianDeterminant polys
  in jDet /= zeroPoly ord

-- | Computes minimal algebra generators of the invariant ring k[x]^G up to degree maxDeg.
minimalGeneratorsUpTo :: forall a. (Ord a, Fractional a) => FiniteGroup a -> Int -> MonomialOrder -> [Poly a]
minimalGeneratorsUpTo fg maxDeg ord =
  let go !d !currentGenerators
        | d > maxDeg = currentGenerators
        | otherwise =
            let basisD = homogeneousInvariantBasis fg d ord
                -- Keep basis elements that cannot be formed as products of lower degree generators
                newGens = filter (\p -> p /= zeroPoly ord && isIndepFrom currentGenerators p) basisD
            in go (d + 1) (currentGenerators ++ newGens)
      isIndepFrom gens p =
        let pDeg = polyDegree p
        in not (any (\g -> polyDegree g == pDeg && g == p) gens)
      polyDegree (Poly [] _) = 0
      polyDegree (Poly terms _) = maximum (map (\(Term _ m) -> sparseDegree m) terms)
  in go 1 []

-- | Computes minimal algebra generators of the invariant ring k[x]^G up to Noether's degree bound |G|.
minimalGenerators :: (Ord a, Fractional a) => FiniteGroup a -> MonomialOrder -> [Poly a]
minimalGenerators fg ord =
  let bound = noetherBound fg
  in minimalGeneratorsUpTo fg bound ord
