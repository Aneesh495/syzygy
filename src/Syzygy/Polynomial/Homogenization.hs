{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Homogenization
  ( homogenizePoly
  , dehomogenizePoly
  , homogeneousComponents
  , isHomogeneous
  , projectivizePoly
  , affineSlicePoly
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse
  ( sparseDegree
  , sparseSetExp
  , sparseSupport
  , sparseGetExp
  , sparseOne
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , isZeroPoly
  , polyTerms
  , polyOrder
  , polyDegree
  , fromTerms
  )
import Syzygy.Monomial.Operations (monomialShift)

-- | Homogenizes a multivariate polynomial f with respect to a new homogenization variable x_0:
-- f^h(x_0, x_1, ..., x_n) = x_0^d * f(x_1/x_0, ..., x_n/x_0) where d = deg(f).
-- Every variable index in f is shifted by +1 (x_i -> x_{i+1}),
-- and x_0 is inserted with exponent d - deg(m) so that all terms have identical total degree d.
homogenizePoly :: forall c. (Eq c, Num c) => Poly c -> Poly c
homogenizePoly p
  | isZeroPoly p = p
  | otherwise =
      let ord = polyOrder p
          d = polyDegree p
          hTerms = [ Term c (makeHomogeneousMonomial d m) | Term c m <- polyTerms p ]
      in fromTerms ord hTerms
  where
    makeHomogeneousMonomial targetDeg m =
      let shiftedM = monomialShift 1 m
          degM = sparseDegree m
          exp0 = targetDeg - degM
      in if exp0 == 0
         then shiftedM
         else sparseSetExp 0 exp0 shiftedM

-- | Dehomogenizes a homogeneous polynomial by evaluating x_0 = 1
-- and shifting variable indices back: x_{i+1} -> x_i.
dehomogenizePoly :: forall c. (Eq c, Num c) => Poly c -> Poly c
dehomogenizePoly (Poly terms ord) =
  let dehomogTerms = [ Term c (unshiftMonomial m) | Term c m <- terms ]
  in fromTerms ord dehomogTerms
  where
    unshiftMonomial m =
      let supp = sparseSupport m
          suppWithout0 = filter (> 0) supp
      in foldl (\acc varId ->
                  let expVal = sparseGetExp varId m
                  in sparseSetExp (varId - 1) expVal acc
               ) sparseOne suppWithout0

-- | Decomposes a polynomial f into its graded homogeneous components:
-- f = f_0 + f_1 + ... + f_d where each f_k contains terms of total degree exactly k.
-- Returns a map from degree k to the homogeneous polynomial component f_k.
homogeneousComponents :: forall c. (Eq c, Num c) => Poly c -> Map.Map Int (Poly c)
homogeneousComponents p@(Poly terms ord)
  | isZeroPoly p = Map.empty
  | otherwise =
      let grouped = foldl insertTerm Map.empty terms
      in Map.map (fromTerms ord) grouped
  where
    insertTerm accMap t@(Term _ m) =
      let d = sparseDegree m
      in Map.insertWith (++) d [t] accMap

-- | Tests if a polynomial is homogeneous (all non-zero terms have identical total degree).
isHomogeneous :: Poly c -> Bool
isHomogeneous p
  | isZeroPoly p = True
  | otherwise =
      let terms = polyTerms p
          degs = [ sparseDegree m | Term _ m <- terms ]
      in case degs of
           [] -> True
           (d0:ds) -> all (== d0) ds

-- | Projectivizes an affine polynomial into projective space P^n.
projectivizePoly :: (Eq c, Num c) => Poly c -> Poly c
projectivizePoly = homogenizePoly

-- | Affine slice of a projective variety: sets variable x_k = 1.
affineSlicePoly :: forall c. (Eq c, Num c) => Int -> Poly c -> Poly c
affineSlicePoly varId (Poly terms ord) =
  let slicedTerms = [ Term c (sparseSetExp varId 0 m) | Term c m <- terms ]
  in fromTerms ord slicedTerms
