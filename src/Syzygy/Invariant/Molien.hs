{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Invariant.Molien
Description : Molien's formula and Hilbert-Poincare series for invariant rings
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Molien's formula (1897) expresses the Hilbert-Poincare series of the graded
ring of invariants k[x_0, ..., x_{n-1}]^G of a finite matrix group G in GL(n, k)
as an average of rational functions over the group:

  Phi(t) = sum_{d=0}^infinity dim(k[x]^G_d) * t^d
         = (1 / |G|) * sum_{g in G} (1 / det(I - t * g))

By the Hochster-Eagon Theorem (1971), k[x]^G is a Cohen-Macaulay graded ring,
meaning the Molien series can always be written as a rational function
P(t) / product_{i=1}^n (1 - t^{d_i}) where d_i are the degrees of primary invariants.
-}
module Syzygy.Invariant.Molien
  ( molienTerm
  , molienSeries
  , molienDimension
  , molienDimensions
  ) where

import Syzygy.Linear.Matrix (Matrix(..))
import Syzygy.Solver.Multiplication (characteristicPolynomial)
import Syzygy.Invariant.FiniteGroup (FiniteGroup(..), groupOrder)

-- | Computes the Taylor series expansion of 1 / det(I - t * g) up to degree maxDeg.
-- Returns the coefficients [a_0, a_1, ..., a_{maxDeg}].
--
-- Since det(I - t * g) = 1 + p_1 * t + ... + p_n * t^n, the reciprocal
-- satisfies a_0 = 1 and a_m = - sum_{k=1}^{min(n, m)} p_k * a_{m-k}.
molienTerm :: forall a. (Eq a, Fractional a) => Matrix a -> Int -> [a]
molienTerm g maxDeg =
  let -- Char poly gives [c_0, c_1, ..., c_n] for det(lambda * I - g) where c_n = 1
      charCoeffs = characteristicPolynomial g
      -- Reverse gives coefficients [p_0, p_1, ..., p_n] of det(I - t * g)
      pCoeffs = reverse charCoeffs
      n = length pCoeffs - 1
      -- Expand reciprocal power series up to maxDeg
      go !m !acc
        | m > maxDeg = reverse acc
        | m == 0     = go 1 [1]
        | otherwise  =
            let -- acc is [a_{m-1}, a_{m-2}, ..., a_0]
                recipStep k
                  | k > n || k > m = 0
                  | otherwise =
                      let pk = pCoeffs !! k
                          amk = acc !! (k - 1)
                      in pk * amk
                val = - sum [ recipStep k | k <- [1 .. min n m] ]
            in go (m + 1) (val : acc)
  in if maxDeg < 0 then [] else go 0 []

-- | Computes the Molien series of the invariant ring k[x]^G up to degree maxDeg:
-- Phi(t) = (1 / |G|) * sum_{g in G} (1 / det(I - t * g))
-- Returns the list of dimensions [dim(k[x]^G_0), dim(k[x]^G_1), ..., dim(k[x]^G_{maxDeg})].
molienSeries :: (Eq a, Fractional a) => FiniteGroup a -> Int -> [a]
molienSeries fg maxDeg =
  let gList = fgElements fg
      orderG = groupOrder fg
  in if orderG <= 0 || maxDeg < 0
     then []
     else
       let terms = [ molienTerm g maxDeg | g <- gList ]
           -- Sum coefficient-wise across all group elements
           sumTerms = foldl1 (zipWith (+)) terms
           scaleFactor = recip (fromIntegral orderG)
       in map (* scaleFactor) sumTerms

-- | Computes the dimension of the space of invariants of degree d for group G.
molienDimension :: (Eq a, Fractional a) => FiniteGroup a -> Int -> a
molienDimension fg d =
  let series = molienSeries fg d
  in if d < length series then series !! d else 0

-- | Computes the list of dimensions of invariant subspaces from degree 0 to maxDeg.
molienDimensions :: (Eq a, Fractional a) => FiniteGroup a -> Int -> [a]
molienDimensions = molienSeries
