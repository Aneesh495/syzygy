{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Singularities
  ( jacobianMatrix
  , singularLocusIdeal
  , isSmoothVariety
  , singularLocusDimension
  , matrixMinors
  , polynomialDeterminant
  ) where

import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Polynomial.Base
  ( Poly(..)
  , zeroPoly
  , isZeroPoly
  )
import Syzygy.Polynomial.Arithmetic
  ( addPoly
  , subPoly
  , mulPoly
  , polyDerivative
  )
import Syzygy.Geometry.Ideal
  ( Ideal(..)
  , idealFromGenerators
  , idealSum
  , isProperIdeal
  )
import Syzygy.Geometry.Dimension (idealDimension, idealCodimension)

-- | Computes the m x n Jacobian matrix of partial derivatives df_i / dx_j
-- for a list of m polynomials [f_0, ..., f_{m-1}] in n variables [0 .. n-1].
jacobianMatrix :: (Eq a, Num a) => Int -> [Poly a] -> [[Poly a]]
jacobianMatrix numVars polys =
  [ [ polyDerivative varId f | varId <- [0 .. numVars - 1] ]
  | f <- polys
  ]

-- | Computes the determinant of a square k x k matrix of polynomials using Leibniz expansion.
polynomialDeterminant :: forall a. (Eq a, Num a) => MonomialOrder -> [[Poly a]] -> Poly a
polynomialDeterminant polyOrder mat =
  case mat of
    [] -> zeroPoly polyOrder
    [[p]] -> p
    [[a, b], [c, d]] ->
      subPoly (mulPoly a d) (mulPoly b c)
    (firstRow:restRows) ->
      let n = length mat
          cols = [0 .. n - 1]
      in foldl addPoly (zeroPoly polyOrder)
           [ let subMat = [ [ row !! c | c <- cols, c /= j ] | row <- restRows ]
                 subDet = polynomialDeterminant polyOrder subMat
                 term = mulPoly (firstRow !! j) subDet
             in if even j then term else subPoly (zeroPoly polyOrder) term
           | j <- cols
           ]

-- | Generates all k x k minors (determinants of k x k submatrices) of an m x n polynomial matrix.
matrixMinors :: (Eq a, Num a) => Int -> MonomialOrder -> [[Poly a]] -> [Poly a]
matrixMinors k polyOrder mat
  | k <= 0 = []
  | otherwise =
      let m = length mat
          n = case mat of (r:_) -> length r; [] -> 0
      in if k > m || k > n
         then []
         else
           let rowCombos = combinations k [0 .. m - 1]
               colCombos = combinations k [0 .. n - 1]
           in [ let subMatrix = [ [ (mat !! r) !! c | c <- cols ] | r <- rows ]
                in polynomialDeterminant polyOrder subMatrix
              | rows <- rowCombos
              , cols <- colCombos
              ]
  where
    combinations 0 _ = [[]]
    combinations _ [] = []
    combinations r (x:xs) =
      map (x :) (combinations (r - 1) xs) ++ combinations r xs

-- | Computes the singular locus ideal of an affine variety V(I):
-- Sing(I) = I + <all c x c minors of the Jacobian matrix J(I)>,
-- where c = codim(I) = n - dim(V(I)).
-- By the Jacobian Criterion, p in V(I) is singular iff all c x c minors vanish at p.
singularLocusIdeal :: forall a. (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
singularLocusIdeal numVars ideal@(Ideal gens order) =
  let nonZeroGens = filter (not . isZeroPoly) gens
      c = idealCodimension numVars nonZeroGens
  in if c <= 0
     then ideal
     else
       let jac = jacobianMatrix numVars nonZeroGens
           minors = matrixMinors c order jac
           nonZeroMinors = filter (not . isZeroPoly) minors
           minorId = idealFromGenerators order nonZeroMinors
       in idealSum ideal minorId

-- | Checks whether an affine variety V(I) is smooth (non-singular):
-- V(I) is smooth iff the singular locus is empty, which holds iff Sing(I) == <1>.
isSmoothVariety :: (Ord a, Fractional a) => Int -> Ideal a -> Bool
isSmoothVariety numVars ideal =
  let singId = singularLocusIdeal numVars ideal
  in not (isProperIdeal singId)

-- | Computes the Krull dimension of the singular locus of V(I).
-- Returns -1 if the variety is completely smooth (singular locus is empty).
singularLocusDimension :: (Ord a, Fractional a) => Int -> Ideal a -> Int
singularLocusDimension numVars ideal =
  let singId = singularLocusIdeal numVars ideal
  in idealDimension numVars (idealGens singId)
