{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Dimension
  ( idealDimension
  , idealCodimension
  , independentVariableSubsets
  , maximallyIndependentSubsets
  , isIndependentSet
  ) where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Syzygy.Monomial.Sparse (SparseMonomial(..))
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  , leadingMonomial
  )

-- | Checks whether a subset of variables U subset {0, ..., n-1} is independent
-- modulo the leading ideal LT(G).
-- U is independent iff no leading monomial of G is contained entirely in the polynomial
-- subring k[U].
isIndependentSet :: [SparseMonomial] -> [Int] -> Bool
isIndependentSet lms uVars =
  let uSet = Set.fromList uVars
  in not (any (inSubring uSet) lms)
  where
    inSubring uSet (SparseMonomial m)
      | Map.null m = False
      | otherwise  = all (`Set.member` uSet) (Map.keys m)

-- | Finds all independent subsets of variables modulo the Groebner basis G.
independentVariableSubsets :: (Eq a, Fractional a) => Int -> [Poly a] -> [[Int]]
independentVariableSubsets numVars gb =
  let nonZeroGB = filter (not . isZeroPoly) gb
      lms = [ m | g <- nonZeroGB, Just m <- [leadingMonomial g] ]
      allSubsets = subsets [0 .. numVars - 1]
  in filter (isIndependentSet lms) allSubsets
  where
    subsets [] = [[]]
    subsets (x:xs) =
      let subs = subsets xs
      in subs ++ map (x :) subs

-- | Computes the Krull dimension of the quotient ring R/I (dimension of the affine variety V(I)).
-- By Groebner basis dimension theory:
-- dim(R/I) = max { |U| | U is an independent set of variables modulo LT(G) }.
-- For the zero ideal <0>, returns n. For the unit ideal <1>, returns -1.
idealDimension :: (Eq a, Fractional a) => Int -> [Poly a] -> Int
idealDimension numVars gb =
  let nonZeroGB = filter (not . isZeroPoly) gb
  in if any isOne nonZeroGB
     then -1
     else
       let indepSets = independentVariableSubsets numVars nonZeroGB
       in if null indepSets then 0 else maximum (map length indepSets)
  where
    isOne (Poly [Term c (SparseMonomial m)] _) = c /= 0 && Map.null m
    isOne _ = False

-- | Computes the codimension of an ideal: codim(I) = n - dim(R/I).
idealCodimension :: (Eq a, Fractional a) => Int -> [Poly a] -> Int
idealCodimension numVars gb =
  let d = idealDimension numVars gb
  in if d < 0 then numVars + 1 else numVars - d

-- | Returns the maximally independent sets of variables (those achieving the maximal dimension).
maximallyIndependentSubsets :: (Eq a, Fractional a) => Int -> [Poly a] -> [[Int]]
maximallyIndependentSubsets numVars gb =
  let d = idealDimension numVars gb
      allIndep = independentVariableSubsets numVars gb
  in filter (\u -> length u == d) allIndep
