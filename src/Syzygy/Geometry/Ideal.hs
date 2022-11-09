{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Ideal
  ( Ideal(..)
  , idealFromGenerators
  , idealSum
  , idealProduct
  , idealPower
  , idealIntersection
  , idealContains
  , idealEqual
  , isProperIdeal
  , isZeroIdeal
  , idealGroebnerBasis
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  , sparseMaxVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  )
import Syzygy.Polynomial.Arithmetic (mulPoly)
import Syzygy.Polynomial.Division (normalForm)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)

-- | Represents an ideal I = <f_1, ..., f_r> in the polynomial ring R = k[x_0, ..., x_{n-1}].
data Ideal a = Ideal
  { idealGens  :: ![Poly a]
  , idealOrder :: !MonomialOrder
  } deriving (Eq, Show, Generic)

-- | Constructs an ideal from a list of generators.
idealFromGenerators :: (Ord a, Fractional a) => MonomialOrder -> [Poly a] -> Ideal a
idealFromGenerators order gens =
  let nonZeroGens = filter (not . isZeroPoly) gens
  in Ideal nonZeroGens order

-- | Computes the reduced Groebner basis of the ideal.
idealGroebnerBasis :: (Ord a, Fractional a) => Ideal a -> [Poly a]
idealGroebnerBasis (Ideal gens _) =
  reducedGrobnerBasis (buchberger gens)

-- | Sum of two ideals: I + J = <f_1, ..., f_r, g_1, ..., g_s>.
idealSum :: (Ord a, Fractional a) => Ideal a -> Ideal a -> Ideal a
idealSum (Ideal gens1 order) (Ideal gens2 _) =
  idealFromGenerators order (gens1 ++ gens2)

-- | Product of two ideals: I * J = <f_i * g_j | f_i in I, g_j in J>.
idealProduct :: (Ord a, Fractional a) => Ideal a -> Ideal a -> Ideal a
idealProduct (Ideal gens1 order) (Ideal gens2 _) =
  let products = [ mulPoly f g | f <- gens1, g <- gens2 ]
  in idealFromGenerators order products

-- | Power of an ideal: I^k.
idealPower :: (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
idealPower k ideal@(Ideal _ order)
  | k <= 0 = Ideal [Poly [Term 1 sparseOne] order] order
  | k == 1 = ideal
  | otherwise = idealProduct ideal (idealPower (k - 1) ideal)

-- | Intersection of two ideals I cap J using the auxiliary variable t:
-- I cap J = (t * I + (1 - t) * J) cap k[x_0, ..., x_{n-1}].
idealIntersection :: forall a. (Ord a, Fractional a) => Ideal a -> Ideal a -> Ideal a
idealIntersection (Ideal gens1 order) (Ideal gens2 _) =
  let maxVar1 = maximum (-1 : [ sparseMaxVar m | f <- gens1, Term _ m <- polyTerms f ])
      maxVar2 = maximum (-1 : [ sparseMaxVar m | g <- gens2, Term _ m <- polyTerms g ])
      tVar = max maxVar1 maxVar2 + 1
      -- Construct t and (1 - t)
      tPoly = Poly [Term 1 (sparseVar tVar 1)] elimOrder
      oneMinusTPoly = Poly [Term 1 sparseOne, Term (-1) (sparseVar tVar 1)] elimOrder
      -- Embed generators into block elimination order with t in block 1
      elimOrder = BlockOrder 1 Lex order
      tGens1 = [ mulPoly tPoly (embedPoly f) | f <- gens1 ]
      oneMinusTGens2 = [ mulPoly oneMinusTPoly (embedPoly g) | g <- gens2 ]
      bigGB = buchberger (tGens1 ++ oneMinusTGens2)
      -- Eliminate t: retain polynomials whose terms do not contain tVar
      eliminated = [ projectPoly p
                   | p <- bigGB
                   , not (containsVar tVar p)
                   ]
  in idealFromGenerators order (filter (not . isZeroPoly) eliminated)
  where
    polyTerms (Poly terms _) = terms

    embedPoly (Poly terms _) =
      Poly terms (BlockOrder 1 Lex order)

    projectPoly (Poly terms _) =
      Poly terms order

    containsVar v (Poly terms _) =
      any (\(Term _ (SparseMonomial m)) -> Map.member v m) terms

-- | Tests ideal membership: f in I iff normalForm(f, GB(I)) == 0.
idealContains :: (Ord a, Fractional a) => Ideal a -> Poly a -> Bool
idealContains ideal f =
  let gb = idealGroebnerBasis ideal
  in isZeroPoly (normalForm f gb)

-- | Tests ideal equality: I == J iff reducedGB(I) == reducedGB(J).
idealEqual :: (Ord a, Fractional a) => Ideal a -> Ideal a -> Bool
idealEqual ideal1 ideal2 =
  let gb1 = idealGroebnerBasis ideal1
      gb2 = idealGroebnerBasis ideal2
  in gb1 == gb2

-- | Checks whether an ideal is proper: I /= <1>.
isProperIdeal :: (Ord a, Fractional a) => Ideal a -> Bool
isProperIdeal ideal@(Ideal _ order) =
  not (idealContains ideal (Poly [Term 1 sparseOne] order))

-- | Checks whether an ideal is the zero ideal <0>.
isZeroIdeal :: (Ord a, Fractional a) => Ideal a -> Bool
isZeroIdeal (Ideal gens _) =
  all isZeroPoly gens
