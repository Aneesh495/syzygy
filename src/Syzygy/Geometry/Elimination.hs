{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Elimination
  ( eliminateVariables
  , eliminationIdeal
  , eliminateK
  , freeOfVariables
  , substituteVariable
  , projectVarietyClosure
  ) where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (SparseMonomial(..))
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , reorderPoly
  )
import Syzygy.Polynomial.Arithmetic (addPoly)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Geometry.Ideal (Ideal(..), idealFromGenerators)

-- | Checks whether a polynomial is free of the specified variables.
freeOfVariables :: [Int] -> Poly a -> Bool
freeOfVariables vars (Poly terms _) =
  let varSet = Set.fromList vars
  in not (any (\(Term _ (SparseMonomial m)) -> any (`Set.member` varSet) (Map.keys m)) terms)

-- | Eliminates the first k variables [0 .. k-1] from an ideal, computing:
-- I_k = I cap k[x_k, ..., x_{n-1}].
-- Uses the Elimination Theorem: a Lex Groebner basis or BlockOrder Groebner basis
-- restricted to polynomials not containing x_0, ..., x_{k-1} forms a Groebner basis of I_k.
eliminationIdeal :: (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
eliminationIdeal k (Ideal gens order)
  | k <= 0 = Ideal gens order
  | otherwise =
      let elimOrder = BlockOrder k Lex order
          reorderedGens = [ reorderPoly elimOrder g | g <- gens ]
          gb = reducedGrobnerBasis (buchberger reorderedGens)
          elimVars = [0 .. k - 1]
          filtered = [ reorderPoly order p | p <- gb, freeOfVariables elimVars p ]
      in idealFromGenerators order filtered

-- | Eliminates an arbitrary subset of variables from an ideal.
eliminateVariables :: (Ord a, Fractional a) => [Int] -> Ideal a -> Ideal a
eliminateVariables vars ideal@(Ideal gens order)
  | null vars = ideal
  | otherwise =
      let sortedVars = Set.toList (Set.fromList vars)
          elimOrder = Lex
          reorderedGens = [ reorderPoly elimOrder g | g <- gens ]
          gb = reducedGrobnerBasis (buchberger reorderedGens)
          filtered = [ reorderPoly order p | p <- gb, freeOfVariables sortedVars p ]
      in idealFromGenerators order filtered

-- | Low-level helper: eliminates the first k variables directly from a list of polynomials.
eliminateK :: (Ord a, Fractional a) => Int -> MonomialOrder -> [Poly a] -> [Poly a]
eliminateK k targetOrder polys =
  let elimOrder = BlockOrder k Lex targetOrder
      reordered = [ reorderPoly elimOrder p | p <- polys ]
      gb = reducedGrobnerBasis (buchberger reordered)
      elimVars = [0 .. k - 1]
  in [ reorderPoly targetOrder p | p <- gb, freeOfVariables elimVars p ]

-- | Evaluates / specializes a variable x_{varId} at a scalar value c in k: p(..., x_{varId} = c, ...).
substituteVariable :: forall a. (Eq a, Num a) => Int -> a -> Poly a -> Poly a
substituteVariable varId val (Poly terms order) =
  let evaluatedTerms = [ evalTerm coeff m | Term coeff m <- terms ]
  in foldl addPoly (Poly [] order) evaluatedTerms
  where
    evalTerm c (SparseMonomial m) =
      case Map.lookup varId m of
        Nothing -> Poly [Term c (SparseMonomial m)] order
        Just e  ->
          let c' = c * (val ^ e)
              m' = SparseMonomial (Map.delete varId m)
          in if c' == 0 then Poly [] order else Poly [Term c' m'] order

-- | Computes the Zariski closure of the projection of affine variety V(I)
-- onto the coordinate subspace (x_k, ..., x_{n-1}):
-- \overline{\pi_k(V(I))} = V(I_k).
projectVarietyClosure :: (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
projectVarietyClosure = eliminationIdeal
