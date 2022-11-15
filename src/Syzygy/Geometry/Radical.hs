{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Radical
  ( isInRadical
  , radicalZeroDimensional
  , isRadicalIdeal
  , squareFreePart
  ) where

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
  , reorderPoly
  )
import Syzygy.Polynomial.Arithmetic
  ( subPoly
  , mulPoly
  , polyDerivative
  )
import Syzygy.Polynomial.MultivariateGCD (multivariateGCD)
import Syzygy.Polynomial.Division (multivariateDiv)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Geometry.Ideal
  ( Ideal(..)
  , idealFromGenerators
  , idealEqual
  , idealSum
  )

-- | Tests radical ideal membership f in sqrt(I) using Rabinowitsch's Trick (1929).
-- Introduces a slack variable w and computes the Groebner basis of <I, 1 - w * f>.
-- By Hilbert's Nullstellensatz, f in sqrt(I) iff 1 in <I, 1 - w * f>.
isInRadical :: forall a. (Ord a, Fractional a) => Ideal a -> Poly a -> Bool
isInRadical (Ideal gens order) f
  | isZeroPoly f = True
  | null gens = False
  | otherwise =
      let maxVarGens = maximum (-1 : [ sparseMaxVar m | g <- gens, Term _ m <- polyTerms g ])
          maxVarF = maximum (-1 : [ sparseMaxVar m | Term _ m <- polyTerms f ])
          wVar = max maxVarGens maxVarF + 1
          elimOrder = BlockOrder 1 Lex order

          wPoly = Poly [Term 1 (sparseVar wVar 1)] elimOrder
          onePoly = Poly [Term 1 sparseOne] elimOrder
          fEmbed = reorderPoly elimOrder f
          rabinowitschGen = subPoly onePoly (mulPoly wPoly fEmbed)
          gensEmbed = [ reorderPoly elimOrder g | g <- gens ]

          gb = reducedGrobnerBasis (buchberger (rabinowitschGen : gensEmbed))
      in any isOne gb
  where
    polyTerms (Poly terms _) = terms
    isOne (Poly [Term c (SparseMonomial m)] _) = c /= 0 && Map.null m
    isOne _ = False

-- | Computes the square-free part of a polynomial: f_red = f / gcd(f, df/dx).
squareFreePart :: forall a. (Ord a, Fractional a) => Int -> Poly a -> Poly a
squareFreePart varId f =
  let df = polyDerivative varId f
      g = multivariateGCD f df
  in case multivariateDiv f [g] of
       ([q], remF) | isZeroPoly remF -> q
       _ -> f

-- | Computes the radical sqrt(I) of a zero-dimensional ideal using Seidenberg's Lemma 92:
-- For each coordinate x_i, let g_i be the monic generator of I cap k[x_i].
-- Let h_i = g_i / gcd(g_i, g_i') be its square-free reduction.
-- Then sqrt(I) = I + <h_0(x_0), ..., h_{n-1}(x_{n-1})>.
radicalZeroDimensional :: forall a. (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
radicalZeroDimensional numVars ideal@(Ideal _ order) =
  let -- For each variable i, eliminate all other variables to find univariate generator in k[x_i]
      univariateGens = [ findUnivariateGen i | i <- [0 .. numVars - 1] ]
      sqFreeGens = [ squareFreePart i g | (i, g) <- zip [0 .. numVars - 1] univariateGens ]
      sqFreeIdeal = idealFromGenerators order sqFreeGens
  in idealSum ideal sqFreeIdeal
  where
    findUnivariateGen varId =
      let elimOrder = Lex
          reorderedGens = [ reorderPoly elimOrder g | g <- idealGens ideal ]
          gb = reducedGrobnerBasis (buchberger reorderedGens)
          -- Find element purely in x_{varId}
          univariateCands = filter (isPureInVar varId) gb
      in case univariateCands of
           (cand:_) -> reorderPoly order cand
           []       -> Poly [Term 1 (sparseVar varId 1)] order -- fallback

    isPureInVar vId (Poly terms _) =
      all (\(Term _ (SparseMonomial m)) ->
            case Map.keys m of
              []  -> True
              [v] -> v == vId
              _   -> False
          ) terms

-- | Checks whether a zero-dimensional ideal is radical: I == sqrt(I).
isRadicalIdeal :: (Ord a, Fractional a) => Int -> Ideal a -> Bool
isRadicalIdeal numVars ideal =
  let rad = radicalZeroDimensional numVars ideal
  in idealEqual ideal rad
