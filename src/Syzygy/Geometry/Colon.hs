{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Colon
  ( colonIdeal
  , colonPoly
  , saturation
  , saturationPoly
  , saturateIrrelevant
  ) where

import Syzygy.Monomial.Sparse (sparseVar)
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  )
import Syzygy.Polynomial.Division (multivariateDiv)
import Syzygy.Geometry.Ideal
  ( Ideal(..)
  , idealFromGenerators
  , idealIntersection
  , idealEqual
  )

-- | Computes the colon ideal I : <g> for a single polynomial g:
-- I : <g> = (1 / g) * (I cap <g>).
-- If g == 0, returns <1>.
colonPoly :: forall a. (Ord a, Fractional a) => Ideal a -> Poly a -> Ideal a
colonPoly (Ideal gens order) g
  | isZeroPoly g = Ideal [Poly [Term 1 mempty] order] order
  | otherwise =
      let gIdeal = Ideal [g] order
          interIdeal = idealIntersection (Ideal gens order) gIdeal
          interGB = idealGens interIdeal
          -- Divide each generator in I cap <g> by g:
          quotients = [ case multivariateDiv h [g] of
                          ([q], remH) | isZeroPoly remH -> q
                          _ -> h
                      | h <- interGB
                      ]
      in idealFromGenerators order (filter (not . isZeroPoly) quotients)

-- | Computes the ideal quotient (colon ideal) I : J for two ideals I and J:
-- I : J = bigcap_{g in gens(J)} (I : <g>).
colonIdeal :: (Ord a, Fractional a) => Ideal a -> Ideal a -> Ideal a
colonIdeal idealI (Ideal gensJ _)
  | null gensJ = idealI
  | otherwise =
      let colons = [ colonPoly idealI g | g <- gensJ, not (isZeroPoly g) ]
      in case colons of
           []     -> idealI
           (c:cs) -> foldl idealIntersection c cs

-- | Computes the saturation of an ideal I with respect to a polynomial g:
-- I : g^infty = cup_{k >= 1} (I : g^k).
-- Iteratively computes I_0 = I, I_{k+1} = I_k : <g> until I_{k+1} == I_k.
saturationPoly :: (Ord a, Fractional a) => Ideal a -> Poly a -> Ideal a
saturationPoly initialIdeal g =
  go initialIdeal
  where
    go currentIdeal =
      let nextIdeal = colonPoly currentIdeal g
      in if idealEqual currentIdeal nextIdeal
         then currentIdeal
         else go nextIdeal

-- | Computes the ideal saturation I : J^infty:
-- I : J^infty = cup_{k >= 1} (I : J^k).
-- Iteratively computes I_{k+1} = I_k : J until stabilization (by Hilbert basis theorem).
saturation :: (Ord a, Fractional a) => Ideal a -> Ideal a -> Ideal a
saturation initialIdeal idealJ =
  go initialIdeal
  where
    go currentIdeal =
      let nextIdeal = colonIdeal currentIdeal idealJ
      in if idealEqual currentIdeal nextIdeal
         then currentIdeal
         else go nextIdeal

-- | Computes the projective saturation of a homogeneous ideal with respect to
-- the irrelevant ideal m = <x_0, ..., x_{n-1}>:
-- I : m^infty.
-- Strips out the embedded component supported at the vertex of the affine cone.
saturateIrrelevant :: (Ord a, Fractional a) => Int -> Ideal a -> Ideal a
saturateIrrelevant numVars ideal@(Ideal _ order) =
  let irrelevantGens = [ Poly [Term 1 (sparseVar v 1)] order | v <- [0 .. numVars - 1] ]
      irrelevantIdeal = Ideal irrelevantGens order
  in saturation ideal irrelevantIdeal
