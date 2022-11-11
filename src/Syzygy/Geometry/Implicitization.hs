{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.Implicitization
  ( implicitizePolynomial
  , implicitizeRational
  , twistedCubicIdeal
  , whitneyUmbrellaIdeal
  , enneperSurfaceIdeal
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  , reorderPoly
  )
import Syzygy.Polynomial.Arithmetic (subPoly, mulPoly)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Geometry.Ideal (Ideal(..), idealFromGenerators)
import Syzygy.Geometry.Elimination (freeOfVariables)

-- | Implicitization of polynomial parametrizations:
-- Given coordinate parametrizations [f_0(t), f_1(t), ..., f_{m-1}(t)] in parameters t_0, ..., t_{k-1},
-- computes the implicit polynomial relations F(x_0, ..., x_{m-1}) = 0 defining the variety.
-- Shifts coordinate variables x_i to indices [numParams .. numParams + m - 1]
-- and eliminates the parameter variables [0 .. numParams - 1].
implicitizePolynomial :: forall a. (Ord a, Fractional a)
                      => Int
                      -> [Poly a]
                      -> [Poly a]
implicitizePolynomial numParams paramPolys =
  let m = length paramPolys
      elimOrder = BlockOrder numParams Lex DegRevLex
      -- Coordinate variable x_i is mapped to variable ID (numParams + i)
      coordVar i = Poly [Term 1 (sparseVar (numParams + i) 1)] elimOrder
      -- Ideals relations: x_i - f_i(t) = 0
      relations = [ subPoly (coordVar i) (reorderPoly elimOrder (paramPolys !! i))
                  | i <- [0 .. m - 1]
                  ]
      gb = reducedGrobnerBasis (buchberger relations)
      paramVars = [0 .. numParams - 1]
      -- Retain polynomials containing only coordinate variables
      eliminated = [ reorderPoly DegRevLex (shiftVarsDown numParams p)
                   | p <- gb
                   , freeOfVariables paramVars p
                   , not (isZeroPoly p)
                   ]
  in eliminated
  where
    shiftVarsDown shift (Poly terms order) =
      let shifted = [ Term c (shiftMonomial shift m) | Term c m <- terms ]
      in Poly shifted order

    shiftMonomial shift (SparseMonomial m) =
      let shiftedEntries = [ (v - shift, e) | (v, e) <- Map.toList m ]
      in SparseMonomial (Map.fromList shiftedEntries)

-- | Implicitization of rational parametrizations:
-- Given rational functions [(p_0, q_0), ..., (p_{m-1}, q_{m-1})] representing
-- x_i = p_i(t) / q_i(t), clears denominators via q_i * x_i - p_i = 0
-- and adds the Rabinowitsch variable w to eliminate the denominator locus:
-- 1 - w * (q_0 * ... * q_{m-1}) = 0.
implicitizeRational :: forall a. (Ord a, Fractional a)
                    => Int
                    -> [(Poly a, Poly a)]
                    -> [Poly a]
implicitizeRational numParams rationalParams =
  let m = length rationalParams
      wVar = numParams + m
      elimOrder = BlockOrder (numParams + 1) Lex DegRevLex
      coordVar i = Poly [Term 1 (sparseVar (numParams + i) 1)] elimOrder
      wPoly = Poly [Term 1 (sparseVar wVar 1)] elimOrder
      onePoly = Poly [Term 1 sparseOne] elimOrder

      -- Relations: q_i(t) * x_i - p_i(t) = 0
      denomRelations = [ subPoly (mulPoly (reorderPoly elimOrder qI) (coordVar i))
                                 (reorderPoly elimOrder pI)
                       | (i, (pI, qI)) <- zip [0 :: Int ..] rationalParams
                       ]
      -- Common denominator product Q(t) = q_0 * ... * q_{m-1}
      commonDenom = foldl mulPoly onePoly [ reorderPoly elimOrder q | (_, q) <- rationalParams ]
      rabinowitschRel = subPoly onePoly (mulPoly wPoly commonDenom)

      allRels = rabinowitschRel : denomRelations
      gb = reducedGrobnerBasis (buchberger allRels)
      elimVars = [0 .. numParams - 1] ++ [wVar]

      eliminated = [ reorderPoly DegRevLex (shiftVarsDown numParams p)
                   | p <- gb
                   , freeOfVariables elimVars p
                   , not (isZeroPoly p)
                   ]
  in eliminated
  where
    shiftVarsDown shift (Poly terms order) =
      let shifted = [ Term c (shiftMonomial shift m) | Term c m <- terms ]
      in Poly shifted order

    shiftMonomial shift (SparseMonomial m) =
      let shiftedEntries = [ (v - shift, e) | (v, e) <- Map.toList m ]
      in SparseMonomial (Map.fromList shiftedEntries)

-- | Classic twisted cubic curve in A^3:
-- x = t, y = t^2, z = t^3.
-- Implicit ideal is <y - x^2, z - x*y> (or equivalently <y - x^2, z - x^3>).
twistedCubicIdeal :: (Ord a, Fractional a) => Ideal a
twistedCubicIdeal =
  let t = Poly [Term 1 (sparseVar 0 1)] DegRevLex
      t2 = Poly [Term 1 (sparseVar 0 2)] DegRevLex
      t3 = Poly [Term 1 (sparseVar 0 3)] DegRevLex
      implicitPolys = implicitizePolynomial 1 [t, t2, t3]
  in idealFromGenerators DegRevLex implicitPolys

-- | Whitney umbrella surface in A^3:
-- x = u * v, y = v, z = u^2.
-- Implicit equation is x^2 - y^2 * z = 0.
whitneyUmbrellaIdeal :: (Ord a, Fractional a) => Ideal a
whitneyUmbrellaIdeal =
  let u = Poly [Term 1 (sparseVar 0 1)] DegRevLex
      v = Poly [Term 1 (sparseVar 1 1)] DegRevLex
      u2 = Poly [Term 1 (sparseVar 0 2)] DegRevLex
      uv = mulPoly u v
      implicitPolys = implicitizePolynomial 2 [uv, v, u2]
  in idealFromGenerators DegRevLex implicitPolys

-- | Enneper minimal surface in A^3:
-- x = u - u^3/3 + u*v^2, y = v - v^3/3 + v*u^2, z = u^2 - v^2.
enneperSurfaceIdeal :: forall a. (Ord a, Fractional a) => Ideal a
enneperSurfaceIdeal =
  let u = Poly [Term 1 (sparseVar 0 1)] DegRevLex
      v = Poly [Term 1 (sparseVar 1 1)] DegRevLex
      u2 = Poly [Term 1 (sparseVar 0 2)] DegRevLex
      v2 = Poly [Term 1 (sparseVar 1 2)] DegRevLex
      u3 = Poly [Term 1 (sparseVar 0 3)] DegRevLex
      v3 = Poly [Term 1 (sparseVar 1 3)] DegRevLex

      xParam = subPoly (subPoly u (scalePoly (1/3) u3)) (mulPoly u v2)
      yParam = subPoly (subPoly v (scalePoly (1/3) v3)) (mulPoly v u2)
      zParam = subPoly u2 v2
      implicitPolys = implicitizePolynomial 2 [xParam, yParam, zParam]
  in idealFromGenerators DegRevLex implicitPolys
  where
    scalePoly c (Poly terms order) =
      Poly [ Term (c * coeff) m | Term coeff m <- terms ] order
