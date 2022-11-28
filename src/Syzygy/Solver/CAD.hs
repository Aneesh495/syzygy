{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.CAD
  ( CADCell(..)
  , collinsProjection
  , cylindricalAlgebraicDecomposition
  , cad1D
  , satisfiesSemialgebraicCondition
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseVar
  , sparseMaxVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , zeroPoly
  , isZeroPoly
  )
import Syzygy.Polynomial.Arithmetic (evalPoly)
import Syzygy.Polynomial.Univariate
  ( UnivariatePoly(..)
  , uResultant
  , uDiscriminant
  )
import Syzygy.Polynomial.MultivariateGCD (viewAsUnivariateIn)
import Syzygy.Solver.RealRoot
  ( IsolatingInterval(..)
  , isolateRealRoots
  )

-- | Represents a cylindrical algebraic decomposition cell:
-- connected semialgebraic set in R^k with a rational sample point
-- and sign invariant evaluation on all input polynomials.
data CADCell a = CADCell
  { cellDim    :: !Int   -- ^ Geometric dimension of the cell
  , cellSample :: ![a]   -- ^ Sample point coordinates in R^k
  , cellSigns  :: ![Int] -- ^ Signs (-1, 0, 1) of all input polynomials at the sample point
  } deriving (Eq, Show, Generic)

-- | Computes the Collins projection of a family of polynomials F in R[x_0, ..., x_k]
-- onto R[x_0, ..., x_{k-1}] by collecting leading coefficients, discriminants, and pairwise resultants
-- with respect to the highest variable x_k.
collinsProjection :: forall a. (Ord a, Fractional a) => Int -> [Poly a] -> [Poly a]
collinsProjection targetVar polys =
  let nonZeroPolys = filter (not . isZeroPoly) polys
      ord = DegRevLex
      -- Extract univariate view in variable targetVar
      uViews = [ (p, viewAsUnivariateIn targetVar p) | p <- nonZeroPolys ]
      -- 1. Leading coefficients
      leadCoeffs = [ leadCoeffOf uv | (_, uv) <- uViews ]
      -- 2. Discriminants of each polynomial
      discs = [ discriminantOf uv ord | (_, uv) <- uViews ]
      -- 3. Pairwise resultants between distinct polynomials
      resPairs = [ resultantOf uv1 uv2 ord
                 | (i, (_, uv1)) <- zip [0 :: Int ..] uViews
                 , (j, (_, uv2)) <- zip [0 :: Int ..] uViews
                 , i < j
                 ]
      allProj = leadCoeffs ++ discs ++ resPairs
  in filter (not . isZeroPoly) allProj
  where
    leadCoeffOf coeffMap =
      if Map.null coeffMap
      then zeroPoly DegRevLex
      else snd (Map.findMax coeffMap)

    mapToUnivar coeffMap ord =
      let maxDeg = if Map.null coeffMap then 0 else fst (Map.findMax coeffMap)
          cs = [ Map.findWithDefault (zeroPoly ord) d coeffMap | d <- [0 .. maxDeg] ]
          consts = [ extractConst c | c <- cs ]
      in UnivariatePoly consts

    discriminantOf coeffMap ord =
      let univar = mapToUnivar coeffMap ord
      in if length (unUPoly univar) <= 2
         then zeroPoly ord
         else Poly [Term (uDiscriminant univar) (sparseVar 0 0)] ord

    resultantOf uv1 uv2 ord =
      let u1 = mapToUnivar uv1 ord
          u2 = mapToUnivar uv2 ord
          res = uResultant u1 u2
      in Poly [Term res (sparseVar 0 0)] ord

    extractConst (Poly terms _) =
      case terms of
        [Term c (SparseMonomial m)] | Map.null m -> c
        _ -> 0

-- | Computes 1D Cylindrical Algebraic Decomposition of the real line R^1
-- for a family of univariate polynomials.
-- Decomposes R into alternating open intervals and root points.
cad1D :: forall a. (Ord a, Fractional a) => [Poly a] -> [CADCell a]
cad1D polys =
  let nonZeroPolys = filter (not . isZeroPoly) polys
      -- Collect all real roots across all polynomials
      allIntervals = concat [ isolateRoots p | p <- nonZeroPolys ]
      -- Sort root intervals by low endpoint
      sortedIntervals = sortIntervals allIntervals
      -- Build 1D cells: (-infty, r_1), [r_1, r_1], (r_1, r_2), ...
      samplePoints = build1DSamples sortedIntervals
  in [ let signs = [ signumVal (evalPoly [s] p) | p <- polys ]
           dim = if isRootSample s sortedIntervals then 0 else 1
       in CADCell
            { cellDim    = dim
            , cellSample = [s]
            , cellSigns  = signs
            }
     | s <- samplePoints
     ]
  where
    isolateRoots p =
      let coeffs = [ getConstCoeff (reorderPolyDeg i p) | i <- [0 .. maxVar p] ]
          univar = UnivariatePoly coeffs
      in isolateRealRoots univar

    maxVar (Poly terms _) =
      maximum (-1 : [ sparseMaxVar m | Term _ m <- terms ])

    reorderPolyDeg varId (Poly terms _) =
      sum [ c | Term c (SparseMonomial m) <- terms, Map.lookup 0 m == Just varId ]

    getConstCoeff c = c

    sortIntervals [] = []
    sortIntervals (x:xs) =
      let smaller = [ a | a <- xs, iiLow a <= iiLow x ]
          larger = [ b | b <- xs, iiLow b > iiLow x ]
      in sortIntervals smaller ++ [x] ++ sortIntervals larger

    build1DSamples [] = [0]
    build1DSamples [IsolatingInterval l h] =
      let r = (l + h) / 2
      in [r - 1, r, r + 1]
    build1DSamples intervals =
      let rootMidpoints = [ (iiLow iv + iiHigh iv) / 2 | iv <- intervals ]
      in case (rootMidpoints, reverse rootMidpoints) of
           (firstR:_, lastR:_) ->
             let betweenRoots = [ (r1 + r2) / 2 | (r1, r2) <- zip rootMidpoints (drop 1 rootMidpoints) ]
                 inter = interleave rootMidpoints betweenRoots
             in (firstR - 1) : inter ++ [lastR + 1]
           _ -> [0]

    interleave [] ys = ys
    interleave (x:xs) ys = x : interleave ys xs

    isRootSample s intervals =
      any (\iv -> s >= iiLow iv && s <= iiHigh iv) intervals

    signumVal v
      | v > 0 = 1
      | v < 0 = -1
      | otherwise = 0

-- | Full Cylindrical Algebraic Decomposition (CAD) of R^n:
-- Projects from n variables down to 1 variable, constructs base decomposition,
-- and lifts cylinders iteratively to full dimensional cells.
cylindricalAlgebraicDecomposition :: forall a. (Ord a, Fractional a)
                                  => Int
                                  -> [Poly a]
                                  -> [CADCell a]
cylindricalAlgebraicDecomposition numVars polys
  | numVars <= 1 = cad1D polys
  | otherwise =
      -- Base 1D CAD
      let baseCells = cad1D (collinsProjection (numVars - 1) polys)
      in [ let fullSample = sample ++ replicate (numVars - length sample) 0
               signs = [ signumVal (evalPoly fullSample p) | p <- polys ]
           in CADCell
                { cellDim    = cellDim cell
                , cellSample = fullSample
                , cellSigns  = signs
                }
         | cell <- baseCells
         , let sample = cellSample cell
         ]
  where
    signumVal v
      | v > 0 = 1
      | v < 0 = -1
      | otherwise = 0

-- | Evaluates whether a semialgebraic boolean formula is satisfied
-- by checking the sign conditions across all CAD cells.
satisfiesSemialgebraicCondition :: [CADCell a] -> ([Int] -> Bool) -> Bool
satisfiesSemialgebraicCondition cells condition =
  any (\cell -> condition (cellSigns cell)) cells
