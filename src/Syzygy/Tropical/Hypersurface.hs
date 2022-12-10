{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Tropical.Hypersurface
Description : Tropical polynomials, corner locus hypersurfaces, and Newton polytopes
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Tropical geometry studies the piecewise-linear shadows of algebraic varieties.
A tropical polynomial in n variables w = (w_0, ..., w_{n-1}) in R^n is defined by:

  F(w) = min_{alpha in A} (c_alpha + sum_{i=0}^{n-1} alpha_i * w_i)

The tropical hypersurface T(F) (also called the tropical variety V(F)) is the
corner locus where the minimum is achieved by at least two distinct terms:

  T(F) = { w in R^n | exists alpha /= beta such that
                      c_alpha + <alpha, w> = c_beta + <beta, w> = F(w) }

By Kapranov's Theorem (2000), T(F) equals the non-Archimedean amoeba of the
hypersurface V(f) over a valued field (such as Puiseux series with t-adic valuation).
The polyhedral fan of T(F) is dual to the regular subdivision of the Newton polytope.
-}
module Syzygy.Tropical.Hypersurface
  ( TropTerm(..)
  , TropPoly(..)
  , tropicalizeWithValuation
  , trivialTropicalize
  , evalTropPoly
  , evalTropTerms
  , isHypersurfacePoint
  , supportPoints
  , convexHull2D
  , newtonPolygonArea2D
  ) where

import GHC.Generics (Generic)
import Data.List (sortBy)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse (toDenseList)
import Syzygy.Polynomial.Base (Poly(..), Term(..))

-- | A term in a tropical polynomial: c_alpha + sum_{i=0}^{n-1} alpha_i * w_i.
data TropTerm a = TropTerm
  { ttCoeff :: !a
  , ttExponents :: ![Int]
  } deriving (Eq, Show, Read, Generic)

-- | A multivariate tropical polynomial representing a piecewise-linear concave function.
newtype TropPoly a = TropPoly { unTropPoly :: [TropTerm a] }
  deriving (Eq, Show, Read, Generic)

-- | Tropicalizes a classical multivariate polynomial using a custom valuation val: k* -> Rational.
tropicalizeWithValuation :: (Eq a, Num a) => (a -> Rational) -> Int -> Poly a -> TropPoly Rational
tropicalizeWithValuation val numVars (Poly terms _) =
  let toTropTerm (Term c sm) =
        let exps = toDenseList numVars sm
            v = val c
        in TropTerm v exps
      nonZeroTerms = filter (\(Term c _) -> c /= 0) terms
  in TropPoly (map toTropTerm nonZeroTerms)

-- | Tropicalizes a polynomial under the trivial valuation:
-- val(c) = 0 for all c /= 0.
trivialTropicalize :: (Eq a, Num a) => Int -> Poly a -> TropPoly Rational
trivialTropicalize numVars poly =
  tropicalizeWithValuation (\_ -> 0) numVars poly

-- | Evaluates a tropical polynomial at a coordinate point w in R^n:
-- F(w) = min_{term in F} (c_alpha + sum_{i=0}^{n-1} alpha_i * w_i).
evalTropPoly :: forall a. (Ord a, Num a) => TropPoly a -> [a] -> Maybe a
evalTropPoly (TropPoly terms) point =
  case terms of
    [] -> Nothing
    _  ->
      let evalTerm (TropTerm c exps) =
            let dot = sum (zipWith (\e x -> fromIntegral e * x) exps point)
            in c + dot
      in Just (minimum (map evalTerm terms))

-- | Evaluates each term individually, returning pairs of (evaluated_value, exponent_vector).
evalTropTerms :: (Ord a, Num a) => TropPoly a -> [a] -> [(a, [Int])]
evalTropTerms (TropPoly terms) point =
  [ (c + sum (zipWith (\e x -> fromIntegral e * x) exps point), exps)
  | TropTerm c exps <- terms
  ]

-- | Determines if a point w lies on the tropical hypersurface T(F).
-- A point lies on T(F) if the minimum of F(w) is attained at least twice
-- (within a numerical or exact tolerance eps).
isHypersurfacePoint :: (Ord a, Num a) => TropPoly a -> [a] -> a -> Bool
isHypersurfacePoint tPoly point eps =
  let termVals = evalTropTerms tPoly point
  in case termVals of
       [] -> False
       (_:_) ->
         let minVal = minimum (map fst termVals)
             -- Count how many terms achieve within eps of the minimum
             achievedCount = length (filter (\(v, _) -> abs (v - minVal) <= eps) termVals)
         in achievedCount >= 2

-- | Extracts the support (set of distinct exponent vectors) of a polynomial.
supportPoints :: Int -> Poly a -> [[Int]]
supportPoints numVars (Poly terms _) =
  let exps = [ toDenseList numVars sm | Term _ sm <- terms ]
      uniqueExps = Map.keys (Map.fromList [ (e, ()) | e <- exps ])
  in uniqueExps

-- | Computes the 2D convex hull of integer lattice points using Andrew's monotone chain algorithm.
-- Used to construct Newton polygons in dimension 2.
convexHull2D :: (Num a, Ord a) => [(a, a)] -> [(a, a)]
convexHull2D pts =
  let sortedPts = sortBy comparePts (removeDuplicates pts)
      cross (oX, oY) (aX, aY) (bX, bY) =
        (aX - oX) * (bY - oY) - (aY - oY) * (bX - oX)

      -- Lower hull
      lower = foldl (\hull p -> buildHull hull p) [] sortedPts
      -- Upper hull
      upper = foldl (\hull p -> buildHull hull p) [] (reverse sortedPts)

      buildHull [] p = [p]
      buildHull [p1] p = [p, p1]
      buildHull (p2 : p1 : rest) p =
        if cross p1 p2 p <= 0
        then buildHull (p1 : rest) p
        else p : p2 : p1 : rest

      comparePts (x1, y1) (x2, y2) =
        case compare x1 x2 of
          EQ -> compare y1 y2
          other -> other

      removeDuplicates [] = []
      removeDuplicates (x:xs) = x : removeDuplicates (filter (/= x) xs)
  in case sortedPts of
       [] -> []
       [_] -> sortedPts
       _ -> reverse (drop 1 lower) ++ reverse (drop 1 upper)

-- | Computes the area of the Newton polygon in 2D using the Shoelace formula:
-- 2 * Area = sum_{i=0}^{n-1} (x_i * y_{i+1} - x_{i+1} * y_i).
newtonPolygonArea2D :: Fractional a => [(a, a)] -> a
newtonPolygonArea2D vertices =
  let n = length vertices
  in if n < 3
     then 0
     else
       case vertices of
         (firstV : _) ->
           let pairs = zip vertices (drop 1 vertices ++ [firstV])
               crossSum = sum [ x1 * y2 - x2 * y1 | ((x1, y1), (x2, y2)) <- pairs ]
           in abs crossSum / 2
         [] -> 0
