{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.RealRoot
  ( IsolatingInterval(..)
  , sturmSequence
  , countSignVariations
  , realRootsCountInInterval
  , totalRealRootsCount
  , cauchyBound
  , isolateRealRoots
  , refineInterval
  , uEval
  ) where

import GHC.Generics (Generic)
import Syzygy.Polynomial.Univariate
  ( UnivariatePoly(..)
  , uPoly
  , uDegree
  , uDerivative
  , uQuotRem
  , uMul
  , uSquareFree
  )

-- | Isolating interval [low, high] guaranteed to contain exactly one real root.
data IsolatingInterval a = IsolatingInterval
  { iiLow  :: !a
  , iiHigh :: !a
  } deriving (Eq, Show, Generic)

-- | Evaluates a univariate polynomial at a scalar value x via Horner's rule.
uEval :: Num c => UnivariatePoly c -> c -> c
uEval (UnivariatePoly cs) x = foldr (\c acc -> c + x * acc) 0 cs

-- | Negates a univariate polynomial.
uNeg :: Num c => UnivariatePoly c -> UnivariatePoly c
uNeg (UnivariatePoly cs) = UnivariatePoly (map negate cs)

-- | Square-free part of a univariate polynomial: product of its square-free factors.
uSquareFreePart :: (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c
uSquareFreePart p =
  case uSquareFree p of
    [] -> p
    factors -> foldl uMul (uPoly [1]) factors

-- | Constructs the classic Sturm sequence for a univariate polynomial f:
-- f_0 = f, f_1 = f', f_{k+1} = -rem(f_{k-1}, f_k).
sturmSequence :: forall a. (Eq a, Fractional a) => UnivariatePoly a -> [UnivariatePoly a]
sturmSequence p =
  let p0 = p
      p1 = uDerivative p
  in if uDegree p1 < 0
     then [p0]
     else p0 : p1 : buildSturm p0 p1
  where
    buildSturm f0 f1 =
      let (_, remPoly) = uQuotRem f0 f1
          negRem = uNeg remPoly
      in if uDegree negRem < 0
         then []
         else negRem : buildSturm f1 negRem

-- | Counts the number of sign variations in a sequence of scalar values, ignoring zeros:
-- e.g. [+, -, +, +] has 2 variations.
countSignVariations :: (Ord a, Num a) => [a] -> Int
countSignVariations vals =
  let nonZeros = filter (/= 0) vals
  in case nonZeros of
       [] -> 0
       [_] -> 0
       (x:xs) ->
         let (_, count) = foldl (\(prev, acc) curr ->
                                  if (prev > 0 && curr < 0) || (prev < 0 && curr > 0)
                                  then (curr, acc + 1)
                                  else (curr, acc)
                                ) (x, 0) xs
         in count

-- | By Sturm's Theorem (1829):
-- The exact number of real roots of f in the half-open interval (a, b] is:
-- N(a, b] = V(a) - V(b)
-- where V(x) is the number of sign variations in the Sturm sequence evaluated at x.
realRootsCountInInterval :: forall a. (Ord a, Fractional a)
                         => [UnivariatePoly a]
                         -> a
                         -> a
                         -> Int
realRootsCountInInterval sturmSeq a b =
  let valsA = [ uEval p a | p <- sturmSeq ]
      valsB = [ uEval p b | p <- sturmSeq ]
      vA = countSignVariations valsA
      vB = countSignVariations valsB
  in max 0 (vA - vB)

-- | Computes the Cauchy root bound:
-- all real (and complex) roots of f(x) = sum_{i=0}^n a_i * x^i satisfy |z| <= 1 + max_{i < n} |a_i / a_n|.
cauchyBound :: forall a. (Ord a, Fractional a) => UnivariatePoly a -> a
cauchyBound (UnivariatePoly cs) =
  let trimmed = dropWhile (== 0) (reverse cs)
  in case trimmed of
       [] -> 1
       [_] -> 1
       (lead:restRev) ->
         let maxRatio = maximum [ abs c / abs lead | c <- restRev ]
         in 1 + maxRatio

-- | Returns the total count of real roots in (-infty, infty).
totalRealRootsCount :: (Ord a, Fractional a) => UnivariatePoly a -> Int
totalRealRootsCount p =
  let sqFree = uSquareFreePart p
      sturm = sturmSequence sqFree
      bound = cauchyBound sqFree
  in realRootsCountInInterval sturm (-bound - 1) (bound + 1)

-- | Isolates all real roots of a univariate polynomial into disjoint isolating intervals [low, high]
-- each containing exactly one real root, using Sturm bisection.
isolateRealRoots :: forall a. (Ord a, Fractional a)
                 => UnivariatePoly a
                 -> [IsolatingInterval a]
isolateRealRoots p =
  let sqFree = uSquareFreePart p
      sturm = sturmSequence sqFree
      bound = cauchyBound sqFree
      initLow = -bound - 1
      initHigh = bound + 1
  in bisect sturm initLow initHigh
  where
    bisect sturm low high =
      let nRoots = realRootsCountInInterval sturm low high
      in if nRoots == 0
         then []
         else if nRoots == 1
              then [IsolatingInterval low high]
              else let mid = (low + high) / 2
                   in bisect sturm low mid ++ bisect sturm mid high

-- | Refines an isolating interval to width at most epsilon using bisection.
refineInterval :: forall a. (Ord a, Fractional a)
               => UnivariatePoly a
               -> a
               -> IsolatingInterval a
               -> IsolatingInterval a
refineInterval p eps (IsolatingInterval low high)
  | high - low <= eps = IsolatingInterval low high
  | otherwise =
      let mid = (low + high) / 2
          sturm = sturmSequence (uSquareFreePart p)
          leftCount = realRootsCountInInterval sturm low mid
      in if leftCount > 0
         then refineInterval p eps (IsolatingInterval low mid)
         else refineInterval p eps (IsolatingInterval mid high)
