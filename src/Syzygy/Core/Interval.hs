{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Core.Interval
  ( Interval(..)
  , interval
  , pointInterval
  , intervalLower
  , intervalUpper
  , intervalMid
  , intervalWidth
  , intervalRadius
  , intervalIntersect
  , intervalHull
  , inInterval
  , isSubsetOf
  , isDisjoint
  , intervalPow
  , intervalNewtonStep
  , bisectInterval
  ) where

import GHC.Generics (Generic)
import Syzygy.Core.Rational (Rat(..), rat)

-- | Exact closed interval [lower, upper] with rational endpoints where lower <= upper.
data Interval = Interval {-# UNPACK #-} !Rat {-# UNPACK #-} !Rat
  deriving (Eq, Generic)

-- | Smart constructor for closed interval [lo, hi].
-- Enforces lo <= hi by swapping if necessary.
interval :: Rat -> Rat -> Interval
interval lo hi
  | lo <= hi  = Interval lo hi
  | otherwise = Interval hi lo

-- | Degenerate point interval [x, x].
pointInterval :: Rat -> Interval
pointInterval x = Interval x x

-- | Lower bound of the interval.
intervalLower :: Interval -> Rat
intervalLower (Interval lo _) = lo

-- | Upper bound of the interval.
intervalUpper :: Interval -> Rat
intervalUpper (Interval _ hi) = hi

-- | Midpoint of the interval: (lower + upper) / 2.
intervalMid :: Interval -> Rat
intervalMid (Interval lo hi) = (lo + hi) / 2

-- | Width of the interval: upper - lower.
intervalWidth :: Interval -> Rat
intervalWidth (Interval lo hi) = hi - lo

-- | Radius of the interval: (upper - lower) / 2.
intervalRadius :: Interval -> Rat
intervalRadius (Interval lo hi) = (hi - lo) / 2

-- | Tests if a rational point x lies inside the interval [lo, hi].
inInterval :: Rat -> Interval -> Bool
inInterval x (Interval lo hi) = lo <= x && x <= hi

-- | Tests if the first interval is completely contained within the second interval.
isSubsetOf :: Interval -> Interval -> Bool
isSubsetOf (Interval l1 h1) (Interval l2 h2) = l2 <= l1 && h1 <= h2

-- | Tests if two intervals are strictly disjoint (empty intersection).
isDisjoint :: Interval -> Interval -> Bool
isDisjoint (Interval l1 h1) (Interval l2 h2) = h1 < l2 || h2 < l1

-- | Exact intersection of two intervals.
-- Returns Just [max(l1, l2), min(h1, h2)] if they overlap, or Nothing if disjoint.
intervalIntersect :: Interval -> Interval -> Maybe Interval
intervalIntersect (Interval l1 h1) (Interval l2 h2)
  | maxLo <= minHi = Just (Interval maxLo minHi)
  | otherwise      = Nothing
  where
    maxLo = max l1 l2
    minHi = min h1 h2

-- | Interval hull (the smallest single interval containing both intervals).
intervalHull :: Interval -> Interval -> Interval
intervalHull (Interval l1 h1) (Interval l2 h2) =
  Interval (min l1 l2) (max h1 h2)

-- | Exact non-negative integer power of an interval.
intervalPow :: Interval -> Integer -> Interval
intervalPow _ 0 = Interval 1 1
intervalPow (Interval lo hi) p
  | p < 0 = error "intervalPow: negative power requires non-zero interval"
  | odd p = Interval (lo ^ p) (hi ^ p)
  | lo >= 0 = Interval (lo ^ p) (hi ^ p)
  | hi <= 0 = Interval (hi ^ p) (lo ^ p)
  | otherwise =
      let maxVal = max (abs lo ^ p) (hi ^ p)
      in Interval 0 maxVal

-- | Bisects an interval into two equal halves at the midpoint.
bisectInterval :: Interval -> (Interval, Interval)
bisectInterval iv@(Interval lo hi) =
  let m = intervalMid iv
  in (Interval lo m, Interval m hi)

-- | Interval Newton step: N(X) = m - f(m) / f'(X) where m = mid(X).
-- If N(X) is a subset of X, then f has a unique real root in X by Moore's theorem.
intervalNewtonStep :: (Rat -> Rat) -> (Interval -> Interval) -> Interval -> Maybe Interval
intervalNewtonStep f fPrime iv@(Interval _ _) =
  let m = intervalMid iv
      fm = f m
      fPrimeX = fPrime iv
  in if inInterval 0 fPrimeX
     then Nothing
     else
       let step = pointInterval fm / fPrimeX
           nX = pointInterval m - step
       in intervalIntersect iv nX

instance Show Interval where
  show (Interval lo hi) = "[" ++ show lo ++ ", " ++ show hi ++ "]"

instance Num Interval where
  (Interval l1 h1) + (Interval l2 h2) =
    Interval (l1 + l2) (h1 + h2)

  (Interval l1 h1) - (Interval l2 h2) =
    Interval (l1 - h2) (h1 - l2)

  (Interval l1 h1) * (Interval l2 h2) =
    let p1 = l1 * l2
        p2 = l1 * h2
        p3 = h1 * l2
        p4 = h1 * h2
        minP = min (min p1 p2) (min p3 p4)
        maxP = max (max p1 p2) (max p3 p4)
    in Interval minP maxP

  negate (Interval lo hi) = Interval (-hi) (-lo)

  abs (Interval lo hi)
    | lo >= 0   = Interval lo hi
    | hi <= 0   = Interval (-hi) (-lo)
    | otherwise = Interval 0 (max (-lo) hi)

  signum (Interval lo hi)
    | lo > 0    = Interval 1 1
    | hi < 0    = Interval (-1) (-1)
    | lo == 0 && hi == 0 = Interval 0 0
    | otherwise = Interval (-1) 1

  fromInteger n =
    let r = rat n 1
    in Interval r r

instance Fractional Interval where
  fromRational r =
    let q = fromRational r
    in Interval q q

  (Interval l1 h1) / (Interval l2 h2)
    | l2 <= 0 && h2 >= 0 = error "Interval division by an interval containing zero"
    | otherwise =
        let q1 = l1 / l2
            q2 = l1 / h2
            q3 = h1 / l2
            q4 = h1 / h2
            minQ = min (min q1 q2) (min q3 q4)
            maxQ = max (max q1 q2) (max q3 q4)
        in Interval minQ maxQ

  recip (Interval lo hi)
    | lo <= 0 && hi >= 0 = error "Interval recip: interval contains zero"
    | otherwise = Interval (1 / hi) (1 / lo)
