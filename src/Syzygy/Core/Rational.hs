{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Core.Rational
  ( Rat(..)
  , rat
  , numerator
  , denominator
  , recipRat
  , approximate
  , continuedFraction
  , fromContinuedFraction
  , convergents
  , toDecimalString
  , mediant
  , powRat
  , floorRat
  , ceilingRat
  , isZeroRat
  , isOneRat
  ) where

import GHC.Generics (Generic)
import qualified Data.Ratio as R
import Syzygy.Core.BigInt (extendedGCD)

-- | Exact arbitrary precision rational number representing numerator / denominator
-- where denominator > 0 and gcd(numerator, denominator) == 1.
data Rat = Rat {-# UNPACK #-} !Integer {-# UNPACK #-} !Integer
  deriving (Eq, Generic)

-- | Smart constructor for Rat enforcing canonical representation:
-- positive denominator and coprime numerator and denominator.
rat :: Integer -> Integer -> Rat
rat _ 0 = error "rat: division by zero denominator"
rat n d
  | d < 0     = rat (-n) (-d)
  | n == 0    = Rat 0 1
  | otherwise =
      let (g, _, _) = extendedGCD (abs n) d
      in Rat (n `div` g) (d `div` g)

-- | Extracts the numerator.
numerator :: Rat -> Integer
numerator (Rat n _) = n

-- | Extracts the denominator (always > 0).
denominator :: Rat -> Integer
denominator (Rat _ d) = d

-- | Exact reciprocal of a rational number.
recipRat :: Rat -> Rat
recipRat (Rat 0 _) = error "recipRat: division by zero"
recipRat (Rat n d)
  | n < 0     = Rat (-d) (-n)
  | otherwise = Rat d n

-- | Exact integer power of a rational number.
powRat :: Rat -> Integer -> Rat
powRat _ 0 = Rat 1 1
powRat (Rat 0 _) p
  | p > 0     = Rat 0 1
  | otherwise = error "powRat: division by zero"
powRat (Rat n d) p
  | p < 0     = powRat (recipRat (Rat n d)) (-p)
  | otherwise = Rat (n ^ p) (d ^ p)

-- | Tests if a rational number is zero.
isZeroRat :: Rat -> Bool
isZeroRat (Rat 0 _) = True
isZeroRat _         = False

-- | Tests if a rational number is one.
isOneRat :: Rat -> Bool
isOneRat (Rat 1 1) = True
isOneRat _         = False

-- | Computes the floor of a rational number.
floorRat :: Rat -> Integer
floorRat (Rat n d) = n `div` d

-- | Computes the ceiling of a rational number.
ceilingRat :: Rat -> Integer
ceilingRat (Rat n d) =
  let (q, r) = n `divMod` d
  in if r == 0 then q else q + 1

-- | Mediant of two rational numbers: (n1 + n2) / (d1 + d2).
mediant :: Rat -> Rat -> Rat
mediant (Rat n1 d1) (Rat n2 d2) = rat (n1 + n2) (d1 + d2)

-- | Computes the exact finite continued fraction expansion of a rational number:
-- [a0; a1, a2, ..., an].
continuedFraction :: Rat -> [Integer]
continuedFraction (Rat n d) = go n d
  where
    go _ 0 = []
    go num den =
      let (q, r) = num `divMod` den
      in q : go den r

-- | Reconstructs a rational number from its continued fraction coefficients.
fromContinuedFraction :: [Integer] -> Rat
fromContinuedFraction [] = Rat 0 1
fromContinuedFraction [a0] = Rat a0 1
fromContinuedFraction (a0:as) =
  Rat a0 1 + recipRat (fromContinuedFraction as)

-- | Computes the list of successive rational convergents p_k / q_k
-- from a continued fraction expansion.
convergents :: [Integer] -> [Rat]
convergents [] = []
convergents (a0:as) = go as (Rat a0 1) (Rat 1 0)
  where
    go [] curr _ = [curr]
    go (x:xs) (Rat p1 q1) (Rat p0 q0) =
      let p2 = x * p1 + p0
          q2 = x * q1 + q0
          nextConv = Rat p2 q2
      in (Rat p1 q1) : go xs nextConv (Rat p1 q1)

-- | Computes the best rational approximation to r with denominator bounded by maxDen
-- using the Stern-Brocot tree / continued fraction convergents.
approximate :: Integer -> Rat -> Rat
approximate maxDen r@(Rat _ d)
  | d <= maxDen = r
  | maxDen <= 0 = error "approximate: maxDen must be positive"
  | otherwise   =
      let cf = continuedFraction r
          convs = convergents cf
          validConvs = filter (\(Rat _ cd) -> cd <= maxDen) convs
      in if null validConvs then Rat (floorRat r) 1 else last validConvs

-- | Formats a rational number as a decimal string up to a specified number of digits.
toDecimalString :: Int -> Rat -> String
toDecimalString precision (Rat n d)
  | d == 1    = show n
  | otherwise =
      let (intPart, remPart) = n `divMod` d
          sign = if n < 0 && intPart == 0 then "-" else ""
          remAbs = abs remPart
          go 0 _ = ""
          go !k !currRem =
            let (digit, nextRem) = (currRem * 10) `divMod` d
            in show digit ++ go (k - 1) nextRem
          fracPart = go precision remAbs
      in sign ++ show intPart ++ "." ++ fracPart

instance Ord Rat where
  compare (Rat n1 d1) (Rat n2 d2) = compare (n1 * d2) (n2 * d1)

instance Show Rat where
  show (Rat n 1) = show n
  show (Rat n d) = show n ++ "/" ++ show d

instance Num Rat where
  (Rat n1 d1) + (Rat n2 d2)
    | d1 == d2  = rat (n1 + n2) d1
    | otherwise = rat (n1 * d2 + n2 * d1) (d1 * d2)

  (Rat n1 d1) - (Rat n2 d2)
    | d1 == d2  = rat (n1 - n2) d1
    | otherwise = rat (n1 * d2 - n2 * d1) (d1 * d2)

  (Rat n1 d1) * (Rat n2 d2) = rat (n1 * n2) (d1 * d2)

  negate (Rat n d) = Rat (-n) d

  abs (Rat n d) = Rat (abs n) d

  signum (Rat n _)
    | n > 0     = Rat 1 1
    | n < 0     = Rat (-1) 1
    | otherwise = Rat 0 1

  fromInteger n = Rat n 1

instance Fractional Rat where
  fromRational r =
    rat (R.numerator r) (R.denominator r)

  (Rat n1 d1) / (Rat n2 d2)
    | n2 == 0   = error "division by zero Rat"
    | otherwise = rat (n1 * d2) (d1 * n2)

  recip = recipRat

instance Real Rat where
  toRational (Rat n d) = n R.% d

instance RealFrac Rat where
  properFraction (Rat n d) =
    let (q, r) = n `quotRem` d
    in (fromInteger q, Rat r d)
