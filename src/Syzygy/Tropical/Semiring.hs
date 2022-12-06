{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Tropical.Semiring
Description : Min-plus and max-plus idempotent semirings for tropical geometry
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

The tropical semiring T = (R union {+infinity}, oplus, odot) linearizes
algebraic geometry via the logarithm with base tending to infinity (Maslov dequantization):

  a oplus b = min(a, b)          (Tropical addition, identity: +infinity)
  a odot b  = a + b              (Tropical multiplication, identity: 0)

Properties:
  - Idempotent addition: a oplus a = a
  - Distributivity: a odot (b oplus c) = (a odot b) oplus (a odot c)
  - Annihilator: (+infinity) odot a = +infinity
  - Tropical power: a^{odot n} = n * a

Also provides the dual max-plus semiring (R union {negative infinity}, max, +)
frequently used in discrete event systems and optimal control.
-}
module Syzygy.Tropical.Semiring
  ( Tropical(..)
  , tropZero
  , tropOne
  , tropPlus
  , tropTimes
  , tropPower
  , isTropZero
  , fromTropical
  , MaxPlus(..)
  , maxZero
  , maxOne
  , maxPlus
  , maxTimes
  , maxPower
  , tropDotProduct
  , tropMatMul
  , tropMatVecMul
  ) where

import GHC.Generics (Generic)

-- | Min-plus tropical number: either a finite real value or positive infinity.
data Tropical a
  = TropicalVal !a
  | TropicalInfinity
  deriving (Eq, Show, Read, Generic)

-- | Ordering on tropical numbers: finite numbers compare normally,
-- while TropicalInfinity is strictly greater than all finite values.
instance Ord a => Ord (Tropical a) where
  compare (TropicalVal x) (TropicalVal y) = compare x y
  compare (TropicalVal _) TropicalInfinity = LT
  compare TropicalInfinity (TropicalVal _) = GT
  compare TropicalInfinity TropicalInfinity = EQ

-- | Additive identity: 0_T = +infinity.
tropZero :: Tropical a
tropZero = TropicalInfinity

-- | Multiplicative identity: 1_T = 0.
tropOne :: Num a => Tropical a
tropOne = TropicalVal 0

-- | Tropical addition: a oplus b = min(a, b).
tropPlus :: Ord a => Tropical a -> Tropical a -> Tropical a
tropPlus TropicalInfinity y = y
tropPlus x TropicalInfinity = x
tropPlus (TropicalVal x) (TropicalVal y) = TropicalVal (min x y)

-- | Tropical multiplication: a odot b = a + b.
-- TropicalInfinity is absorbing: infinity odot a = infinity.
tropTimes :: Num a => Tropical a -> Tropical a -> Tropical a
tropTimes TropicalInfinity _ = TropicalInfinity
tropTimes _ TropicalInfinity = TropicalInfinity
tropTimes (TropicalVal x) (TropicalVal y) = TropicalVal (x + y)

-- | Tropical power: a^{odot n} = n * a.
tropPower :: (Integral n, Num a) => Tropical a -> n -> Tropical a
tropPower TropicalInfinity _ = TropicalInfinity
tropPower (TropicalVal x) n = TropicalVal (fromIntegral n * x)

-- | Tests if a tropical number is the additive identity (infinity).
isTropZero :: Tropical a -> Bool
isTropZero TropicalInfinity = True
isTropZero (TropicalVal _) = False

-- | Extracts the underlying value, or returns a default if infinity.
fromTropical :: a -> Tropical a -> a
fromTropical def TropicalInfinity = def
fromTropical _ (TropicalVal x) = x

-- | Max-plus tropical number: either a finite value or negative infinity.
data MaxPlus a
  = MaxVal !a
  | NegInfinity
  deriving (Eq, Show, Read, Generic)

instance Ord a => Ord (MaxPlus a) where
  compare (MaxVal x) (MaxVal y) = compare x y
  compare (MaxVal _) NegInfinity = GT
  compare NegInfinity (MaxVal _) = LT
  compare NegInfinity NegInfinity = EQ

-- | Max-plus additive identity: 0_M = -infinity.
maxZero :: MaxPlus a
maxZero = NegInfinity

-- | Max-plus multiplicative identity: 1_M = 0.
maxOne :: Num a => MaxPlus a
maxOne = MaxVal 0

-- | Max-plus addition: a oplus b = max(a, b).
maxPlus :: Ord a => MaxPlus a -> MaxPlus a -> MaxPlus a
maxPlus NegInfinity y = y
maxPlus x NegInfinity = x
maxPlus (MaxVal x) (MaxVal y) = MaxVal (max x y)

-- | Max-plus multiplication: a odot b = a + b.
maxTimes :: Num a => MaxPlus a -> MaxPlus a -> MaxPlus a
maxTimes NegInfinity _ = NegInfinity
maxTimes _ NegInfinity = NegInfinity
maxTimes (MaxVal x) (MaxVal y) = MaxVal (x + y)

-- | Max-plus power: a^{odot n} = n * a.
maxPower :: (Integral n, Num a) => MaxPlus a -> n -> MaxPlus a
maxPower NegInfinity _ = NegInfinity
maxPower (MaxVal x) n = MaxVal (fromIntegral n * x)

-- | Tropical inner product of two vectors:
-- <u, v>_T = oplus_{i=1}^n (u_i odot v_i) = min_{i=1}^n (u_i + v_i)
tropDotProduct :: (Ord a, Num a) => [Tropical a] -> [Tropical a] -> Tropical a
tropDotProduct u v = foldl tropPlus tropZero (zipWith tropTimes u v)

-- | Tropical matrix-vector multiplication:
-- (M odot v)_i = oplus_j (M_{ij} odot v_j)
tropMatVecMul :: (Ord a, Num a) => [[Tropical a]] -> [Tropical a] -> [Tropical a]
tropMatVecMul rows v = [ tropDotProduct row v | row <- rows ]

-- | Tropical matrix-matrix multiplication:
-- (A odot B)_{ik} = oplus_j (A_{ij} odot B_{jk})
-- Used for all-pairs shortest paths and Kleene star closures.
tropMatMul :: (Ord a, Num a) => [[Tropical a]] -> [[Tropical a]] -> [[Tropical a]]
tropMatMul a b =
  let numColsB = case b of [] -> 0; (row:_) -> length row
      colsB = [ [ (b !! r) !! c | r <- [0 .. length b - 1] ] | c <- [0 .. numColsB - 1] ]
  in [ [ tropDotProduct row col | col <- colsB ] | row <- a ]
