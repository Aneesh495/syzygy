{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Core.FiniteField
  ( FF(..)
  , ff
  , ffCharacteristic
  , ffValue
  , ffRecip
  , ffDiv
  , ffPow
  , batchInverse
  , isQuadraticResidue
  , ffSqrt
  , primitiveRoot
  , discreteLog
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Core.BigInt (extendedGCD, modPow, factorize, integerSqrt)

-- | Element of the finite field F_p where value is canonically 0 <= value < characteristic.
data FF = FF {-# UNPACK #-} !Integer {-# UNPACK #-} !Integer
  deriving (Eq, Generic)

-- | Smart constructor for finite field element val mod p.
ff :: Integer -> Integer -> FF
ff _ p | p <= 1 = error "ff: characteristic must be a prime > 1"
ff val p = FF ((val `mod` p + p) `mod` p) p

-- | Returns the characteristic p of the finite field.
ffCharacteristic :: FF -> Integer
ffCharacteristic (FF _ p) = p

-- | Returns the canonical integer representative 0 <= val < p.
ffValue :: FF -> Integer
ffValue (FF val _) = val

-- | Multiplicative inverse in F_p.
ffRecip :: FF -> FF
ffRecip (FF 0 _) = error "ffRecip: division by zero in finite field"
ffRecip (FF val p) =
  let (g, x, _) = extendedGCD val p
  in if g == 1
     then FF ((x `mod` p + p) `mod` p) p
     else error "ffRecip: non-invertible element (characteristic not prime)"

-- | Division in F_p.
ffDiv :: FF -> FF -> FF
ffDiv a b = a * ffRecip b

-- | Exponentiation in F_p using binary exponentiation.
ffPow :: FF -> Integer -> FF
ffPow (FF _ p) 0 = FF 1 p
ffPow (FF 0 p) expVal
  | expVal > 0 = FF 0 p
  | otherwise  = error "ffPow: 0 raised to non-positive power"
ffPow (FF val p) expVal
  | expVal < 0 =
      let inv = ffRecip (FF val p)
      in ffPow inv (-expVal)
  | otherwise  = FF (modPow val expVal p) p

-- | Batch inversion of a list of finite field elements using Montgomery's trick.
-- Inverts k elements simultaneously with only 1 full modular inversion and 3k multiplications.
batchInverse :: [FF] -> [FF]
batchInverse [] = []
batchInverse (firstEl:rest)
  | any (\(FF v _) -> v == 0) (firstEl:rest) =
      error "batchInverse: contains zero element"
  | otherwise =
      let p = ffCharacteristic firstEl
          vals = map ffValue (firstEl:rest)
          -- Prefix products: p_0 = 1, p_i = p_{i-1} * x_i mod p
          prefixes = scanl (\acc x -> (acc * x) `mod` p) 1 vals
          totalProd = last prefixes
          -- Single inversion of total product
          (g, totalInv, _) = extendedGCD totalProd p
          invAll = if g == 1 then (totalInv `mod` p + p) `mod` p else error "batchInverse: product not coprime"
          -- Backward pass to compute individual inverses
          go [] _ _ = []
          go (x:xs) (pref:prefs) !currInv =
            let invX = (currInv * pref) `mod` p
                nextInv = (currInv * x) `mod` p
            in FF invX p : go xs prefs nextInv
          go _ _ _ = []
      in reverse (go (reverse vals) (reverse (init prefixes)) invAll)

-- | Tests if a non-zero element is a quadratic residue modulo p using Euler's criterion.
isQuadraticResidue :: FF -> Bool
isQuadraticResidue (FF 0 _) = True
isQuadraticResidue (FF val p)
  | p == 2    = True
  | otherwise = modPow val ((p - 1) `div` 2) p == 1

-- | Computes the square root of a finite field element using the Tonelli-Shanks algorithm.
-- Returns Just r where r*r == a, or Nothing if a is a quadratic non-residue.
ffSqrt :: FF -> Maybe FF
ffSqrt (FF 0 p) = Just (FF 0 p)
ffSqrt a@(FF val p)
  | not (isQuadraticResidue a) = Nothing
  | p `mod` 4 == 3 =
      let root = modPow val ((p + 1) `div` 4) p
      in Just (FF root p)
  | otherwise =
      -- Tonelli-Shanks for p - 1 = q * 2^s
      let (q, (s :: Int)) = factorTwos (p - 1) (0 :: Int)
          -- Find a quadratic non-residue z
          z = findNonResidue 2
          findNonResidue candidate
            | candidate >= p = error "ffSqrt: no non-residue found"
            | not (isQuadraticResidue (FF candidate p)) = candidate
            | otherwise = findNonResidue (candidate + 1)
          c = modPow z q p
          rInit = modPow val ((q + 1) `div` 2) p
          tInit = modPow val q p
      in go s c rInit tInit
  where
    factorTwos m !(acc :: Int)
      | even m    = factorTwos (m `div` 2) (acc + 1)
      | otherwise = (m, acc)

    go _ _ !r 1 = Just (FF r p)
    go 0 _ _ _  = Nothing
    go !m !c !r !t =
      let i = findSmallestI 1 (t * t `mod` p)
          findSmallestI !idx !currT
            | idx >= m = idx
            | currT == 1 = idx
            | otherwise = findSmallestI (idx + 1) ((currT * currT) `mod` p)
      in if i >= m
         then Nothing
         else
           let b = modPow c (modPow 2 (toInteger (m - i - 1)) (p - 1)) p
               m' = i
               c' = (b * b) `mod` p
               t' = (t * c') `mod` p
               r' = (r * b) `mod` p
           in go m' c' r' t'

-- | Finds a primitive root modulo prime p (generator of the multiplicative group).
primitiveRoot :: Integer -> Maybe Integer
primitiveRoot p
  | p < 2     = Nothing
  | p == 2    = Just 1
  | otherwise =
      let phi = p - 1
          factors = map fst (factorize phi)
          candidates = [2 .. p - 1]
          isPrim g = all (\q -> modPow g (phi `div` q) p /= 1) factors
      in case filter isPrim candidates of
           (g:_) -> Just g
           []    -> Nothing

-- | Discrete logarithm: finds x such that base^x == target mod p
-- using Shanks' baby-step giant-step algorithm in O(sqrt(p)) time.
discreteLog :: FF -> FF -> Maybe Integer
discreteLog (FF b p) (FF target p2)
  | p /= p2 = error "discreteLog: field characteristics must match"
  | target == 1 = Just 0
  | b == 0 = if target == 0 then Just 1 else Nothing
  | otherwise =
      let m = integerSqrt p + 1
          -- Baby steps: compute b^j mod p for 0 <= j < m and store in a lookup map
          babySteps = goBaby 0 1 Map.empty
          goBaby !j !currVal !accMap
            | j >= m = accMap
            | otherwise =
                let acc' = Map.insert currVal j accMap
                    nextVal = (currVal * b) `mod` p
                in goBaby (j + 1) nextVal acc'
          -- Giant step multiplier: b^(-m) mod p
          invB = ffValue (ffRecip (FF b p))
          factorM = modPow invB m p
          -- Giant steps: search for target * (b^(-m))^i in baby steps
          goGiant !i !gamma
            | i >= m = Nothing
            | otherwise = case Map.lookup gamma babySteps of
                Just j  -> Just (i * m + j)
                Nothing -> goGiant (i + 1) ((gamma * factorM) `mod` p)
      in goGiant 0 target

instance Ord FF where
  compare (FF v1 p1) (FF v2 p2) =
    case compare p1 p2 of
      EQ -> compare v1 v2
      other -> other

instance Show FF where
  show (FF val p) = show val ++ " mod " ++ show p

instance Num FF where
  (FF v1 p1) + (FF v2 p2)
    | p1 == p2  = FF ((v1 + v2) `mod` p1) p1
    | otherwise = error "FF (+): characteristics do not match"

  (FF v1 p1) - (FF v2 p2)
    | p1 == p2  = FF ((v1 - v2 + p1) `mod` p1) p1
    | otherwise = error "FF (-): characteristics do not match"

  (FF v1 p1) * (FF v2 p2)
    | p1 == p2  = FF ((v1 * v2) `mod` p1) p1
    | otherwise = error "FF (*): characteristics do not match"

  negate (FF 0 p) = FF 0 p
  negate (FF v p) = FF (p - v) p

  abs x = x
  signum (FF 0 p) = FF 0 p
  signum (FF _ p) = FF 1 p

  fromInteger n = error $ "FF fromInteger requires characteristic: use (ff " ++ show n ++ " p)"

instance Fractional FF where
  fromRational _ = error "FF fromRational requires characteristic: use (ff n p / ff d p)"
  (/) = ffDiv
  recip = ffRecip
