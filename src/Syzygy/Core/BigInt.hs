{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Core.BigInt
  ( extendedGCD
  , modInverse
  , modPow
  , chineseRemainder
  , isMillerRabinPrime
  , isPrime
  , nextPrime
  , smallPrimes
  , eulerPhi
  , jacobiSymbol
  , integerSqrt
  , factorize
  , factorial
  , binomial
  ) where

import Data.Bits (shiftR, testBit)

-- | Extended Euclidean Algorithm.
-- Given integers a and b, computes (g, x, y) such that:
-- a * x + b * y = g = gcd(a, b) with g >= 0.
extendedGCD :: Integer -> Integer -> (Integer, Integer, Integer)
extendedGCD a b
  | b == 0 =
      if a < 0 then (-a, -1, 0) else (a, 1, 0)
  | otherwise =
      let (q, r) = a `quotRem` b
          (g, x1, y1) = extendedGCD b r
          x = y1
          y = x1 - q * y1
      in if g < 0 then (-g, -x, -y) else (g, x, y)

-- | Computes the modular inverse of a modulo m.
-- Returns Just inv if gcd(a, m) == 1, or Nothing otherwise.
modInverse :: Integer -> Integer -> Maybe Integer
modInverse _ 0 = Nothing
modInverse _ 1 = Just 0
modInverse a m
  | m < 0     = modInverse a (-m)
  | otherwise =
      let (g, x, _) = extendedGCD a m
      in if g == 1
         then Just ((x `mod` m + m) `mod` m)
         else Nothing

-- | Modular exponentiation: computes (base ^ exponent) mod modulus
-- using the binary square and multiply algorithm.
modPow :: Integer -> Integer -> Integer -> Integer
modPow _ _ 0 = error "modPow: division by zero modulus"
modPow _ _ 1 = 0
modPow base expVal modulus
  | modulus < 0 = modPow base expVal (-modulus)
  | expVal < 0  = case modInverse base modulus of
      Just inv -> modPow inv (-expVal) modulus
      Nothing  -> error "modPow: base not invertible for negative exponent"
  | otherwise   = go (base `mod` modulus) expVal 1
  where
    go !_ 0 !acc = acc
    go !b !e !acc =
      let !acc' = if testBit e 0 then (acc * b) `mod` modulus else acc
          !b'   = (b * b) `mod` modulus
          !e'   = e `shiftR` 1
      in go b' e' acc'

-- | Chinese Remainder Theorem.
-- Given a list of (remainder, modulus) pairs where moduli are pairwise coprime,
-- finds (x, M) such that x is congruent to each remainder modulo its modulus,
-- where M is the product of all moduli and 0 <= x < M.
chineseRemainder :: [(Integer, Integer)] -> Maybe (Integer, Integer)
chineseRemainder [] = Just (0, 1)
chineseRemainder pairs
  | any (\(_, m) -> m <= 0) pairs = Nothing
  | otherwise = foldl combine (Just (0, 1)) pairs
  where
    combine Nothing _ = Nothing
    combine (Just (x1, m1)) (x2, m2) =
      let (g, s, _) = extendedGCD m1 m2
      in if (x2 - x1) `mod` g /= 0
         then Nothing
         else
           let m = (m1 * m2) `div` g
               step = ((x2 - x1) `div` g) * s `mod` (m2 `div` g)
               x = (x1 + step * m1) `mod` m
           in Just ((x + m) `mod` m, m)

-- | Miller-Rabin probabilistic primality test.
-- The integer k specifies the number of randomized bases to test.
isMillerRabinPrime :: Int -> Integer -> Bool
isMillerRabinPrime k n
  | n < 2     = False
  | n == 2    = True
  | n == 3    = True
  | even n    = False
  | otherwise =
      let (d, s) = factorTwo (n - 1)
          bases  = deterministicBases n
      in all (checkBase n d s) (take k bases)
  where
    factorTwo m
      | even m    = let (d, s) = factorTwo (m `shiftR` 1) in (d, s + 1)
      | otherwise = (m, 0 :: Int)

    checkBase _target d s a
      | a `mod` _target == 0 = True
      | otherwise =
          let x0 = modPow a d _target
          in if x0 == 1 || x0 == _target - 1
             then True
             else loop (s - 1) x0
      where
        loop 0 _  = False
        loop !r !x =
          let x' = (x * x) `mod` _target
          in if x' == _target - 1
             then True
             else if x' == 1
                  then False
                  else loop (r - 1) x'

    deterministicBases val
      | val < 2047                 = [2]
      | val < 1373653              = [2, 3]
      | val < 9080191              = [31, 73]
      | val < 25326001             = [2, 3, 5]
      | val < 3215031751           = [2, 3, 5, 7]
      | val < 4759123141           = [2, 7, 61]
      | val < 1122004669633        = [2, 13, 23, 1662803]
      | val < 2152302898747        = [2, 3, 5, 7, 11]
      | val < 3474749660383        = [2, 3, 5, 7, 11, 13]
      | val < 341550071728321      = [2, 3, 5, 7, 11, 13, 17]
      | otherwise                  = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

-- | Deterministic primality test using verified Miller-Rabin bases.
isPrime :: Integer -> Bool
isPrime = isMillerRabinPrime 12

-- | Computes the smallest prime strictly greater than n.
nextPrime :: Integer -> Integer
nextPrime n
  | n < 2     = 2
  | otherwise =
      let candidate = if even (n + 1) then n + 2 else n + 1
      in findPrime candidate
  where
    findPrime c
      | isPrime c = c
      | otherwise = findPrime (c + 2)

-- | Sieve of Eratosthenes producing the infinite list of primes.
smallPrimes :: [Integer]
smallPrimes = 2 : 3 : filter isPrimeCandidates [5, 7..]
  where
    isPrimeCandidates x = all (\p -> x `rem` p /= 0) (takeWhile (\p -> p * p <= x) smallPrimes)

-- | Euler's totient function phi(n).
-- Counts the number of positive integers up to n that are relatively prime to n.
eulerPhi :: Integer -> Integer
eulerPhi n
  | n <= 0    = 0
  | n == 1    = 1
  | otherwise =
      let factors = factorize n
      in foldl (\acc (p, _) -> acc * (p - 1) `div` p) n factors

-- | Computes the Jacobi symbol (a / n) for odd n > 0.
-- Returns -1, 0, or 1.
jacobiSymbol :: Integer -> Integer -> Int
jacobiSymbol a n
  | n <= 0 || even n = error "jacobiSymbol: n must be positive and odd"
  | otherwise = go (a `mod` n) n 1
  where
    go 0 _ _ = 0
    go 1 _ !acc = acc
    go a' n' !acc
      | even a' =
          let (a'', (k :: Int)) = removeTwos a' (0 :: Int)
              acc' = if even k || n' `mod` 8 == 1 || n' `mod` 8 == 7
                     then acc
                     else -acc
          in go a'' n' acc'
      | otherwise =
          let acc' = if a' `mod` 4 == 3 && n' `mod` 4 == 3 then -acc else acc
          in go (n' `mod` a') a' acc'

    removeTwos x !(k :: Int)
      | even x    = removeTwos (x `shiftR` 1) (k + 1)
      | otherwise = (x, k)

-- | Computes the integer square root floor(sqrt(n)) using Newton's method.
integerSqrt :: Integer -> Integer
integerSqrt n
  | n < 0     = error "integerSqrt: negative argument"
  | n == 0    = 0
  | n < 4     = 1
  | otherwise =
      let initGuess = 1 `shiftR` ((integerLog2 n + 1) `div` 2)
      in newton initGuess
  where
    newton x =
      let nextX = (x + n `div` x) `div` 2
      in if nextX >= x then x else newton nextX

    integerLog2 :: Integer -> Int
    integerLog2 0 = 0
    integerLog2 m = goLog m (0 :: Int)
      where
        goLog 1 !acc = acc
        goLog v !acc = goLog (v `shiftR` 1) (acc + 1)

-- | Prime factorization of an integer n.
-- Returns a list of (prime, exponent) pairs sorted in ascending order of primes.
factorize :: Integer -> [(Integer, Int)]
factorize n
  | n <= 1    = []
  | otherwise = trialDivide n smallPrimes
  where
    trialDivide 1 _ = []
    trialDivide m (p:ps)
      | p * p > m = [(m, 1)]
      | m `rem` p == 0 =
          let (m', count) = extractPower m p 0
          in (p, count) : trialDivide m' ps
      | otherwise = trialDivide m ps
    trialDivide m [] = [(m, 1)]

    extractPower m p !acc
      | m `rem` p == 0 = extractPower (m `div` p) p (acc + 1)
      | otherwise      = (m, acc)

-- | Factorial of a non-negative integer n!
factorial :: Integer -> Integer
factorial n
  | n < 0     = error "factorial: negative integer"
  | n == 0    = 1
  | otherwise = product [1..n]

-- | Binomial coefficient (n choose k).
binomial :: Integer -> Integer -> Integer
binomial n k
  | k < 0 || k > n = 0
  | k == 0 || k == n = 1
  | k > n `div` 2  = binomial n (n - k)
  | otherwise      =
      let num = product [(n - k + 1) .. n]
          den = factorial k
      in num `div` den
