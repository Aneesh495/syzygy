{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Core.GaloisField
  ( GF(..)
  , gf
  , gfDegree
  , gfCharacteristic
  , gfCoeffs
  , gfModulus
  , gfRecip
  , gfPow
  , frobenius
  , gfTrace
  , gfNorm
  ) where

import GHC.Generics (Generic)
import Syzygy.Core.FiniteField (FF(..), ff, ffCharacteristic, ffValue, ffRecip)

-- | Element of the Galois Field extension F_{p^k} represented as a polynomial
-- modulo an irreducible polynomial f(x) over F_p.
-- The first list represents the polynomial coefficients [c_0, c_1, ..., c_{k-1}]
-- in ascending degree order, and the second list is the monic modulus [m_0, ..., m_k].
data GF = GF ![FF] ![FF]
  deriving (Eq, Generic)

-- | Trims trailing zero coefficients from a polynomial list.
trimPoly :: [FF] -> [FF]
trimPoly [] = []
trimPoly coeffs =
  let rev = dropWhile (\(FF v _) -> v == 0) (reverse coeffs)
  in reverse rev

-- | Degree of a polynomial represented by coefficients in ascending degree.
polyDegree :: [FF] -> Int
polyDegree coeffs =
  case trimPoly coeffs of
    [] -> -1
    xs -> length xs - 1

-- | Polynomial addition over F_p.
polyAdd :: [FF] -> [FF] -> [FF]
polyAdd [] ys = ys
polyAdd xs [] = xs
polyAdd (x:xs) (y:ys) = (x + y) : polyAdd xs ys

-- | Polynomial subtraction over F_p.
polySub :: [FF] -> [FF] -> [FF]
polySub [] ys = map negate ys
polySub xs [] = xs
polySub (x:xs) (y:ys) = (x - y) : polySub xs ys

-- | Polynomial multiplication over F_p.
polyMul :: [FF] -> [FF] -> [FF]
polyMul [] _ = []
polyMul _ [] = []
polyMul (x:xs) ys =
  let term = map (x *) ys
      rest = FF 0 (ffCharacteristic x) : polyMul xs ys
  in polyAdd term rest

-- | Polynomial division with remainder over F_p.
-- Given a(x) and b(x) with b /= 0, returns (quotient, remainder) such that
-- a(x) = quotient(x) * b(x) + remainder(x) with deg(remainder) < deg(b).
polyQuotRem :: [FF] -> [FF] -> ([FF], [FF])
polyQuotRem a b
  | null (trimPoly b) = error "polyQuotRem: division by zero polynomial"
  | otherwise = go (trimPoly a) []
  where
    degB = polyDegree b
    leadB = case trimPoly b of
      [] -> error "polyQuotRem: empty divisor"
      xs -> last xs
    invLeadB = ffRecip leadB

    go currA accQ
      | polyDegree currA < degB = (trimPoly accQ, trimPoly currA)
      | otherwise =
          let degCurr = polyDegree currA
              shift = degCurr - degB
              leadA = case trimPoly currA of
                [] -> error "polyQuotRem: unexpected empty polynomial"
                xs -> last xs
              factor = leadA * invLeadB
              termQ = replicate shift (FF 0 (ffCharacteristic leadA)) ++ [factor]
              subtractor = polyMul termQ b
              nextA = polySub currA subtractor
              nextAccQ = polyAdd accQ termQ
          in go nextA nextAccQ

-- | Extended Euclidean Algorithm for polynomials over F_p.
-- Computes (g, s, t) such that s*a + t*b = g = gcd(a, b) with g monic.
polyExtendedGCD :: [FF] -> [FF] -> ([FF], [FF], [FF])
polyExtendedGCD a b
  | null (trimPoly b) =
      case trimPoly a of
        [] -> ([], [], [])
        leadList ->
          let lead = last leadList
              invLead = ffRecip lead
          in (map (* invLead) a, [invLead], [FF 0 (ffCharacteristic lead)])
  | otherwise =
      let (q, r) = polyQuotRem a b
          (g, s1, t1) = polyExtendedGCD b r
          s = t1
          t = polySub s1 (polyMul q t1)
      in (g, s, t)

-- | Reduces a polynomial modulo the irreducible modulus.
reduceMod :: [FF] -> [FF] -> [FF]
reduceMod poly modulus =
  let (_, remPoly) = polyQuotRem poly modulus
  in trimPoly remPoly

-- | Smart constructor for Galois Field elements.
-- Arguments: integer coefficients, integer modulus coefficients, prime characteristic p.
gf :: [Integer] -> [Integer] -> Integer -> GF
gf coeffInts modInts p
  | p <= 1 = error "gf: characteristic must be a prime > 1"
  | polyDegree modFF <= 0 = error "gf: modulus must have degree >= 1"
  | otherwise =
      let redCoeffs = reduceMod coeffFF modFF
      in GF redCoeffs modFF
  where
    coeffFF = map (`ff` p) coeffInts
    modFF   = map (`ff` p) modInts

-- | Returns the extension degree k of F_{p^k}.
gfDegree :: GF -> Int
gfDegree (GF _ modFF) = polyDegree modFF

-- | Returns the prime characteristic p of the base field.
gfCharacteristic :: GF -> Integer
gfCharacteristic (GF _ modFF) =
  case modFF of
    (m:_) -> ffCharacteristic m
    []    -> error "gfCharacteristic: empty modulus"

-- | Returns the polynomial coefficients in ascending degree order.
gfCoeffs :: GF -> [FF]
gfCoeffs (GF c _) = c

-- | Returns the modulus polynomial coefficients.
gfModulus :: GF -> [FF]
gfModulus (GF _ m) = m

-- | Multiplicative inverse in F_{p^k} using the extended Euclidean algorithm.
gfRecip :: GF -> GF
gfRecip (GF c modFF)
  | null (trimPoly c) = error "gfRecip: division by zero in Galois field"
  | otherwise =
      let (g, s, _) = polyExtendedGCD c modFF
      in case trimPoly g of
           [lead] | lead == FF 1 (ffCharacteristic lead) ->
             GF (reduceMod s modFF) modFF
           _ -> error "gfRecip: modulus is reducible, element not invertible"

-- | Exponentiation in F_{p^k} using binary exponentiation.
gfPow :: GF -> Integer -> GF
gfPow (GF _ m) 0 = GF [FF 1 (case m of (x:_) -> ffCharacteristic x; [] -> 2)] m
gfPow a expVal
  | expVal < 0 = gfPow (gfRecip a) (-expVal)
  | otherwise  = go a expVal (gf [1] (map ffValue (gfModulus a)) (gfCharacteristic a))
  where
    go !_ 0 !acc = acc
    go !base e !acc =
      let !acc' = if odd e then acc * base else acc
          !base' = base * base
          !e' = e `div` 2
      in go base' e' acc'

-- | Frobenius endomorphism: sigma(a) = a^p.
frobenius :: GF -> GF
frobenius a = gfPow a (gfCharacteristic a)

-- | Field trace Tr_{F_{p^k} / F_p}(a) = sum_{i=0}^{k-1} a^{p^i}.
-- The trace always evaluates to an element of the base field F_p.
gfTrace :: GF -> FF
gfTrace a =
  let k = gfDegree a
      powers = take k (iterate frobenius a)
      GF traceCoeffs _ = sum powers
  in case traceCoeffs of
       []    -> FF 0 (gfCharacteristic a)
       (c:_) -> c

-- | Field norm N_{F_{p^k} / F_p}(a) = prod_{i=0}^{k-1} a^{p^i}.
-- The norm always evaluates to an element of the base field F_p.
gfNorm :: GF -> FF
gfNorm a =
  let k = gfDegree a
      powers = take k (iterate frobenius a)
      GF normCoeffs _ = product powers
  in case normCoeffs of
       []    -> FF 0 (gfCharacteristic a)
       (c:_) -> c

instance Ord GF where
  compare (GF c1 m1) (GF c2 m2) =
    case compare m1 m2 of
      EQ -> compare c1 c2
      other -> other

instance Show GF where
  show g@(GF [] _) = "0 in GF(" ++ show (gfCharacteristic g) ++ "^" ++ show (gfDegree g) ++ ")"
  show g@(GF c _) =
    let terms = [ show (ffValue coeff) ++ (if d == 0 then "" else if d == 1 then "*alpha" else "*alpha^" ++ show d)
                | (d, coeff) <- zip [0 :: Int ..] c
                , ffValue coeff /= 0
                ]
        polyStr = if null terms then "0" else unwords (reverse terms)
        p = gfCharacteristic g
        k = gfDegree g
    in polyStr ++ " in GF(" ++ show p ++ "^" ++ show k ++ ")"

instance Num GF where
  (GF c1 m1) + (GF c2 m2)
    | m1 == m2  = GF (reduceMod (polyAdd c1 c2) m1) m1
    | otherwise = error "GF (+): modulus polynomials do not match"

  (GF c1 m1) - (GF c2 m2)
    | m1 == m2  = GF (reduceMod (polySub c1 c2) m1) m1
    | otherwise = error "GF (-): modulus polynomials do not match"

  (GF c1 m1) * (GF c2 m2)
    | m1 == m2  = GF (reduceMod (polyMul c1 c2) m1) m1
    | otherwise = error "GF (*): modulus polynomials do not match"

  negate (GF c m) = GF (map negate c) m
  abs x = x
  signum x = x

  fromInteger n = error $ "GF fromInteger requires modulus: use (gf [" ++ show n ++ "] modPoly p)"

instance Fractional GF where
  fromRational _ = error "GF fromRational requires modulus context"
  a / b = a * gfRecip b
  recip = gfRecip
