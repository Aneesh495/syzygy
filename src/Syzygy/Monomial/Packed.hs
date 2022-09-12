{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Monomial.Packed
  ( PackedMonomial(..)
  , packMonomial
  , unpackMonomial
  , packedDegree
  , packedMul
  , packedDivides
  , packedDiv
  , packedLCM
  , packedGCD
  , packedCoprime
  , packedZero
  ) where

import Data.Bits (shiftL, shiftR, (.&.), (.|.))
import Data.Word (Word64)
import GHC.Generics (Generic)

-- | Bit-packed 64-bit exponent vector supporting up to 7 variables with degrees up to 255.
-- Layout:
-- Bits [63:56]: Total degree (8 bits)
-- Bits [55:48]: Exponent x_0 (8 bits)
-- Bits [47:40]: Exponent x_1 (8 bits)
-- Bits [39:32]: Exponent x_2 (8 bits)
-- Bits [31:24]: Exponent x_3 (8 bits)
-- Bits [23:16]: Exponent x_4 (8 bits)
-- Bits [15:8]:  Exponent x_5 (8 bits)
-- Bits [7:0]:   Exponent x_6 (8 bits)
newtype PackedMonomial = PackedMonomial { unPackedMonomial :: Word64 }
  deriving (Eq, Ord, Generic)

-- | Zero monomial 1 = x_0^0 ... x_6^0.
packedZero :: PackedMonomial
packedZero = PackedMonomial 0

-- | Packs an exponent list of length <= 7 into a 64-bit word.
packMonomial :: [Int] -> PackedMonomial
packMonomial exps =
  let padded = take 7 (exps ++ repeat 0)
      deg = sum padded
      degWord = fromIntegral (deg .&. 0xFF) :: Word64
      packVar (idx, val) =
        let w = fromIntegral (val .&. 0xFF) :: Word64
            shiftBits = (6 - idx) * 8
        in w `shiftL` shiftBits
      varWords = foldl (.|.) 0 (map packVar (zip [0 :: Int .. 6] padded))
      totalWord = (degWord `shiftL` 56) .|. varWords
  in PackedMonomial totalWord

-- | Unpacks a PackedMonomial back to an exponent list of length n (n <= 7).
unpackMonomial :: Int -> PackedMonomial -> [Int]
unpackMonomial n (PackedMonomial w) =
  let numVars = min 7 (max 0 n)
      extractVar idx =
        let shiftBits = (6 - idx) * 8
            val = (w `shiftR` shiftBits) .&. 0xFF
        in fromIntegral val
  in [ extractVar i | i <- [0 .. numVars - 1] ]

-- | Extracts the total degree directly from the top byte in O(1).
packedDegree :: PackedMonomial -> Int
packedDegree (PackedMonomial w) =
  fromIntegral ((w `shiftR` 56) .&. 0xFF)

-- | Multiplies two packed monomials: adds total degree and component exponents.
packedMul :: PackedMonomial -> PackedMonomial -> PackedMonomial
packedMul (PackedMonomial w1) (PackedMonomial w2) =
  -- Since degrees and exponents fit in 8-bit fields without carry into adjacent fields,
  -- packed addition can be performed directly on unpacked components or masked additions.
  let deg1 = (w1 `shiftR` 56) .&. 0xFF
      deg2 = (w2 `shiftR` 56) .&. 0xFF
      degSum = min 255 (deg1 + deg2) `shiftL` 56
      addByte shiftBits =
        let b1 = (w1 `shiftR` shiftBits) .&. 0xFF
            b2 = (w2 `shiftR` shiftBits) .&. 0xFF
        in min 255 (b1 + b2) `shiftL` shiftBits
      varsSum = foldl (.|.) 0 [ addByte (s * 8) | s <- [0 .. 6] ]
  in PackedMonomial (degSum .|. varsSum)

-- | Tests if monomial u divides monomial v.
-- u divides v iff for all variables x_i, exp(u, x_i) <= exp(v, x_i).
packedDivides :: PackedMonomial -> PackedMonomial -> Bool
packedDivides (PackedMonomial w1) (PackedMonomial w2) =
  let deg1 = (w1 `shiftR` 56) .&. 0xFF
      deg2 = (w2 `shiftR` 56) .&. 0xFF
  in if deg1 > deg2
     then False
     else all checkByte [0 .. 6]
  where
    checkByte s =
      let shiftBits = s * 8
          b1 = (w1 `shiftR` shiftBits) .&. 0xFF
          b2 = (w2 `shiftR` shiftBits) .&. 0xFF
      in b1 <= b2

-- | Exact quotient of v / u when u divides v.
packedDiv :: PackedMonomial -> PackedMonomial -> PackedMonomial
packedDiv (PackedMonomial wU) (PackedMonomial wV) =
  let degU = (wU `shiftR` 56) .&. 0xFF
      degV = (wV `shiftR` 56) .&. 0xFF
      degDiff = (degV - degU) `shiftL` 56
      subByte shiftBits =
        let bU = (wU `shiftR` shiftBits) .&. 0xFF
            bV = (wV `shiftR` shiftBits) .&. 0xFF
        in (bV - bU) `shiftL` shiftBits
      varsDiff = foldl (.|.) 0 [ subByte (s * 8) | s <- [0 .. 6] ]
  in PackedMonomial (degDiff .|. varsDiff)

-- | Computes the least common multiple (LCM) of two packed monomials:
-- component-wise maximum of exponents.
packedLCM :: PackedMonomial -> PackedMonomial -> PackedMonomial
packedLCM (PackedMonomial w1) (PackedMonomial w2) =
  let maxByte shiftBits =
        let b1 = (w1 `shiftR` shiftBits) .&. 0xFF
            b2 = (w2 `shiftR` shiftBits) .&. 0xFF
        in max b1 b2
      bytes = [ maxByte (s * 8) | s <- [0 .. 6] ]
      newDeg = sum bytes `shiftL` 56
      newVars = foldl (.|.) 0 [ b `shiftL` (s * 8) | (s, b) <- zip [0 .. 6] bytes ]
  in PackedMonomial (newDeg .|. newVars)

-- | Computes the greatest common divisor (GCD) of two packed monomials:
-- component-wise minimum of exponents.
packedGCD :: PackedMonomial -> PackedMonomial -> PackedMonomial
packedGCD (PackedMonomial w1) (PackedMonomial w2) =
  let minByte shiftBits =
        let b1 = (w1 `shiftR` shiftBits) .&. 0xFF
            b2 = (w2 `shiftR` shiftBits) .&. 0xFF
        in min b1 b2
      bytes = [ minByte (s * 8) | s <- [0 .. 6] ]
      newDeg = sum bytes `shiftL` 56
      newVars = foldl (.|.) 0 [ b `shiftL` (s * 8) | (s, b) <- zip [0 .. 6] bytes ]
  in PackedMonomial (newDeg .|. newVars)

-- | Tests if two monomials are coprime: gcd(u, v) == 1.
-- Gebauer-Moller Criterion 1 states that S(f, g) reduces to 0 without division if LM(f) and LM(g) are coprime!
packedCoprime :: PackedMonomial -> PackedMonomial -> Bool
packedCoprime u v =
  packedGCD u v == packedZero

instance Show PackedMonomial where
  show pm =
    let exps = unpackMonomial 7 pm
        deg = packedDegree pm
    in "PackedMonomial " ++ show exps ++ " (deg " ++ show deg ++ ")"
