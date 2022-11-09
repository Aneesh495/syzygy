{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.HilbertSeries
  ( HilbertSeries(..)
  , hilbertSeriesFromResolution
  , hilbertSeriesFromMonomialIdeal
  , hilbertFunction
  , hilbertNumeratorPolynomial
  , krullDimension
  , varietyDegree
  , formatHilbertSeries
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseDegree
  , sparseLCM
  )
import Syzygy.Resolution.FreeResolution (ChainComplex(..))

-- | Rational representation of the Hilbert-Poincare series of a graded module M:
-- HS_M(t) = P(t) / (1 - t)^d
-- where d is the Krull dimension of M, and P(1) = deg(M) is the degree (multiplicity).
data HilbertSeries = HilbertSeries
  { hsNumeratorCoeffs :: ![Integer]  -- ^ Coefficients [a_0, a_1, ..., a_m] of P(t) = sum a_k * t^k
  , hsDenominatorDim  :: !Int        -- ^ Power d of (1 - t)^d in the denominator (Krull dimension)
  , hsTotalVars       :: !Int        -- ^ Number of variables n in the polynomial ring
  } deriving (Eq, Show, Generic)

-- | Computes the Hilbert series directly from a graded chain complex free resolution:
-- HS_M(t) = (sum_{i=0}^p (-1)^i sum_j t^{d_{i,j}}) / (1 - t)^n.
-- Simplifies the rational function by canceling common factors of (1 - t).
hilbertSeriesFromResolution :: ChainComplex c -> Int -> HilbertSeries
hilbertSeriesFromResolution cc numVars =
  let shifts = ccShifts cc
      -- Raw K-polynomial coefficients: sum_{i} (-1)^i sum_j t^{d_{i,j}}
      rawTerms = [ (deg, if even colIdx then 1 else -1)
                 | (colIdx, degList) <- zip [0 :: Int ..] shifts
                 , deg <- degList
                 ]
      rawCoeffMap = foldl (\acc (deg, sign) -> Map.insertWith (+) deg sign acc) Map.empty rawTerms
      maxDeg = if Map.null rawCoeffMap then 0 else fst (Map.findMax rawCoeffMap)
      rawCoeffs = [ Map.findWithDefault 0 d rawCoeffMap | d <- [0 .. maxDeg] ]
  in cancelOneMinusT rawCoeffs numVars numVars

-- | Repeatedly divides numerator polynomial by (1 - t) if P(1) == sum(coeffs) == 0 and dim > 0.
cancelOneMinusT :: [Integer] -> Int -> Int -> HilbertSeries
cancelOneMinusT coeffs dim n
  | dim > 0 && sum coeffs == 0 && not (null coeffs) =
      let qCoeffs = scanl1 (+) coeffs
          qTrimmed = reverse (dropWhile (== 0) (reverse (take (length coeffs - 1) qCoeffs)))
          finalQ = if null qTrimmed then [0] else qTrimmed
      in cancelOneMinusT finalQ (dim - 1) n
  | otherwise =
      let trimmed = reverse (dropWhile (== 0) (reverse coeffs))
          finalCoeffs = if null trimmed then [0] else trimmed
      in HilbertSeries finalCoeffs dim n

-- | Computes the Hilbert series of a monomial ideal <m_1, ..., m_s>
-- using the Taylor resolution / inclusion-exclusion on monomial LCMs:
-- K(t) = 1 - sum t^{deg m_i} + sum_{i < j} t^{deg LCM(m_i, m_j)} - ...
hilbertSeriesFromMonomialIdeal :: Int -> [SparseMonomial] -> HilbertSeries
hilbertSeriesFromMonomialIdeal numVars monos =
  let distinctMonos = filter (/= SparseMonomial Map.empty) monos
      n = length distinctMonos
      -- Generate all non-empty subsets up to size min(n, 10) for inclusion-exclusion
      subsets = generateSubsets (min n 10) distinctMonos
      terms = [ (sparseDegree (foldl1 sparseLCM s), if odd (length s) then -1 else 1)
              | s <- subsets
              ]
      baseMap = Map.singleton (0 :: Int) (1 :: Integer) -- the constant 1
      rawCoeffMap = foldl (\acc (deg, sign) -> Map.insertWith (+) deg sign acc) baseMap terms
      maxDeg = if Map.null rawCoeffMap then 0 else fst (Map.findMax rawCoeffMap)
      rawCoeffs = [ Map.findWithDefault 0 d rawCoeffMap | d <- [0 .. maxDeg] ]
  in cancelOneMinusT rawCoeffs numVars numVars
  where
    generateSubsets maxK xs =
      concat [ combinations k xs | k <- [1 .. maxK] ]

    combinations 0 _ = [[]]
    combinations _ [] = []
    combinations k (y:ys) =
      map (y :) (combinations (k - 1) ys) ++ combinations k ys

-- | Evaluates the Hilbert function HF_M(k) = dim_k (M_k) at degree k.
-- Derived by expanding P(t) / (1 - t)^d:
-- 1 / (1 - t)^d = sum_{j >= 0} binom(j + d - 1, d - 1) * t^j.
hilbertFunction :: HilbertSeries -> Int -> Integer
hilbertFunction (HilbertSeries numCoeffs d _) k
  | k < 0 = 0
  | d == 0 =
      if k < length numCoeffs then numCoeffs !! k else 0
  | otherwise =
      sum [ (numCoeffs !! i) * binomial (k - i + d - 1) (d - 1)
          | i <- [0 .. min k (length numCoeffs - 1)]
          ]
  where
    binomial n r
      | r < 0 || r > n = 0
      | r == 0 || r == n = 1
      | otherwise =
          let r' = min r (n - r)
              num = product [ fromIntegral (n - j) | j <- [0 .. r' - 1] ]
              den = product [ fromIntegral (j + 1) | j <- [0 .. r' - 1] ]
          in num `div` den

-- | Returns the numerator polynomial coefficients [a_0, ..., a_m] of P(t).
hilbertNumeratorPolynomial :: HilbertSeries -> [Integer]
hilbertNumeratorPolynomial = hsNumeratorCoeffs

-- | Returns the Krull dimension of the module (dimension of the algebraic variety).
krullDimension :: HilbertSeries -> Int
krullDimension = hsDenominatorDim

-- | Returns the algebraic degree (multiplicity) of the projective variety: deg(M) = P(1).
varietyDegree :: HilbertSeries -> Int
varietyDegree hs = fromIntegral (sum (hsNumeratorCoeffs hs))

-- | Formats the Hilbert series as a human-readable rational function string:
-- e.g., "(1 + 2*t + t^2) / (1 - t)^3"
formatHilbertSeries :: HilbertSeries -> String
formatHilbertSeries (HilbertSeries coeffs d _) =
  let numStr = formatPoly coeffs
      denStr = case d of
                 0 -> ""
                 1 -> " / (1 - t)"
                 _ -> " / (1 - t)^" ++ show d
  in if d == 0 then numStr else "(" ++ numStr ++ ")" ++ denStr
  where
    formatPoly [] = "0"
    formatPoly [c] = show c
    formatPoly cs =
      let terms = [ formatTerm c p | (p, c) <- zip [0 :: Int ..] cs, c /= 0 ]
      in if null terms then "0" else unwords (joinWithPlus terms)

    formatTerm c 0 = show c
    formatTerm 1 1 = "t"
    formatTerm (-1) 1 = "-t"
    formatTerm c 1 = show c ++ "*t"
    formatTerm 1 p = "t^" ++ show p
    formatTerm (-1) p = "-t^" ++ show p
    formatTerm c p = show c ++ "*t^" ++ show p

    joinWithPlus [] = []
    joinWithPlus (x:xs) = x : concatMap (\t -> ["+", t]) xs
