{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Core.AlgebraicNumber
  ( AlgNum(..)
  , algRoot
  , fromRatAlg
  , algDegree
  , algCoeffs
  , algMinPoly
  , algInterval
  , algRecip
  , refineRoot
  , algApprox
  , evalPolyAt
  ) where

import GHC.Generics (Generic)
import Syzygy.Core.Rational (Rat(..), rat, recipRat)

-- | Exact algebraic number representing an element of the number field Q(alpha)
-- where alpha is the unique real root of the minimal polynomial P(x) in the isolating interval [lo, hi].
-- An element is expressed as c_0 + c_1*alpha + ... + c_{d-1}*alpha^{d-1}.
data AlgNum = AlgNum
  { anCoeffs   :: ![Rat]        -- ^ Ascending coefficients of polynomial in alpha
  , anMinPoly  :: ![Rat]        -- ^ Monic minimal polynomial coefficients [p_0, ..., p_d]
  , anInterval :: !(Rat, Rat)   -- ^ Isolating interval [lo, hi] isolating alpha
  } deriving (Generic)

-- | Trims trailing zero coefficients from a rational polynomial list.
trimRatPoly :: [Rat] -> [Rat]
trimRatPoly [] = []
trimRatPoly coeffs =
  let rev = dropWhile (== 0) (reverse coeffs)
  in reverse rev

-- | Degree of a rational polynomial. Returns -1 for the zero polynomial.
polyDegreeRat :: [Rat] -> Int
polyDegreeRat coeffs =
  case trimRatPoly coeffs of
    [] -> -1
    xs -> length xs - 1

-- | Evaluates a polynomial with rational coefficients at a rational point.
evalPolyAt :: [Rat] -> Rat -> Rat
evalPolyAt [] _ = 0
evalPolyAt coeffs x = foldr (\c acc -> c + x * acc) 0 coeffs

-- | Polynomial addition over Q.
polyAddRat :: [Rat] -> [Rat] -> [Rat]
polyAddRat [] ys = ys
polyAddRat xs [] = xs
polyAddRat (x:xs) (y:ys) = (x + y) : polyAddRat xs ys

-- | Polynomial subtraction over Q.
polySubRat :: [Rat] -> [Rat] -> [Rat]
polySubRat [] ys = map negate ys
polySubRat xs [] = xs
polySubRat (x:xs) (y:ys) = (x - y) : polySubRat xs ys

-- | Polynomial multiplication over Q.
polyMulRat :: [Rat] -> [Rat] -> [Rat]
polyMulRat [] _ = []
polyMulRat _ [] = []
polyMulRat (x:xs) ys =
  let term = map (x *) ys
      rest = 0 : polyMulRat xs ys
  in polyAddRat term rest

-- | Polynomial division with remainder over Q.
-- Given a(x) and b(x) with b /= 0, computes (quotient, remainder)
-- such that a = quotient * b + remainder with deg(remainder) < deg(b).
polyQuotRemRat :: [Rat] -> [Rat] -> ([Rat], [Rat])
polyQuotRemRat a b
  | null (trimRatPoly b) = error "polyQuotRemRat: division by zero polynomial"
  | otherwise = go (trimRatPoly a) []
  where
    degB = polyDegreeRat b
    leadB = case trimRatPoly b of
      [] -> error "polyQuotRemRat: empty divisor"
      xs -> last xs
    invLeadB = recipRat leadB

    go currA accQ
      | polyDegreeRat currA < degB = (trimRatPoly accQ, trimRatPoly currA)
      | otherwise =
          let degCurr = polyDegreeRat currA
              shift = degCurr - degB
              leadA = case trimRatPoly currA of
                [] -> error "polyQuotRemRat: unexpected empty polynomial"
                xs -> last xs
              factor = leadA * invLeadB
              termQ = replicate shift 0 ++ [factor]
              subtractor = polyMulRat termQ b
              nextA = polySubRat currA subtractor
              nextAccQ = polyAddRat accQ termQ
          in go nextA nextAccQ

-- | Extended Euclidean Algorithm for rational polynomials.
-- Computes (g, s, t) such that s*a + t*b = g with g monic.
polyExtendedGCDRat :: [Rat] -> [Rat] -> ([Rat], [Rat], [Rat])
polyExtendedGCDRat a b
  | null (trimRatPoly b) =
      case trimRatPoly a of
        [] -> ([], [], [])
        leadList ->
          let lead = last leadList
              invLead = recipRat lead
          in (map (* invLead) a, [invLead], [0])
  | otherwise =
      let (q, r) = polyQuotRemRat a b
          (g, s1, t1) = polyExtendedGCDRat b r
          s = t1
          t = polySubRat s1 (polyMulRat q t1)
      in (g, s, t)

-- | Reduces a polynomial modulo the minimal polynomial.
reduceModRat :: [Rat] -> [Rat] -> [Rat]
reduceModRat poly minP =
  let (_, remPoly) = polyQuotRemRat poly minP
  in trimRatPoly remPoly

-- | Creates an algebraic number representing the primitive root alpha.
-- Takes the minimal polynomial coefficients [p_0, ..., p_d] and an isolating interval (lo, hi).
algRoot :: [Rat] -> (Rat, Rat) -> AlgNum
algRoot minPoly interval@(lo, hi)
  | polyDegreeRat minPoly <= 0 = error "algRoot: minimal polynomial degree must be >= 1"
  | lo >= hi = error "algRoot: isolating interval must have lo < hi"
  | otherwise =
      let monicP = case trimRatPoly minPoly of
            [] -> error "algRoot: zero minimal polynomial"
            xs -> let lead = last xs in map (/ lead) xs
      in AlgNum [0, 1] monicP interval

-- | Embeds an exact rational number into AlgNum.
fromRatAlg :: Rat -> AlgNum
fromRatAlg r = AlgNum [r] [-r, 1] (r, r + 1)

-- | Degree of the minimal polynomial of the extension.
algDegree :: AlgNum -> Int
algDegree (AlgNum _ p _) = polyDegreeRat p

-- | Extracts the polynomial coefficients in alpha.
algCoeffs :: AlgNum -> [Rat]
algCoeffs (AlgNum c _ _) = c

-- | Extracts the minimal polynomial.
algMinPoly :: AlgNum -> [Rat]
algMinPoly (AlgNum _ p _) = p

-- | Extracts the isolating interval.
algInterval :: AlgNum -> (Rat, Rat)
algInterval (AlgNum _ _ iv) = iv

-- | Multiplicative inverse of an algebraic number.
algRecip :: AlgNum -> AlgNum
algRecip (AlgNum c minP iv)
  | null (trimRatPoly c) = error "algRecip: division by zero"
  | otherwise =
      let (g, s, _) = polyExtendedGCDRat c minP
      in case trimRatPoly g of
           [lead] | lead == 1 -> AlgNum (reduceModRat s minP) minP iv
           _ -> error "algRecip: element not invertible (minimal polynomial reducible)"

-- | Refines the isolating interval of the primitive generator alpha
-- by performing n bisection steps.
refineRoot :: Int -> AlgNum -> AlgNum
refineRoot steps (AlgNum c minP (lo, hi)) =
  let newIv = go steps lo hi
  in AlgNum c minP newIv
  where
    go 0 l h = (l, h)
    go !n l h =
      let mid = (l + h) / 2
          valMid = evalPolyAt minP mid
          valLo  = evalPolyAt minP l
      in if valMid == 0
         then (mid, mid + rat 1 1000000000)
         else if (valLo > 0 && valMid < 0) || (valLo < 0 && valMid > 0)
              then go (n - 1) l mid
              else go (n - 1) mid h

-- | Numerical floating point approximation of an algebraic number.
algApprox :: Int -> AlgNum -> Double
algApprox refinementSteps an@(AlgNum c _ _) =
  let AlgNum _ _ (refLo, refHi) = refineRoot refinementSteps an
      mid = (refLo + refHi) / 2
      Rat n d = evalPolyAt c mid
  in fromInteger n / fromInteger d

instance Eq AlgNum where
  (AlgNum c1 p1 _) == (AlgNum c2 p2 _)
    | p1 == p2  = trimRatPoly c1 == trimRatPoly c2
    | otherwise = False

instance Ord AlgNum where
  compare an1@(AlgNum c1 p1 iv1) an2@(AlgNum c2 p2 iv2)
    | an1 == an2 = EQ
    | p1 == p2 =
        -- Compare by evaluating at a sufficiently refined point in isolating interval
        let refined1 = refineRoot 30 an1
            refined2 = refineRoot 30 an2
            (_, hi1) = anInterval refined1
            (_, hi2) = anInterval refined2
            v1 = evalPolyAt c1 hi1
            v2 = evalPolyAt c2 hi2
        in compare v1 v2
    | otherwise =
        let val1 = algApprox 40 (AlgNum c1 p1 iv1)
            val2 = algApprox 40 (AlgNum c2 p2 iv2)
        in compare val1 val2

instance Show AlgNum where
  show (AlgNum [] _ _) = "0"
  show an@(AlgNum c _ _) =
    let terms = [ show (Rat num den) ++ (if d == 0 then "" else if d == 1 then "*alpha" else "*alpha^" ++ show d)
                | (d, Rat num den) <- zip [0 :: Int ..] c
                , num /= 0
                ]
        polyStr = if null terms then "0" else unwords (reverse terms)
        approxVal = algApprox 20 an
    in polyStr ++ " (~ " ++ show approxVal ++ ")"

instance Num AlgNum where
  (AlgNum c1 p1 iv1) + (AlgNum c2 p2 _)
    | p1 == p2  = AlgNum (reduceModRat (polyAddRat c1 c2) p1) p1 iv1
    | otherwise = error "AlgNum (+): generators belong to distinct algebraic extensions"

  (AlgNum c1 p1 iv1) - (AlgNum c2 p2 _)
    | p1 == p2  = AlgNum (reduceModRat (polySubRat c1 c2) p1) p1 iv1
    | otherwise = error "AlgNum (-): generators belong to distinct algebraic extensions"

  (AlgNum c1 p1 iv1) * (AlgNum c2 p2 _)
    | p1 == p2  = AlgNum (reduceModRat (polyMulRat c1 c2) p1) p1 iv1
    | otherwise = error "AlgNum (*): generators belong to distinct algebraic extensions"

  negate (AlgNum c p iv) = AlgNum (map negate c) p iv
  abs x = if x >= 0 then x else negate x
  signum x
    | x > 0     = fromRatAlg 1
    | x < 0     = fromRatAlg (-1)
    | otherwise = fromRatAlg 0

  fromInteger n = fromRatAlg (Rat n 1)

instance Fractional AlgNum where
  fromRational r = fromRatAlg (fromRational r)
  a / b = a * algRecip b
  recip = algRecip
