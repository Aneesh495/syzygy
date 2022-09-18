{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Univariate
  ( UnivariatePoly(..)
  , uPoly
  , uDegree
  , uLead
  , uIsZero
  , uAdd
  , uSub
  , uMul
  , uQuotRem
  , uGCD
  , uExtendedGCD
  , uDerivative
  , uSubresultantPRS
  , uResultant
  , uDiscriminant
  , uSquareFree
  , uCauchyBound
  ) where

import Syzygy.Core.Rational (Rat)

-- | Dense univariate polynomial c_0 + c_1*x + ... + c_d*x^d
-- represented as a coefficient list in ASCENDING degree order.
newtype UnivariatePoly c = UnivariatePoly { unUPoly :: [c] }
  deriving (Eq)

-- | Trims trailing zero coefficients.
trimU :: (Eq c, Num c) => [c] -> [c]
trimU [] = []
trimU cs =
  let rev = dropWhile (== 0) (reverse cs)
  in reverse rev

-- | Smart constructor for univariate polynomial.
uPoly :: (Eq c, Num c) => [c] -> UnivariatePoly c
uPoly cs = UnivariatePoly (trimU cs)

-- | Degree of univariate polynomial (-1 for zero polynomial).
uDegree :: (Eq c, Num c) => UnivariatePoly c -> Int
uDegree (UnivariatePoly cs) =
  case trimU cs of
    [] -> -1
    xs -> length xs - 1

-- | Leading coefficient of univariate polynomial.
uLead :: (Eq c, Num c) => UnivariatePoly c -> c
uLead (UnivariatePoly cs) =
  case trimU cs of
    [] -> 0
    xs -> last xs

-- | Tests if univariate polynomial is zero.
uIsZero :: (Eq c, Num c) => UnivariatePoly c -> Bool
uIsZero (UnivariatePoly cs) = null (trimU cs)

-- | Univariate polynomial addition.
uAdd :: Num c => UnivariatePoly c -> UnivariatePoly c -> UnivariatePoly c
uAdd (UnivariatePoly xs) (UnivariatePoly ys) =
  UnivariatePoly (addCoeffs xs ys)
  where
    addCoeffs [] b = b
    addCoeffs a [] = a
    addCoeffs (a:as) (b:bs) = (a + b) : addCoeffs as bs

-- | Univariate polynomial subtraction.
uSub :: Num c => UnivariatePoly c -> UnivariatePoly c -> UnivariatePoly c
uSub (UnivariatePoly xs) (UnivariatePoly ys) =
  UnivariatePoly (subCoeffs xs ys)
  where
    subCoeffs [] b = map negate b
    subCoeffs a [] = a
    subCoeffs (a:as) (b:bs) = (a - b) : subCoeffs as bs

-- | Univariate polynomial multiplication.
uMul :: (Eq c, Num c) => UnivariatePoly c -> UnivariatePoly c -> UnivariatePoly c
uMul (UnivariatePoly xs) (UnivariatePoly ys) =
  uPoly (mulCoeffs xs ys)
  where
    mulCoeffs [] _ = []
    mulCoeffs _ [] = []
    mulCoeffs (a:as) b =
      let term = map (a *) b
          rest = 0 : mulCoeffs as b
      in addLists term rest

    addLists [] b = b
    addLists a [] = a
    addLists (a:as) (b:bs) = (a + b) : addLists as bs

-- | Division with remainder over a field: a = q*b + r with deg(r) < deg(b).
uQuotRem :: forall c. (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c -> (UnivariatePoly c, UnivariatePoly c)
uQuotRem a b
  | uIsZero b = error "uQuotRem: division by zero polynomial"
  | otherwise =
      let (qCs, rCs) = go (trimU (unUPoly a)) []
      in (uPoly qCs, uPoly rCs)
  where
    degB = uDegree b
    leadB = uLead b
    invLeadB = recip leadB
    bCoeffs = trimU (unUPoly b)

    go currA accQ
      | length currA - 1 < degB = (accQ, currA)
      | otherwise =
          let degCurr = length currA - 1
              shift = degCurr - degB
              leadA = last currA
              factor = leadA * invLeadB
              termQ = replicate shift 0 ++ [factor]
              subtractor = unUPoly (uMul (UnivariatePoly termQ) (UnivariatePoly bCoeffs))
              nextA = trimU (unUPoly (uSub (UnivariatePoly currA) (UnivariatePoly subtractor)))
              nextAccQ = unUPoly (uAdd (UnivariatePoly accQ) (UnivariatePoly termQ))
          in go nextA nextAccQ

-- | Greatest common divisor (monic) over a field using Euclidean algorithm.
uGCD :: (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c -> UnivariatePoly c
uGCD a b
  | uIsZero b =
      let lc = uLead a
      in if lc == 0 then a else UnivariatePoly (map (/ lc) (unUPoly a))
  | otherwise =
      let (_, r) = uQuotRem a b
      in uGCD b r

-- | Extended Euclidean Algorithm: computes (g, s, t) such that s*a + t*b = g with g monic.
uExtendedGCD :: (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c -> (UnivariatePoly c, UnivariatePoly c, UnivariatePoly c)
uExtendedGCD a b
  | uIsZero b =
      let lc = uLead a
      in if lc == 0
         then (a, UnivariatePoly [1], UnivariatePoly [])
         else (UnivariatePoly (map (/ lc) (unUPoly a)), UnivariatePoly [recip lc], UnivariatePoly [])
  | otherwise =
      let (q, r) = uQuotRem a b
          (g, s1, t1) = uExtendedGCD b r
          s = t1
          t = uSub s1 (uMul q t1)
      in (g, s, t)

-- | Formal derivative of a univariate polynomial: d/dx(sum c_i x^i) = sum i*c_i x^{i-1}.
uDerivative :: Num c => UnivariatePoly c -> UnivariatePoly c
uDerivative (UnivariatePoly []) = UnivariatePoly []
uDerivative (UnivariatePoly (_:cs)) =
  UnivariatePoly [ fromIntegral i * c | (i, c) <- zip [1 :: Int ..] cs ]

-- | Subresultant Polynomial Remainder Sequence (Collins-Brown algorithm)
-- Computes the subresultant PRS over integers / rationals, preventing exponential coefficient growth.
uSubresultantPRS :: forall c. (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c -> [UnivariatePoly c]
uSubresultantPRS p1 p2
  | uDegree p1 < uDegree p2 = uSubresultantPRS p2 p1
  | uIsZero p2 = [p1]
  | otherwise = p1 : p2 : go p1 p2 1 1
  where
    go curr1 curr2 (!g :: c) (!h :: c)
      | uIsZero curr2 = []
      | otherwise =
          let d1 = uDegree curr1
              d2 = uDegree curr2
              delta = d1 - d2
              lc2 = uLead curr2
              (_, remPoly) = uQuotRem (prem lc2 delta curr1) curr2
              nextPoly = if uIsZero remPoly
                         then UnivariatePoly []
                         else UnivariatePoly (map (/ (g * (h ^ delta))) (unUPoly remPoly))
              nextG = lc2
              nextH = if delta == 1
                      then h
                      else (lc2 ^ delta) / (h ^ (delta - 1))
          in if uIsZero nextPoly
             then []
             else nextPoly : go curr2 nextPoly nextG nextH

    prem lc delta (UnivariatePoly cs) =
      UnivariatePoly (map (* (lc ^ (delta + 1))) cs)

-- | Resultant of two polynomials Res(f, g) via the subresultant PRS.
uResultant :: (Eq c, Fractional c) => UnivariatePoly c -> UnivariatePoly c -> c
uResultant f g
  | uIsZero f || uIsZero g = 0
  | uDegree f == 0 = uLead f ^ toInteger (uDegree g)
  | uDegree g == 0 = uLead g ^ toInteger (uDegree f)
  | otherwise =
      let prs = uSubresultantPRS f g
          lastPoly = last prs
      in if uDegree lastPoly > 0
         then 0  -- Shared factor exists, resultant is zero!
         else uLead lastPoly

-- | Discriminant of a univariate polynomial:
-- Disc(f) = (-1)^(n*(n-1)/2) * (1 / a_n) * Res(f, f').
uDiscriminant :: (Eq c, Fractional c) => UnivariatePoly c -> c
uDiscriminant f
  | n <= 0 = 0
  | otherwise =
      let sign = if even (n * (n - 1) `div` 2) then 1 else -1
          an = uLead f
          res = uResultant f (uDerivative f)
      in sign * (res / an)
  where
    n = uDegree f

-- | Yun's algorithm for square-free factorization of a monic polynomial:
-- Computes the list of square-free factors [p_1, p_2, ..., p_k]
-- such that f = p_1^1 * p_2^2 * ... * p_k^k where each p_i is square-free and pairwise coprime.
uSquareFree :: (Eq c, Fractional c) => UnivariatePoly c -> [UnivariatePoly c]
uSquareFree f
  | uDegree f <= 1 = [f]
  | otherwise =
      let fPrime = uDerivative f
          a0 = uGCD f fPrime
          b1 = fst (uQuotRem f a0)
          c1 = fst (uQuotRem fPrime a0)
          d1 = uSub c1 (uDerivative b1)
      in loop b1 d1
  where
    loop !b !d
      | uDegree b <= 0 = []
      | otherwise =
          let p = uGCD b d
              nextB = fst (uQuotRem b p)
              nextC = fst (uQuotRem d p)
              nextD = uSub nextC (uDerivative nextB)
          in p : loop nextB nextD

-- | Cauchy upper bound on the absolute value of all complex/real roots of f:
-- |x| <= 1 + max_{i=0}^{n-1} |a_i / a_n|.
uCauchyBound :: UnivariatePoly Rat -> Rat
uCauchyBound (UnivariatePoly cs) =
  case trimU cs of
    [] -> 0
    [_] -> 0
    coeffs ->
      let an = abs (last coeffs)
          ai = map abs (init coeffs)
          maxRatio = maximum [ c / an | c <- ai ]
      in 1 + maxRatio

instance (Show c, Eq c, Num c) => Show (UnivariatePoly c) where
  show (UnivariatePoly cs) =
    case trimU cs of
      [] -> "0"
      trimmed ->
        let formatTerm (0, c) = show c
            formatTerm (1, c) = show c ++ "*x"
            formatTerm (i, c) = show c ++ "*x^" ++ show i
            terms = [ formatTerm (i, c) | (i, c) <- zip [0 :: Int ..] trimmed, c /= 0 ]
        in if null terms then "0" else unwords terms
