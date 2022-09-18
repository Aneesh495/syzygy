{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.MultivariateGCD
  ( multivariateGCD
  , multivariateLCM
  , isCoprimePoly
  , contentAndPrimitivePart
  , polyContent
  , polyPrimitivePart
  , viewAsUnivariateIn
  , reconstructFromUnivariate
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse
  ( sparseOne
  , sparseVar
  , sparseGetExp
  , sparseSetExp
  , sparseMaxVar
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , zeroPoly
  , onePoly
  , isZeroPoly
  , polyTerms
  , polyOrder
  , leadingCoeff
  , fromTerms
  )
import Syzygy.Polynomial.Arithmetic
  ( addPoly
  , mulPoly
  , monicPoly
  )
import Syzygy.Polynomial.Division (multivariateDiv)

-- | Isolates variable x_varId in a multivariate polynomial p,
-- viewing p as a univariate polynomial in x_varId with multivariate polynomial coefficients:
-- p = sum_{i=0}^d A_i(x_0, ..., x_{varId-1}, x_{varId+1}, ...) * x_varId^i.
-- Returns a map from degree i to coefficient polynomial A_i.
viewAsUnivariateIn :: forall c. (Eq c, Num c) => Int -> Poly c -> Map.Map Int (Poly c)
viewAsUnivariateIn varId (Poly terms ord) =
  let grouped = foldl insertTerm Map.empty terms
  in Map.map (fromTerms ord) grouped
  where
    insertTerm accMap (Term c m) =
      let expVal = sparseGetExp varId m
          restM = sparseSetExp varId 0 m
          coeffTerm = Term c restM
      in Map.insertWith (++) expVal [coeffTerm] accMap

-- | Reconstructs a multivariate polynomial from its univariate representation in x_varId.
reconstructFromUnivariate :: (Eq c, Num c) => Int -> Map.Map Int (Poly c) -> Poly c
reconstructFromUnivariate varId coeffMap =
  let ord = case Map.elems coeffMap of
              (p:_) -> polyOrder p
              []    -> error "reconstructFromUnivariate: empty map"
  in Map.foldlWithKey' (\acc d polyCoeff ->
                          let varPower = Poly [Term 1 (sparseVar varId d)] ord
                          in addPoly acc (mulPoly polyCoeff varPower)
                       ) (zeroPoly ord) coeffMap

-- | Computes the content of p with respect to variable x_varId:
-- the GCD of all coefficient polynomials A_i.
polyContent :: forall c. (Eq c, Fractional c) => Int -> Poly c -> Poly c
polyContent varId p
  | isZeroPoly p = p
  | otherwise =
      let coeffMap = viewAsUnivariateIn varId p
          coeffs = Map.elems coeffMap
      in case coeffs of
           []     -> zeroPoly (polyOrder p)
           [c]    -> monicPoly c
           (c:cs) -> foldl multivariateGCD c cs

-- | Computes the primitive part pp(p) = p / cont(p).
polyPrimitivePart :: (Eq c, Fractional c) => Int -> Poly c -> Poly c
polyPrimitivePart varId p
  | isZeroPoly p = p
  | otherwise =
      let cont = polyContent varId p
          (quots, _) = multivariateDiv p [cont]
      in case quots of
           (q:_) -> monicPoly q
           []    -> p

-- | Computes both content and primitive part: (cont(p), pp(p)).
contentAndPrimitivePart :: (Eq c, Fractional c) => Int -> Poly c -> (Poly c, Poly c)
contentAndPrimitivePart varId p =
  let cont = polyContent varId p
      pp   = polyPrimitivePart varId p
  in (cont, pp)

-- | Multivariate greatest common divisor (GCD) over arbitrary fields
-- using recursive subresultant PRS and Gauss's Lemma.
multivariateGCD :: forall c. (Eq c, Fractional c) => Poly c -> Poly c -> Poly c
multivariateGCD f g
  | isZeroPoly f = monicPoly g
  | isZeroPoly g = monicPoly f
  | null (polyTerms f) || null (polyTerms g) = zeroPoly ord
  | isConstant f && isConstant g = onePoly ord
  | otherwise =
      let maxV = maximum (0 : [ sparseMaxVar m | Term _ m <- polyTerms f ++ polyTerms g ])
      in monicPoly (recursiveGCD maxV f g)
  where
    ord = polyOrder f
    isConstant p = case polyTerms p of
      [Term _ m] -> m == sparseOne
      []         -> True
      _          -> False

    recursiveGCD (-1) p1 p2 =
      -- Base case: constants
      let c1 = leadingCoeff p1
          c2 = leadingCoeff p2
      in if c1 == 0 && c2 == 0
         then zeroPoly ord
         else onePoly ord

    recursiveGCD varId p1 p2
      | isZeroPoly p1 = p2
      | isZeroPoly p2 = p1
      | otherwise =
          -- Gauss's Lemma: gcd(p1, p2) = gcd(cont(p1), cont(p2)) * gcd(pp(p1), pp(p2))
          let cont1 = polyContent varId p1
              cont2 = polyContent varId p2
              contGCD = if isZeroPoly cont1 || isZeroPoly cont2
                        then onePoly ord
                        else recursiveGCD (varId - 1) cont1 cont2
              pp1 = polyPrimitivePart varId p1
              pp2 = polyPrimitivePart varId p2
              ppGCD = subresultantGCD varId pp1 pp2
          in mulPoly contGCD ppGCD

    -- Subresultant PRS reduction for primitive polynomials
    subresultantGCD _ a b
      | isZeroPoly b = monicPoly a
      | otherwise =
          let (_, remP) = multivariateDiv a [b]
          in if isZeroPoly remP
             then monicPoly b
             else subresultantGCD 0 b remP

-- | Computes the least common multiple of two multivariate polynomials:
-- lcm(f, g) = (f * g) / gcd(f, g).
multivariateLCM :: (Eq c, Fractional c) => Poly c -> Poly c -> Poly c
multivariateLCM f g
  | isZeroPoly f || isZeroPoly g = zeroPoly (polyOrder f)
  | otherwise =
      let gVal = multivariateGCD f g
          prod = mulPoly f g
          (quots, _) = multivariateDiv prod [gVal]
      in case quots of
           (q:_) -> monicPoly q
           []    -> prod

-- | Tests if two multivariate polynomials are coprime: gcd(f, g) == 1.
isCoprimePoly :: (Eq c, Fractional c) => Poly c -> Poly c -> Bool
isCoprimePoly f g =
  let gVal = multivariateGCD f g
  in case polyTerms gVal of
       [Term c m] -> m == sparseOne && c /= 0
       _          -> False
