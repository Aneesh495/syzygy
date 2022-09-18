{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Division
  ( multivariateDiv
  , normalForm
  , fullReduction
  , reduceOneStep
  , pseudoDiv
  , dividesLeadingTerm
  ) where

import Syzygy.Monomial.Sparse
  ( sparseDegree
  , sparseDivides
  , sparseDiv
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , zeroPoly
  , isZeroPoly
  , leadingTerm
  , leadingMonomial
  , leadingCoeff
  , polyOrder
  )
import Syzygy.Polynomial.Arithmetic
  ( addPoly
  , subPoly
  , mulTermPoly
  , scalePoly
  )

-- | Tests if the leading term of g divides the leading term of f.
-- If true, returns Just (c, m) such that (c * m) * LT(g) = LT(f).
dividesLeadingTerm :: (Eq c, Fractional c) => Poly c -> Poly c -> Maybe (Term c)
dividesLeadingTerm f g =
  case (leadingTerm f, leadingTerm g) of
    (Just (Term cF mF), Just (Term cG mG))
      | sparseDivides mG mF ->
          let cQuot = cF / cG
              mQuot = sparseDiv mG mF
          in Just (Term cQuot mQuot)
      | otherwise -> Nothing
    _ -> Nothing

-- | Multivariate division algorithm (Cox-Little-O'Shea, Chapter 2).
-- Given a dividend f and an ordered tuple of divisors [f_1, ..., f_s],
-- computes ([a_1, ..., a_s], r) such that:
-- f = a_1*f_1 + ... + a_s*f_s + r
-- where no monomial in the remainder r is divisible by any of LT(f_1), ..., LT(f_s),
-- and multideg(a_i * f_i) <= multideg(f).
multivariateDiv :: forall c. (Eq c, Fractional c) => Poly c -> [Poly c] -> ([Poly c], Poly c)
multivariateDiv f divisors
  | null divisors = ([], f)
  | otherwise =
      let s = length divisors
          initQuotients = replicate s (zeroPoly ord)
      in go f initQuotients (zeroPoly ord)
  where
    ord = polyOrder f

    go !p !quotients !remainder
      | isZeroPoly p = (quotients, remainder)
      | otherwise =
          case tryDivision p divisors 0 of
            Just (i, termQuot) ->
              let subtrahend = mulTermPoly termQuot (divisors !! i)
                  pNext = subPoly p subtrahend
                  quotientsNext = updateList i (\q -> addPoly q (Poly [termQuot] ord)) quotients
              in go pNext quotientsNext remainder
            Nothing ->
              -- No divisor's leading term divides LT(p): move LT(p) into remainder
              case leadingTerm p of
                Just lt@(Term c m) ->
                  let pNext = subPoly p (Poly [lt] ord)
                      remainderNext = addPoly remainder (Poly [Term c m] ord)
                  in go pNext quotients remainderNext
                Nothing -> (quotients, remainder)

    tryDivision _ [] _ = Nothing
    tryDivision p (g:gs) !idx =
      case dividesLeadingTerm p g of
        Just termQuot -> Just (idx, termQuot)
        Nothing       -> tryDivision p gs (idx + 1)

    updateList 0 fn (x:xs) = fn x : xs
    updateList n fn (x:xs) = x : updateList (n - 1) fn xs
    updateList _ _ [] = []

-- | Computes the standard remainder r = f rem [f_1, ..., f_s] (Normal Form).
normalForm :: (Eq c, Fractional c) => Poly c -> [Poly c] -> Poly c
normalForm f divisors = snd (multivariateDiv f divisors)

-- | Reduces a polynomial by one step against a list of divisors.
-- Returns Just p' if a reduction occurred, or Nothing if f is irreducible.
reduceOneStep :: (Eq c, Fractional c) => Poly c -> [Poly c] -> Maybe (Poly c)
reduceOneStep p divisors = go divisors
  where
    go [] = Nothing
    go (g:gs) =
      case dividesLeadingTerm p g of
        Just termQuot ->
          let subtrahend = mulTermPoly termQuot g
          in Just (subPoly p subtrahend)
        Nothing -> go gs

-- | Full inter-reduction (tail reduction) of a polynomial against divisors.
-- Reduces all terms of p (not just the leading term), ensuring every monomial
-- in the resulting polynomial is completely irreducible against LT(divisors).
fullReduction :: (Eq c, Fractional c) => Poly c -> [Poly c] -> Poly c
fullReduction p divisors
  | isZeroPoly p = p
  | otherwise    = go p (zeroPoly ord)
  where
    ord = polyOrder p

    go !curr !acc
      | isZeroPoly curr = acc
      | otherwise =
          case findReducibleTerm curr divisors of
            Just (termQuot, g) ->
              let subtrahend = mulTermPoly termQuot g
                  nextCurr = subPoly curr subtrahend
              in go nextCurr acc
            Nothing ->
              case leadingTerm curr of
                Just lt@(Term c m) ->
                  let nextCurr = subPoly curr (Poly [lt] ord)
                      nextAcc = addPoly acc (Poly [Term c m] ord)
                  in go nextCurr nextAcc
                Nothing -> acc

    findReducibleTerm curr gs =
      case leadingTerm curr of
        Just _ ->
          let checkDivisors [] = Nothing
              checkDivisors (g:rest) =
                case dividesLeadingTerm curr g of
                  Just t  -> Just (t, g)
                  Nothing -> checkDivisors rest
          in checkDivisors gs
        Nothing -> Nothing

-- | Pseudo-division of f by g over a unique factorization domain (integral domain).
-- Multiplies f by LC(g)^(deg(f) - deg(g) + 1) to avoid fractions:
-- LC(g)^delta * f = q * g + r with deg(r) < deg(g).
-- Returns (factor, quotient, remainder).
pseudoDiv :: forall c. (Eq c, Integral c) => Poly c -> Poly c -> (c, Poly c, Poly c)
pseudoDiv f g
  | isZeroPoly g = error "pseudoDiv: division by zero polynomial"
  | isZeroPoly f = (1, zeroPoly ord, zeroPoly ord)
  | otherwise =
      let degF = termDeg f
          degG = termDeg g
      in if degF < degG
         then (1, zeroPoly ord, f)
         else
           let delta = degF - degG + 1
               lcG = leadingCoeff g
               factor = lcG ^ delta
               scaledF = scalePoly factor f
           in go scaledF (zeroPoly ord) lcG degG
  where
    ord = polyOrder f
    termDeg p = case leadingMonomial p of
      Just m  -> sparseDegree m
      Nothing -> -1

    go !currP !accQ !lcG !degG
      | isZeroPoly currP = (1, accQ, zeroPoly ord)
      | otherwise =
          case (leadingTerm currP, leadingTerm g) of
            (Just (Term cP mP), Just (Term _ mG))
              | sparseDivides mG mP ->
                  let mQuot = sparseDiv mG mP
                      termQ = Term (cP `quot` lcG) mQuot
                      subtrahend = mulTermPoly termQ g
                      nextP = subPoly currP subtrahend
                      nextQ = addPoly accQ (Poly [termQ] ord)
                  in go nextP nextQ lcG degG
              | otherwise -> (1, accQ, currP)
            _ -> (1, accQ, currP)
