{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module Syzygy.Polynomial.Arithmetic
  ( addPoly
  , subPoly
  , scalePoly
  , mulTermPoly
  , mulPoly
  , powPoly
  , monicPoly
  , polyDerivative
  , evalPoly
  , substitutePoly
  , substituteAll
  ) where

import Syzygy.Monomial.Sparse
  ( sparseOne
  , sparseMul
  , sparseGetExp
  , sparseSetExp
  , toDenseList
  , compareSparse
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , zeroPoly
  , onePoly
  , constantPoly
  , polyTerms
  , polyOrder
  , leadingCoeff
  )

-- | Adds two polynomials with the same MonomialOrder in O(N_1 + N_2) time
-- by merging their sorted term lists.
addPoly :: forall c. (Eq c, Num c) => Poly c -> Poly c -> Poly c
addPoly (Poly terms1 ord1) (Poly terms2 ord2)
  | ord1 /= ord2 = error "addPoly: polynomials have distinct MonomialOrders"
  | otherwise    =
      let numVars = max (findMaxVar terms1) (findMaxVar terms2)
          merged = merge terms1 terms2 numVars
      in Poly merged ord1
  where
    findMaxVar ts = maximum (0 : [ maximum (0 : toDenseList 100 m) | Term _ m <- ts ])

    merge [] ys _ = ys
    merge xs [] _ = xs
    merge (t1@(Term c1 m1) : rest1) (t2@(Term c2 m2) : rest2) nVars =
      case compareSparse ord1 nVars m1 m2 of
        GT -> t1 : merge rest1 (t2 : rest2) nVars
        LT -> t2 : merge (t1 : rest1) rest2 nVars
        EQ ->
          let c = c1 + c2
          in if c == 0
             then merge rest1 rest2 nVars
             else Term c m1 : merge rest1 rest2 nVars

-- | Subtracts the second polynomial from the first: p1 - p2.
subPoly :: (Eq c, Num c) => Poly c -> Poly c -> Poly c
subPoly p1 (Poly terms2 ord2) =
  let negTerms = [ Term (-c) m | Term c m <- terms2 ]
  in addPoly p1 (Poly negTerms ord2)

-- | Multiplies a polynomial by a scalar factor c.
scalePoly :: (Eq c, Num c) => c -> Poly c -> Poly c
scalePoly 0 p = zeroPoly (polyOrder p)
scalePoly 1 p = p
scalePoly c (Poly terms ord) =
  let scaled = [ Term (c * coeff) m | Term coeff m <- terms ]
  in Poly scaled ord

-- | Multiplies a polynomial by a single monomial term (c * m).
mulTermPoly :: (Eq c, Num c) => Term c -> Poly c -> Poly c
mulTermPoly (Term 0 _) p = zeroPoly (polyOrder p)
mulTermPoly (Term c1 m1) (Poly terms ord)
  | c1 == 1 && m1 == sparseOne = Poly terms ord
  | otherwise =
      let newTerms = [ Term (c1 * c2) (sparseMul m1 m2) | Term c2 m2 <- terms ]
      in Poly newTerms ord

-- | Multiplies two polynomials p1 * p2 using Cauchy product and k-way merging.
mulPoly :: forall c. (Eq c, Num c) => Poly c -> Poly c -> Poly c
mulPoly (Poly [] ord) _ = zeroPoly ord
mulPoly _ (Poly [] ord) = zeroPoly ord
mulPoly p1 (Poly terms2 _) =
  foldl (\acc t2 -> addPoly acc (mulTermPoly t2 p1)) (zeroPoly (polyOrder p1)) terms2

-- | Exponentiation of a polynomial by a non-negative integer using binary squaring.
powPoly :: (Eq c, Num c) => Poly c -> Integer -> Poly c
powPoly p 0 = onePoly (polyOrder p)
powPoly p 1 = p
powPoly p n
  | n < 0 = error "powPoly: negative exponent requires rational field inversion"
  | otherwise = go p n (onePoly (polyOrder p))
  where
    go !_ 0 !acc = acc
    go !b !e !acc =
      let !acc' = if odd e then mulPoly acc b else acc
          !b'   = mulPoly b b
          !e'   = e `div` 2
      in go b' e' acc'

-- | Normalizes a polynomial to be monic (leading coefficient equals 1)
-- by dividing all coefficients by the leading coefficient LC(p).
monicPoly :: (Eq c, Fractional c) => Poly c -> Poly c
monicPoly p
  | null (polyTerms p) = p
  | otherwise =
      let lc = leadingCoeff p
          invLC = recip lc
      in scalePoly invLC p

-- | Computes the formal partial derivative d(p)/d(x_varId).
polyDerivative :: forall c. (Eq c, Num c) => Int -> Poly c -> Poly c
polyDerivative varId (Poly terms ord) =
  let derivTerms = concatMap diffTerm terms
  in Poly derivTerms ord
  where
    diffTerm (Term c m) =
      let expVal = sparseGetExp varId m
      in if expVal == 0
         then []
         else
           let newCoeff = c * fromIntegral expVal
               newM = sparseSetExp varId (expVal - 1) m
           in [Term newCoeff newM]

-- | Evaluates a multivariate polynomial at a coordinate vector [v_0, v_1, ...].
evalPoly :: forall c. Num c => [c] -> Poly c -> c
evalPoly coords (Poly terms _) =
  sum [ evalTerm t | t <- terms ]
  where
    evalTerm (Term c m) =
      let exps = toDenseList (length coords) m
          factors = [ val ^ e | (val, e) <- zip coords exps, e > 0 ]
      in c * product factors

-- | Substitutes variable x_varId with another polynomial replacementPoly in p.
substitutePoly :: (Eq c, Num c) => Int -> Poly c -> Poly c -> Poly c
substitutePoly varId replacementPoly (Poly terms ord) =
  foldl addPoly (zeroPoly ord) [ evalTerm t | t <- terms ]
  where
    evalTerm (Term c m) =
      let expVal = sparseGetExp varId m
          restM = sparseSetExp varId 0 m
          restPoly = Poly [Term c restM] ord
          substPower = powPoly replacementPoly (fromIntegral expVal)
      in mulPoly restPoly substPower

-- | Simultaneously substitutes multiple variables with corresponding polynomials.
substituteAll :: (Eq c, Num c) => [(Int, Poly c)] -> Poly c -> Poly c
substituteAll substitutions initialPoly =
  foldl (\p (varId, repl) -> substitutePoly varId repl p) initialPoly substitutions

instance (Eq c, Num c) => Num (Poly c) where
  (+) = addPoly
  (-) = subPoly
  (*) = mulPoly
  negate p = scalePoly (-1) p
  abs x = x
  signum _ = onePoly (error "Poly signum: undefined MonomialOrder")
  fromInteger n = constantPoly (fromInteger n) (error "Poly fromInteger: undefined MonomialOrder")
