{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , zeroPoly
  , onePoly
  , constantPoly
  , varPoly
  , termPoly
  , fromTerms
  , polyTerms
  , polyOrder
  , isZeroPoly
  , leadingTerm
  , leadingMonomial
  , leadingCoeff
  , polyDegree
  , polySupport
  , polyVarCount
  , reorderPoly
  ) where

import Data.List (sortBy)
import GHC.Generics (Generic)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  , sparseDegree
  , sparseSupport
  , sparseMaxVar
  , compareSparse
  )

-- | A single polynomial term c * m where c is a non-zero coefficient and m is a sparse monomial.
data Term c = Term !c !SparseMonomial
  deriving (Eq, Generic)

-- | Canonical multivariate polynomial in k[x_0, x_1, ...] represented as a list of terms
-- sorted strictly in DESCENDING order with respect to the specified MonomialOrder.
-- All coefficients are guaranteed to be non-zero (c /= 0).
data Poly c = Poly ![Term c] !MonomialOrder
  deriving (Generic)

-- | Zero polynomial (empty term list).
zeroPoly :: MonomialOrder -> Poly c
zeroPoly ord = Poly [] ord

-- | Identity polynomial 1.
onePoly :: Num c => MonomialOrder -> Poly c
onePoly ord = Poly [Term 1 sparseOne] ord

-- | Constant polynomial c.
constantPoly :: (Eq c, Num c) => c -> MonomialOrder -> Poly c
constantPoly c ord
  | c == 0    = zeroPoly ord
  | otherwise = Poly [Term c sparseOne] ord

-- | Single variable polynomial x_varId.
varPoly :: (Eq c, Num c) => Int -> MonomialOrder -> Poly c
varPoly varId ord =
  Poly [Term 1 (sparseVar varId 1)] ord

-- | Monomial term polynomial c * m.
termPoly :: (Eq c, Num c) => c -> SparseMonomial -> MonomialOrder -> Poly c
termPoly c m ord
  | c == 0    = zeroPoly ord
  | otherwise = Poly [Term c m] ord

-- | Constructs a canonical polynomial from an arbitrary list of terms.
-- Sorts terms descending by MonomialOrder, merges identical monomials by adding coefficients,
-- and filters out all terms with zero coefficients.
fromTerms :: forall c. (Eq c, Num c) => MonomialOrder -> [Term c] -> Poly c
fromTerms ord rawTerms =
  let nonZero = filter (\(Term c _) -> c /= 0) rawTerms
      numVars = maximum (0 : [ sparseMaxVar m + 1 | Term _ m <- nonZero ])
      sorted = sortBy (\(Term _ m1) (Term _ m2) ->
                         case compareSparse ord numVars m1 m2 of
                           LT -> GT  -- Descending order: largest term first!
                           GT -> LT
                           EQ -> EQ
                      ) nonZero
      merged = mergeSorted sorted
  in Poly merged ord
  where
    mergeSorted [] = []
    mergeSorted [t] = [t]
    mergeSorted (Term c1 m1 : Term c2 m2 : rest)
      | m1 == m2 =
          let sumC = c1 + c2
          in if sumC == 0
             then mergeSorted rest
             else mergeSorted (Term sumC m1 : rest)
      | otherwise =
          Term c1 m1 : mergeSorted (Term c2 m2 : rest)

-- | Returns the list of non-zero terms in descending order.
polyTerms :: Poly c -> [Term c]
polyTerms (Poly terms _) = terms

-- | Returns the monomial ordering associated with the polynomial.
polyOrder :: Poly c -> MonomialOrder
polyOrder (Poly _ ord) = ord

-- | Tests if a polynomial is identically zero.
isZeroPoly :: Poly c -> Bool
isZeroPoly (Poly terms _) = null terms

-- | Returns the leading term (head of the term list), or Nothing if zero.
leadingTerm :: Poly c -> Maybe (Term c)
leadingTerm (Poly [] _) = Nothing
leadingTerm (Poly (t:_) _) = Just t

-- | Returns the leading monomial LM(f), or Nothing if zero.
leadingMonomial :: Poly c -> Maybe SparseMonomial
leadingMonomial p = case leadingTerm p of
  Just (Term _ m) -> Just m
  Nothing         -> Nothing

-- | Returns the leading coefficient LC(f), or 0 if zero.
leadingCoeff :: Num c => Poly c -> c
leadingCoeff p = case leadingTerm p of
  Just (Term c _) -> c
  Nothing         -> 0

-- | Total degree of the polynomial: max degree over all terms (or -1 if zero).
polyDegree :: Poly c -> Int
polyDegree (Poly [] _) = -1
polyDegree (Poly terms _) =
  maximum [ sparseDegree m | Term _ m <- terms ]

-- | Returns all variable indices present in any term of the polynomial.
polySupport :: Poly c -> [Int]
polySupport (Poly terms _) =
  foldl (\acc (Term _ m) -> acc ++ filter (`notElem` acc) (sparseSupport m)) [] terms

-- | Effective variable count needed to represent this polynomial.
polyVarCount :: Poly c -> Int
polyVarCount (Poly [] _) = 0
polyVarCount (Poly terms _) =
  maximum (0 : [ sparseMaxVar m + 1 | Term _ m <- terms ])

-- | Converts a polynomial to a new MonomialOrder, re-sorting its terms canonically.
reorderPoly :: (Eq c, Num c) => MonomialOrder -> Poly c -> Poly c
reorderPoly newOrd p@(Poly terms currentOrd)
  | newOrd == currentOrd = p
  | otherwise            = fromTerms newOrd terms

instance Eq c => Eq (Poly c) where
  (Poly t1 ord1) == (Poly t2 ord2)
    | ord1 == ord2 = t1 == t2
    | otherwise    = error "Poly (==): comparing polynomials with distinct MonomialOrders"

instance Show c => Show (Term c) where
  show (Term c m)
    | m == sparseOne = show c
    | otherwise      = show c ++ "*" ++ show m

instance (Show c, Eq c, Ord c, Num c) => Show (Poly c) where
  show (Poly [] _) = "0"
  show (Poly (t0:ts) _) =
    let formatFirst (Term c m)
          | m == sparseOne = show c
          | c == 1         = show m
          | c == -1        = "-" ++ show m
          | otherwise      = show c ++ "*" ++ show m
        formatSubsequent (Term c m)
          | m == sparseOne = if c > 0 then " + " ++ show c else " - " ++ show (abs c)
          | c == 1         = " + " ++ show m
          | c == -1        = " - " ++ show m
          | c > 0          = " + " ++ show c ++ "*" ++ show m
          | otherwise      = " - " ++ show (abs c) ++ "*" ++ show m
    in formatFirst t0 ++ concatMap formatSubsequent ts
