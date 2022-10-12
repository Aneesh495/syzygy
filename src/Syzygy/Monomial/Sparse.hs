{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  , sparseDegree
  , sparseSupport
  , sparseMaxVar
  , sparseGetExp
  , sparseSetExp
  , sparseMul
  , sparseDivides
  , sparseDiv
  , sparseLCM
  , sparseGCD
  , sparseCoprime
  , toDenseList
  , fromDenseList
  , compareSparse
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder, compareMonomial)

-- | Arbitrary-variable sparse exponent vector represented as a strictly positive map
-- from 0-indexed variable ID to non-zero exponent.
-- The empty map corresponds to the identity monomial 1 = x_0^0 ... x_n^0.
newtype SparseMonomial = SparseMonomial { unSparseMonomial :: Map.Map Int Int }
  deriving (Eq, Ord, Generic)

-- | The identity monomial 1.
sparseOne :: SparseMonomial
sparseOne = SparseMonomial Map.empty

-- | Single variable monomial x_var^expVal.
sparseVar :: Int -> Int -> SparseMonomial
sparseVar varId expVal
  | expVal < 0 = error "sparseVar: exponent must be non-negative"
  | expVal == 0 = sparseOne
  | otherwise = SparseMonomial (Map.singleton varId expVal)

-- | Total degree of the sparse monomial: sum of all non-zero exponents.
sparseDegree :: SparseMonomial -> Int
sparseDegree (SparseMonomial m) = Map.foldl' (+) 0 m

-- | Returns the list of 0-indexed variable indices present with positive exponent.
sparseSupport :: SparseMonomial -> [Int]
sparseSupport (SparseMonomial m) = Map.keys m

-- | Returns the maximum variable index present, or -1 if the monomial is 1.
sparseMaxVar :: SparseMonomial -> Int
sparseMaxVar (SparseMonomial m)
  | Map.null m = -1
  | otherwise  = fst (Map.findMax m)

-- | Extracts the exponent of a specific variable x_varId.
sparseGetExp :: Int -> SparseMonomial -> Int
sparseGetExp varId (SparseMonomial m) = Map.findWithDefault 0 varId m

-- | Sets or updates the exponent of a variable x_varId.
sparseSetExp :: Int -> Int -> SparseMonomial -> SparseMonomial
sparseSetExp varId expVal (SparseMonomial m)
  | expVal < 0  = error "sparseSetExp: exponent must be non-negative"
  | expVal == 0 = SparseMonomial (Map.delete varId m)
  | otherwise   = SparseMonomial (Map.insert varId expVal m)

-- | Monomial multiplication: adds exponents component-wise.
sparseMul :: SparseMonomial -> SparseMonomial -> SparseMonomial
sparseMul (SparseMonomial m1) (SparseMonomial m2) =
  SparseMonomial (Map.unionWith (+) m1 m2)

-- | Tests if monomial u divides monomial v.
sparseDivides :: SparseMonomial -> SparseMonomial -> Bool
sparseDivides (SparseMonomial m1) (SparseMonomial m2) =
  Map.isSubmapOfBy (<=) m1 m2

-- | Exact monomial division v / u.
-- Precondition: u divides v.
sparseDiv :: SparseMonomial -> SparseMonomial -> SparseMonomial
sparseDiv (SparseMonomial mU) (SparseMonomial mV) =
  let diffMap = Map.differenceWith subExp mV mU
      subExp v u
        | v > u     = Just (v - u)
        | v == u    = Nothing
        | otherwise = error "sparseDiv: u does not divide v"
  in SparseMonomial diffMap

-- | Least common multiple (LCM) of two monomials: component-wise maximum of exponents.
sparseLCM :: SparseMonomial -> SparseMonomial -> SparseMonomial
sparseLCM (SparseMonomial m1) (SparseMonomial m2) =
  SparseMonomial (Map.unionWith max m1 m2)

-- | Greatest common divisor (GCD) of two monomials: component-wise minimum of exponents.
sparseGCD :: SparseMonomial -> SparseMonomial -> SparseMonomial
sparseGCD (SparseMonomial m1) (SparseMonomial m2) =
  SparseMonomial (Map.filter (> 0) (Map.intersectionWith min m1 m2))

-- | Tests if two monomials are coprime: gcd(u, v) == 1.
sparseCoprime :: SparseMonomial -> SparseMonomial -> Bool
sparseCoprime (SparseMonomial m1) (SparseMonomial m2) =
  Map.null (Map.intersection m1 m2)

-- | Converts a sparse monomial to a dense exponent list of length numVars.
toDenseList :: Int -> SparseMonomial -> [Int]
toDenseList numVars (SparseMonomial m) =
  [ Map.findWithDefault 0 i m | i <- [0 .. numVars - 1] ]

-- | Constructs a sparse monomial from a dense exponent list.
fromDenseList :: [Int] -> SparseMonomial
fromDenseList exps =
  let pairs = [ (idx, e) | (idx, e) <- zip [0 :: Int ..] exps, e > 0 ]
  in SparseMonomial (Map.fromList pairs)

-- | Compares two sparse monomials according to a given MonomialOrder in numVars variables.
compareSparse :: MonomialOrder -> Int -> SparseMonomial -> SparseMonomial -> Ordering
compareSparse order numVars m1 m2 =
  let n = max numVars (max (sparseMaxVar m1 + 1) (sparseMaxVar m2 + 1))
      dense1 = toDenseList n m1
      dense2 = toDenseList n m2
  in compareMonomial order dense1 dense2

instance Semigroup SparseMonomial where
  (<>) = sparseMul

instance Monoid SparseMonomial where
  mempty = sparseOne

instance Show SparseMonomial where
  show (SparseMonomial m)
    | Map.null m = "1"
    | otherwise =
        let formatTerm (varId, e) =
              "x" ++ show varId ++ (if e == 1 then "" else "^" ++ show e)
        in unwords [ formatTerm p | p <- Map.toAscList m ]
