{-# LANGUAGE BangPatterns #-}

module Syzygy.Monomial.Operations
  ( sMonomialMultipliers
  , monomialShift
  , monomialDegreeShift
  , monomialVarCount
  , areDisjointSupport
  ) where

import Syzygy.Monomial.Sparse
  ( SparseMonomial
  , sparseLCM
  , sparseDiv
  , sparseSupport
  , sparseMaxVar
  , sparseOne
  , sparseVar
  , sparseMul
  , sparseGetExp
  , sparseSetExp
  )

-- | Computes the S-monomial multipliers for two leading monomials m1 and m2.
-- Returns (u1, u2) such that:
-- u1 * m1 = u2 * m2 = lcm(m1, m2).
sMonomialMultipliers :: SparseMonomial -> SparseMonomial -> (SparseMonomial, SparseMonomial)
sMonomialMultipliers m1 m2 =
  let lcmM = sparseLCM m1 m2
      u1   = sparseDiv m1 lcmM
      u2   = sparseDiv m2 lcmM
  in (u1, u2)

-- | Shifts all variable indices in a monomial by a constant offset delta:
-- x_i -> x_{i + delta}.
monomialShift :: Int -> SparseMonomial -> SparseMonomial
monomialShift delta m
  | delta == 0 = m
  | otherwise =
      let supp = sparseSupport m
      in foldl (\acc i ->
                  let e = sparseGetExp i m
                  in sparseSetExp (i + delta) e acc
               ) sparseOne supp

-- | Computes the effective variable count needed to represent a collection of monomials.
monomialVarCount :: [SparseMonomial] -> Int
monomialVarCount [] = 0
monomialVarCount ms = maximum (0 : map (\m -> sparseMaxVar m + 1) ms)

-- | Tests if two monomials have completely disjoint variable supports.
-- If disjoint, their leading terms are coprime and S(f, g) reduces to zero by Buchberger's first criterion.
areDisjointSupport :: SparseMonomial -> SparseMonomial -> Bool
areDisjointSupport m1 m2 =
  let s1 = sparseSupport m1
      s2 = sparseSupport m2
  in null [ x | x <- s1, x `elem` s2 ]

-- | Homogenizes a monomial of degree d with respect to homogenization variable x_0:
-- x_0^(targetDeg - d) * m.
monomialDegreeShift :: Int -> Int -> SparseMonomial -> SparseMonomial
monomialDegreeShift targetDeg currentDeg m
  | targetDeg < currentDeg = error "monomialDegreeShift: targetDeg < currentDeg"
  | targetDeg == currentDeg = m
  | otherwise =
      let diffDeg = targetDeg - currentDeg
      in sparseMul (sparseVar 0 diffDeg) (monomialShift 1 m)
