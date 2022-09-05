{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Monomial.Order
  ( MonomialOrder(..)
  , compareMonomial
  , cmpLex
  , cmpDegLex
  , cmpDegRevLex
  , cmpBlock
  , cmpMatrix
  , monomialDegree
  , monomialSupport
  ) where

import GHC.Generics (Generic)

-- | Monomial ordering for multivariate polynomials.
-- Satisfies well-ordering and compatibility with polynomial multiplication:
-- u > v implies u * w > v * w for all monomials w.
data MonomialOrder
  = Lex                                       -- ^ Pure Lexicographical order
  | DegLex                                    -- ^ Degree Lexicographical order
  | DegRevLex                                 -- ^ Graded Reverse Lexicographical order (Grevlex)
  | BlockOrder !Int !MonomialOrder !MonomialOrder -- ^ Elimination block order splitting at variable k
  | MatrixOrder ![[Int]]                      -- ^ Weight matrix order by Robbiano's theorem
  deriving (Eq, Show, Read, Generic)

-- | Computes the total degree of an exponent vector sum(alpha_i).
monomialDegree :: [Int] -> Int
monomialDegree = sum

-- | Returns the 0-indexed variable indices that have strictly positive exponents.
monomialSupport :: [Int] -> [Int]
monomialSupport exps = [ idx | (idx, e) <- zip [0 :: Int ..] exps, e > 0 ]

-- | Pure Lexicographical Order:
-- alpha >_lex beta iff the leftmost non-zero entry of alpha - beta is positive.
cmpLex :: [Int] -> [Int] -> Ordering
cmpLex [] [] = EQ
cmpLex [] (y:ys) = if y > 0 then LT else cmpLex [] ys
cmpLex (x:xs) [] = if x > 0 then GT else cmpLex xs []
cmpLex (x:xs) (y:ys) =
  case compare x y of
    EQ -> cmpLex xs ys
    other -> other

-- | Degree Lexicographical Order:
-- Compares total degree first; breaks ties using Lex.
cmpDegLex :: [Int] -> [Int] -> Ordering
cmpDegLex alpha beta =
  case compare (monomialDegree alpha) (monomialDegree beta) of
    EQ -> cmpLex alpha beta
    other -> other

-- | Graded Reverse Lexicographical Order (Grevlex):
-- Compares total degree first. If degrees are equal,
-- alpha >_grevlex beta iff the RIGHTMOST non-zero entry of alpha - beta is NEGATIVE!
-- This is the most computationally efficient order for Groebner bases.
cmpDegRevLex :: [Int] -> [Int] -> Ordering
cmpDegRevLex alpha beta =
  case compare (monomialDegree alpha) (monomialDegree beta) of
    EQ -> revLex (reverse alpha) (reverse beta)
    other -> other
  where
    revLex [] [] = EQ
    revLex [] (y:ys) = if y > 0 then GT else revLex [] ys
    revLex (x:xs) [] = if x > 0 then LT else revLex xs []
    revLex (x:xs) (y:ys) =
      case compare x y of
        EQ -> revLex xs ys
        LT -> GT  -- Note the reversal: smaller entry at the right means larger monomial!
        GT -> LT

-- | Block (Elimination) Order:
-- Splits variables into two blocks: [x_0 .. x_{k-1}] and [x_k .. x_n].
-- First compares the first block using order1; breaks ties with order2 on the second block.
-- Crucial for computing elimination ideals I \cap k[x_k .. x_n].
cmpBlock :: Int -> MonomialOrder -> MonomialOrder -> [Int] -> [Int] -> Ordering
cmpBlock k order1 order2 alpha beta =
  let (alpha1, alpha2) = splitAt k alpha
      (beta1, beta2)   = splitAt k beta
  in case compareMonomial order1 alpha1 beta1 of
       EQ -> compareMonomial order2 alpha2 beta2
       other -> other

-- | Weight Matrix Order:
-- Given an invertible m x n weight matrix M, alpha >_M beta iff M*alpha >_lex M*beta.
cmpMatrix :: [[Int]] -> [Int] -> [Int] -> Ordering
cmpMatrix [] _ _ = EQ
cmpMatrix (row:rows) alpha beta =
  let wA = sum (zipWith (*) row alpha)
      wB = sum (zipWith (*) row beta)
  in case compare wA wB of
       EQ -> cmpMatrix rows alpha beta
       other -> other

-- | Unified dispatcher for comparing two exponent vectors under any MonomialOrder.
compareMonomial :: MonomialOrder -> [Int] -> [Int] -> Ordering
compareMonomial Lex alpha beta = cmpLex alpha beta
compareMonomial DegLex alpha beta = cmpDegLex alpha beta
compareMonomial DegRevLex alpha beta = cmpDegRevLex alpha beta
compareMonomial (BlockOrder k o1 o2) alpha beta = cmpBlock k o1 o2 alpha beta
compareMonomial (MatrixOrder m) alpha beta = cmpMatrix m alpha beta
