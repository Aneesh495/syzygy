{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Invariant.Reynolds
  ( reynoldsOperator
  , actMatrixOnPoly
  , isInvariantPolynomial
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , zeroPoly
  )
import Syzygy.Polynomial.Arithmetic
  ( addPoly
  , mulPoly
  , powPoly
  , scalePoly
  )
import Syzygy.Linear.Matrix (Matrix(..))
import Syzygy.Invariant.FiniteGroup
  ( FiniteGroup(..)
  , groupOrder
  )

-- | Computes the action of a matrix g in GL(n, k) on a polynomial f:
-- (g * f)(x) = f(g^T * x).
-- Maps each variable x_i to the linear combination sum_{j=0}^{n-1} g_{j, i} * x_j.
actMatrixOnPoly :: forall a. (Eq a, Num a) => Matrix a -> Poly a -> Poly a
actMatrixOnPoly (Matrix numR numC entries) (Poly terms ord) =
  let -- Construct linear polynomial images of each variable x_i under g
      varImages = [ Poly [ Term ((entries !! j) !! i) (sparseVar j 1)
                         | j <- [0 .. numR - 1]
                         , (entries !! j) !! i /= 0
                         ] ord
                  | i <- [0 .. numC - 1]
                  ]
      evalTransformedTerm c (SparseMonomial m) =
        let factors = [ case Map.lookup v m of
                          Just expVal -> powPoly (varImages !! v) (fromIntegral expVal)
                          Nothing     -> Poly [Term 1 (sparseVar 0 0)] ord
                      | v <- [0 .. numC - 1]
                      , Map.member v m
                      ]
            baseProduct = foldl mulPoly (Poly [Term 1 (sparseVar 0 0)] ord) factors
        in scalePoly c baseProduct
      -- Transform each term c * x_0^{e_0} ... x_{n-1}^{e_{n-1}}
      transTerms = [ evalTransformedTerm c m | Term c m <- terms ]
  in foldl addPoly (zeroPoly ord) transTerms

-- | Computes the Reynolds operator R_G(f) = (1 / |G|) * sum_{g in G} (g * f).
-- Projects any polynomial f in k[x_0, ..., x_{n-1}] onto the invariant subring k[x]^G.
reynoldsOperator :: forall a. (Ord a, Fractional a) => FiniteGroup a -> Poly a -> Poly a
reynoldsOperator fg f =
  let gList = fgElements fg
      orderG = groupOrder fg
  in if orderG <= 0
     then f
     else
       let ord = case f of Poly _ o -> o
           invSum = foldl addPoly (zeroPoly ord) [ actMatrixOnPoly g f | g <- gList ]
           scaleFactor = recip (fromIntegral orderG)
       in scalePoly scaleFactor invSum

-- | Checks whether a polynomial f is invariant under the group G:
-- f in k[x]^G iff R_G(f) == f.
isInvariantPolynomial :: (Ord a, Fractional a) => FiniteGroup a -> Poly a -> Bool
isInvariantPolynomial fg f =
  reynoldsOperator fg f == f
