{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.Quotient
  ( QuotientRing(..)
  , makeQuotientRing
  , quotientDimension
  , quotientVector
  , quotientPoly
  , quotientMultiply
  , multiplicationMatrix
  , hermiteTraceMatrix
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (SparseMonomial(..))
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  )
import Syzygy.Polynomial.Arithmetic (mulPoly)
import Syzygy.Polynomial.Division (normalForm)
import Syzygy.Linear.Vector (Vector(..), vecFromList, vecToList)
import Syzygy.Linear.Matrix (Matrix(..), matTrace)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Grobner.FGLM (isZeroDimensional, quotientStandardBasis)

-- | Represents the quotient algebra A = k[x_0, ..., x_{n-1}] / I
-- for a zero-dimensional ideal I.
-- Possesses finite vector space dimension D = dim_k(A).
data QuotientRing a = QuotientRing
  { qrBasis     :: ![SparseMonomial] -- ^ Standard monomial basis [b_0, ..., b_{D-1}]
  , qrGB        :: ![Poly a]         -- ^ Reduced Groebner basis of I
  , qrNumVars   :: !Int              -- ^ Number of variables n in the polynomial ring
  , qrOrder     :: !MonomialOrder    -- ^ Monomial order used for the Groebner basis
  } deriving (Eq, Show, Generic)

-- | Constructs a QuotientRing from a list of ideal generators.
-- Returns Nothing if the ideal is not zero-dimensional.
makeQuotientRing :: forall a. (Ord a, Fractional a)
                 => MonomialOrder
                 -> Int
                 -> [Poly a]
                 -> Maybe (QuotientRing a)
makeQuotientRing order numVars gens =
  let gb = reducedGrobnerBasis (buchberger gens)
  in if not (isZeroDimensional numVars gb)
     then Nothing
     else
       let stdBasis = quotientStandardBasis numVars gb
       in Just QuotientRing
            { qrBasis   = stdBasis
            , qrGB      = gb
            , qrNumVars = numVars
            , qrOrder   = order
            }

-- | Returns the vector space dimension D = dim_k(A).
quotientDimension :: QuotientRing a -> Int
quotientDimension qr = length (qrBasis qr)

-- | Computes the coordinate vector in k^D of a polynomial f modulo I
-- with respect to the standard monomial basis B = [b_0, ..., b_{D-1}].
quotientVector :: (Ord a, Fractional a) => QuotientRing a -> Poly a -> Vector a
quotientVector qr f =
  let nf = normalForm f (qrGB qr)
      basis = qrBasis qr
      termMap = Map.fromList [ (m, c) | Term c m <- polyTerms nf ]
      coords = [ Map.findWithDefault 0 b termMap | b <- basis ]
  in vecFromList coords
  where
    polyTerms (Poly terms _) = terms

-- | Reconstructs a polynomial in normal form from its coordinate vector in k^D.
quotientPoly :: (Eq a, Num a) => QuotientRing a -> Vector a -> Poly a
quotientPoly qr vec =
  let basis = qrBasis qr
      coords = vecToList vec
      terms = [ Term c b | (c, b) <- zip coords basis, c /= 0 ]
  in Poly terms (qrOrder qr)

-- | Multiplies two elements in the quotient ring: f * g mod I = NF(f * g, GB).
quotientMultiply :: (Ord a, Fractional a) => QuotientRing a -> Poly a -> Poly a -> Poly a
quotientMultiply qr f g =
  normalForm (mulPoly f g) (qrGB qr)

-- | Computes the D x D matrix representation of the multiplication operator M_f : A -> A
-- where M_f(g) = f * g mod I.
-- Column j of M_f is the coordinate vector of f * b_j mod I.
multiplicationMatrix :: forall a. (Ord a, Fractional a) => QuotientRing a -> Poly a -> Matrix a
multiplicationMatrix qr f =
  let basis = qrBasis qr
      d = length basis
      -- Column j is NF(f * b_j, GB)
      cols = [ vecToList (quotientVector qr (mulPoly f (Poly [Term 1 b] (qrOrder qr))))
             | b <- basis
             ]
      -- Transpose to row-major format
      rowMajor = [ [ cols !! colIdx !! rowIdx | colIdx <- [0 .. d - 1] ]
                 | rowIdx <- [0 .. d - 1]
                 ]
  in Matrix d d rowMajor

-- | Computes the Hermite trace quadratic form matrix H(A):
-- H_{i, j} = Trace(M_{b_i * b_j}).
-- By Hermite's Theorem:
-- rank(H) = number of distinct complex solutions of I.
-- signature(H) = number of distinct real solutions of I!
hermiteTraceMatrix :: forall a. (Ord a, Fractional a) => QuotientRing a -> Matrix a
hermiteTraceMatrix qr =
  let basis = qrBasis qr
      d = length basis
      entries = [ [ matTrace (multiplicationMatrix qr (mulPoly (Poly [Term 1 bi] (qrOrder qr)) (Poly [Term 1 bj] (qrOrder qr))))
                  | bj <- basis
                  ]
                | bi <- basis
                ]
  in Matrix d d entries
