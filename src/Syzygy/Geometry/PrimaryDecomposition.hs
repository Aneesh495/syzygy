{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Geometry.PrimaryDecomposition
  ( PrimaryComponent(..)
  , primaryDecompositionZeroDim
  , associatedPrimesZeroDim
  , isPrimeIdealZeroDim
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , reorderPoly
  )
import Syzygy.Polynomial.MultivariateGCD (multivariateGCD)
import Syzygy.Polynomial.Division (multivariateDiv)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Geometry.Ideal
  ( Ideal(..)
  , idealFromGenerators
  , idealEqual
  , idealSum
  )
import Syzygy.Geometry.Radical (radicalZeroDimensional, squareFreePart)

-- | Represents a primary component Q_i and its associated prime ideal p_i = sqrt(Q_i).
data PrimaryComponent a = PrimaryComponent
  { pcPrimaryIdeal :: !(Ideal a) -- ^ The primary ideal Q_i
  , pcPrimeIdeal   :: !(Ideal a) -- ^ The associated prime ideal p_i = sqrt(Q_i)
  , pcIsIsolated   :: !Bool      -- ^ True if p_i is a minimal prime, False if embedded
  } deriving (Eq, Show, Generic)

-- | Computes the primary decomposition of a zero-dimensional ideal I:
-- I = Q_1 cap Q_2 cap ... cap Q_r
-- using the Gianni-Trager-Zacharias (GTZ 1988) algorithm via univariate projections.
primaryDecompositionZeroDim :: forall a. (Ord a, Fractional a)
                            => Int
                            -> Ideal a
                            -> [PrimaryComponent a]
primaryDecompositionZeroDim numVars initialIdeal =
  let components = decomposeAlongVars 0 initialIdeal
  in [ let rad = radicalZeroDimensional numVars q
       in PrimaryComponent
            { pcPrimaryIdeal = q
            , pcPrimeIdeal   = rad
            , pcIsIsolated   = True
            }
     | q <- components
     ]
  where
    decomposeAlongVars varId currentIdeal
      | varId >= numVars = [currentIdeal]
      | otherwise =
          let univarPoly = findUnivariate currentIdeal varId
              factors = factorSquareFree univarPoly varId
          in case factors of
               []  -> [currentIdeal]
               [_] -> decomposeAlongVars (varId + 1) currentIdeal
               fs  ->
                 let subIdeals = [ idealSum currentIdeal (idealFromGenerators (idealOrder currentIdeal) [f])
                                 | f <- fs
                                 ]
                 in concatMap (decomposeAlongVars (varId + 1)) subIdeals

    findUnivariate ideal vId =
      let elimOrder = Lex
          reorderedGens = [ reorderPoly elimOrder g | g <- idealGens ideal ]
          gb = reducedGrobnerBasis (buchberger reorderedGens)
          cands = filter (isPureInVar vId) gb
      in case cands of
           (cand:_) -> reorderPoly (idealOrder ideal) cand
           []       -> Poly [Term 1 (sparseVar vId 1)] (idealOrder ideal)

    isPureInVar vId (Poly terms _) =
      all (\(Term _ (SparseMonomial m)) ->
            case Map.keys m of
              []  -> True
              [v] -> v == vId
              _   -> False
          ) terms

    -- Splits a polynomial into coprime square-free parts using derivatives and GCDs
    factorSquareFree poly vId =
      let sqFree = squareFreePart vId poly
          g = multivariateGCD poly sqFree
      in if poly == sqFree || g == poly
         then [poly]
         else case multivariateDiv poly [g] of
                ([q], _) -> [g, q]
                _        -> [poly]

-- | Returns the list of associated prime ideals {p_1, ..., p_r} of a zero-dimensional ideal.
associatedPrimesZeroDim :: (Ord a, Fractional a) => Int -> Ideal a -> [Ideal a]
associatedPrimesZeroDim numVars ideal =
  map pcPrimeIdeal (primaryDecompositionZeroDim numVars ideal)

-- | Checks whether a zero-dimensional ideal is prime:
-- I is prime iff I is primary and radical.
isPrimeIdealZeroDim :: (Ord a, Fractional a) => Int -> Ideal a -> Bool
isPrimeIdealZeroDim numVars ideal =
  case primaryDecompositionZeroDim numVars ideal of
    [comp] -> idealEqual ideal (pcPrimeIdeal comp)
    _      -> False
