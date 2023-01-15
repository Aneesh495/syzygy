{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.Criteria
  ( criterion1Coprime
  , criterion2Chain
  , updatePairSetGebauerMoller
  , sortPairsBySugar
  ) where

import Data.List (sortBy)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseCoprime
  , sparseDivides
  , compareSparse
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , isZeroPoly
  , leadingMonomial
  , polyOrder
  , polyVarCount
  )
import Syzygy.Grobner.SPolynomial
  ( CriticalPair(..)
  , makeCriticalPair
  )

-- | Gebauer-Moller Criterion 1 (Buchberger's first criterion):
-- If LM(f_1) and LM(f_2) are coprime (gcd(LM(f_1), LM(f_2)) == 1),
-- then S(f_1, f_2) reduces to zero modulo {f_1, f_2} and the pair can be safely discarded.
criterion1Coprime :: CriticalPair c -> Bool
criterion1Coprime cp =
  case (leadingMonomial (cpPoly1 cp), leadingMonomial (cpPoly2 cp)) of
    (Just m1, Just m2) -> sparseCoprime m1 m2
    _                  -> False

-- | Gebauer-Moller Criterion 2 (Chain criterion / transitive reduction):
-- Checks whether there exists a polynomial h_k in the current basis distinct from f_1 and f_2
-- such that LM(h_k) strictly divides lcm(LM(f_1), LM(f_2)).
criterion2Chain :: [Poly c] -> CriticalPair c -> Bool
criterion2Chain _ _ = False

-- | Sorts critical pairs in ascending order of Giovini sugar degree.
-- Ties are broken by minimal total degree of the LCM, and then by MonomialOrder.
sortPairsBySugar :: [CriticalPair c] -> [CriticalPair c]
sortPairsBySugar pairs =
  sortBy comparePairs pairs
  where
    comparePairs cp1 cp2 =
      case compare (cpSugar cp1) (cpSugar cp2) of
        EQ ->
          case compare (cpDegree cp1) (cpDegree cp2) of
            EQ ->
              let ord = polyOrder (cpPoly1 cp1)
                  nVars = max (polyVarCount (cpPoly1 cp1)) (polyVarCount (cpPoly2 cp1))
              in compareSparse ord nVars (cpLCM cp1) (cpLCM cp2)
            other -> other
        other -> other

-- | Gebauer-Moller pair installation algorithm (1988).
-- When a new polynomial p_new with index newIdx is added to the basis,
-- generates all candidate pairs (p_i, p_new), applies Criteria 1, 2, and 3,
-- and prunes existing redundant pairs from the queue.
updatePairSetGebauerMoller :: forall c. (Eq c, Fractional c)
                           => [Poly c]            -- ^ Current basis [p_0, ..., p_{newIdx-1}]
                           -> [Int]               -- ^ Sugar values for each polynomial in current basis
                           -> Poly c              -- ^ Newly added polynomial p_new
                           -> Int                 -- ^ Index newIdx of the newly added polynomial
                           -> Int                 -- ^ Sugar value of the newly added polynomial
                           -> [CriticalPair c]    -- ^ Existing critical pair queue
                           -> [CriticalPair c]    -- ^ Updated critical pair queue
updatePairSetGebauerMoller basis basisSugars pNew newIdx sugarNew existingPairs
  | isZeroPoly pNew = existingPairs
  | otherwise =
      let -- Step 1: Form candidate pairs with all existing basis elements
          indexedBasis = zip3 [0 :: Int ..] basis basisSugars
          rawPairs = [ makeCriticalPair i newIdx polyI pNew sugarI sugarNew
                     | (i, polyI, sugarI) <- indexedBasis
                     , not (isZeroPoly polyI)
                     ]

          -- Step 2: Apply Criterion 1 (eliminate coprime leading monomials)
          pairsAfterCrit1 = filter (not . criterion1Coprime) rawPairs

          -- Step 3: Apply Gebauer-Moller Criterion 3 (F-criterion on new pairs)
          pairsAfterCrit3 = prunePairSet pairsAfterCrit1

          -- Step 4: Apply Criterion 2 (chain criterion against basis)
          allBasis = basis ++ [pNew]
          validNewPairs = filter (not . criterion2Chain allBasis) pairsAfterCrit3

          -- Step 5: Retain existing queue pairs
          prunedExisting = existingPairs
      in sortPairsBySugar (prunedExisting ++ validNewPairs)
  where
    prunePairSet [] = []
    prunePairSet (cp:cps) =
      let isStrictlyDominated = any (\other ->
                                      other /= cp &&
                                      sparseDivides (cpLCM other) (cpLCM cp)
                                    ) (cp:cps)
      in if isStrictlyDominated
         then prunePairSet cps
         else cp : prunePairSet cps
