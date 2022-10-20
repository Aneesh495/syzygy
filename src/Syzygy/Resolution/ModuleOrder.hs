{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.ModuleOrder
  ( ModuleOrder(..)
  , cmpModMonomial
  , sortModPoly
  , canonicalizeModPoly
  ) where

import GHC.Generics (Generic)
import qualified Data.List as List
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( sparseDegree
  , sparseMul
  , compareSparse
  )
import Syzygy.Resolution.Module
  ( ModMonomial(..)
  , ModTerm(..)
  , ModPoly(..)
  )

-- | Monomial ordering on a graded free module R^m.
data ModuleOrder
  = TOP !MonomialOrder ![Int]
    -- ^ Term Over Position: compares total shifted degree first (if shifts provided),
    -- then monomial order, and breaks ties by component index (smaller index first).
  | POT !MonomialOrder
    -- ^ Position Over Term: compares component index first (smaller index first),
    -- then breaks ties by polynomial monomial order. Essential for module elimination.
  | SchreyerOrder !ModuleOrder ![ModMonomial]
    -- ^ Schreyer order induced by generators g_0, ..., g_{m-1} with leading monomials LM(g_i).
    -- Compares m_1 * LM(g_i) with m_2 * LM(g_j) using the underlying module order.
    -- By Schreyer's Theorem, syzygies obtained from S-reductions form a Groebner basis
    -- automatically under this order.
  deriving (Eq, Show, Generic)

-- | Compares two module monomials under a given ModuleOrder in numVars variables.
cmpModMonomial :: ModuleOrder -> Int -> ModMonomial -> ModMonomial -> Ordering
cmpModMonomial order numVars (ModMonomial exp1 c1) (ModMonomial exp2 c2) =
  case order of
    TOP polyOrder shifts ->
      let shift1 = if c1 < length shifts then shifts !! c1 else 0
          shift2 = if c2 < length shifts then shifts !! c2 else 0
          deg1 = sparseDegree exp1 + shift1
          deg2 = sparseDegree exp2 + shift2
      in case compare deg1 deg2 of
           LT -> LT
           GT -> GT
           EQ ->
             case compareSparse polyOrder numVars exp1 exp2 of
               LT -> LT
               GT -> GT
               EQ -> compare c2 c1 -- smaller component index is greater: e_0 > e_1 > ...

    POT polyOrder ->
      case compare c2 c1 of -- smaller component index is greater
        LT -> LT
        GT -> GT
        EQ -> compareSparse polyOrder numVars exp1 exp2

    SchreyerOrder baseOrder leadingGens ->
      let lead1 = if c1 < length leadingGens then leadingGens !! c1 else ModMonomial exp1 c1
          lead2 = if c2 < length leadingGens then leadingGens !! c2 else ModMonomial exp2 c2
          shifted1 = ModMonomial (sparseMul exp1 (modMonoExp lead1)) (modMonoComp lead1)
          shifted2 = ModMonomial (sparseMul exp2 (modMonoExp lead2)) (modMonoComp lead2)
      in case cmpModMonomial baseOrder numVars shifted1 shifted2 of
           LT -> LT
           GT -> GT
           EQ -> compare c2 c1 -- tie break on generator index

-- | Sorts the terms of a module polynomial in strictly descending order with respect to
-- the specified ModuleOrder.
sortModPoly :: ModuleOrder -> Int -> ModPoly c -> ModPoly c
sortModPoly order numVars (ModPoly terms) =
  let sortedTerms = List.sortBy (\(ModTerm _ m1 c1) (ModTerm _ m2 c2) ->
                                  cmpModMonomial order numVars (ModMonomial m2 c2) (ModMonomial m1 c1)
                                ) terms
  in ModPoly sortedTerms

-- | Combines like terms and sorts in strictly descending order.
canonicalizeModPoly :: (Eq c, Num c) => ModuleOrder -> Int -> ModPoly c -> ModPoly c
canonicalizeModPoly order numVars (ModPoly terms) =
  let combined = foldl (\acc (ModTerm c m compIdx) ->
                         Map.insertWith (+) (m, compIdx) c acc
                       ) Map.empty terms
      nonZeroTerms = [ ModTerm c m compIdx
                     | ((m, compIdx), c) <- Map.toList combined
                     , c /= 0
                     ]
      sorted = List.sortBy (\(ModTerm _ m1 c1) (ModTerm _ m2 c2) ->
                             cmpModMonomial order numVars (ModMonomial m2 c2) (ModMonomial m1 c1)
                           ) nonZeroTerms
  in ModPoly sorted
