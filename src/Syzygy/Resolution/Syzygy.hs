{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.Syzygy
  ( syzygyModule
  , idealSyzygies
  , verifySyzygy
  , SyzygyResult(..)
  ) where

import GHC.Generics (Generic)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( sparseLCM
  , sparseDiv
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  )
import Syzygy.Resolution.Module
  ( ModMonomial(..)
  , ModTerm(..)
  , ModPoly(..)
  , zeroModPoly
  , isZeroModPoly
  , addModPoly
  , subModPoly
  , modLeadingTerm
  , modLeadingMonomial
  , modPolyFromComponents
  , mulPolyModPoly
  )
import Syzygy.Resolution.ModuleOrder
  ( ModuleOrder(..)
  , canonicalizeModPoly
  )
import Syzygy.Resolution.ModuleGrobner
  ( moduleBuchberger
  , moduleDivision
  , moduleReduce
  )

-- | Detailed result of a syzygy computation.
data SyzygyResult c = SyzygyResult
  { syzGenerators :: ![ModPoly c]
    -- ^ The computed syzygy generators in R^s
  , syzOriginalBasis :: ![ModPoly c]
    -- ^ The Groebner basis used to generate the syzygies
  } deriving (Eq, Show, Generic)

-- | Computes the syzygy module ker(phi) for a submodule M subset R^m given by generators [g_1, ..., g_s].
-- Uses Schreyer's Theorem:
-- 1. Computes a Groebner basis G = [g_1, ..., g_k] of the submodule.
-- 2. For each critical pair (i, j) with i < j sharing a component:
--    S(g_i, g_j) = u_ij * g_i - u_ji * g_j = sum_k q_k * g_k
-- 3. The syzygy vector sigma_ij = u_ij * e_i - u_ji * e_j - sum_k q_k * e_k in R^k satisfies
--    phi(sigma_ij) = 0.
-- 4. Reduces the syzygies to obtain a minimal generating set.
syzygyModule :: forall c. (Eq c, Fractional c)
             => ModuleOrder
             -> Int
             -> [ModPoly c]
             -> [ModPoly c]
syzygyModule order numVars gens =
  let gb = moduleReduce order numVars (moduleBuchberger order numVars gens)
      sRank = length gb
      pOrder = extractPolyOrder order
      pairs = [ (i, j)
              | i <- [0 .. sRank - 1]
              , j <- [i + 1 .. sRank - 1]
              , shareComponent (gb !! i) (gb !! j)
              ]
      syzs = [ computeSchreyerSyzygy i j (gb !! i) (gb !! j) gb pOrder
             | (i, j) <- pairs
             ]
      nonZeroSyzs = filter (not . isZeroModPoly) syzs
      schreyerOrder = TOP pOrder []
  in moduleReduce schreyerOrder numVars nonZeroSyzs
  where
    extractPolyOrder (TOP po _) = po
    extractPolyOrder (POT po) = po
    extractPolyOrder (SchreyerOrder bo _) = extractPolyOrder bo

    shareComponent g1 g2 =
      case (modLeadingMonomial g1, modLeadingMonomial g2) of
        (Just (ModMonomial _ c1), Just (ModMonomial _ c2)) -> c1 == c2
        _ -> False

    computeSchreyerSyzygy i j gI gJ gb pOrder =
      case (modLeadingTerm gI, modLeadingTerm gJ) of
        (Just (ModTerm cI mI _), Just (ModTerm cJ mJ _)) ->
          let lcmM = sparseLCM mI mJ
              uIExp = sparseDiv lcmM mI
              uJExp = sparseDiv lcmM mJ
              uIPoly = Poly [Term (1 / cI) uIExp] pOrder
              uJPoly = Poly [Term (1 / cJ) uJExp] pOrder
              sPoly = canonicalizeModPoly order numVars
                        (subModPoly (mulPolyModPoly uIPoly gI) (mulPolyModPoly uJPoly gJ))
              (quotients, _) = moduleDivision order numVars sPoly gb
              -- Construct syzygy: uI * e_i - uJ * e_j - sum_k (q_k * e_k)
              uITerm = ModTerm (1 / cI) uIExp i
              uJTerm = ModTerm (-(1 / cJ)) uJExp j
              qTerms = [ ModTerm (-c) m k
                       | (k, Poly qP _) <- zip [0 :: Int ..] quotients
                       , Term c m <- qP
                       ]
              allTerms = uITerm : uJTerm : qTerms
          in canonicalizeModPoly (TOP pOrder []) numVars (ModPoly allTerms)
        _ -> zeroModPoly

-- | Computes the syzygy module of an ideal I = <f_1, ..., f_r> in the polynomial ring R.
-- Each polynomial f_i is treated as an element f_i * e_0 in the rank 1 free module R^1.
idealSyzygies :: forall c. (Eq c, Fractional c)
              => MonomialOrder
              -> Int
              -> [Poly c]
              -> [ModPoly c]
idealSyzygies polyOrder numVars polys =
  let nonZeroPolys = filter (not . isZeroPoly) polys
      modGens = [ modPolyFromComponents [p] | p <- nonZeroPolys ]
      mOrder = TOP polyOrder []
  in syzygyModule mOrder numVars modGens

-- | Verifies whether a module vector sigma = (h_1, ..., h_s) in R^s is a valid syzygy
-- of the generator list [g_1, ..., g_s], i.e., sum_{i=1}^s h_i * g_i = 0.
verifySyzygy :: (Eq c, Fractional c)
             => [ModPoly c]
             -> ModPoly c
             -> Bool
verifySyzygy gens (ModPoly sigmaTerms) =
  let groupedTerms = foldl (\acc (ModTerm c m compIdx) ->
                             if compIdx < length gens
                             then let g = gens !! compIdx
                                      pOrder = case modLeadingTerm g of
                                                 _ -> DegRevLex
                                      polyFactor = Poly [Term c m] pOrder
                                      termProduct = mulPolyModPoly polyFactor g
                                  in addModPoly acc termProduct
                             else acc
                           ) zeroModPoly sigmaTerms
  in isZeroModPoly groupedTerms
