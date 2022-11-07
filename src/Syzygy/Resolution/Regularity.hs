{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.Regularity
  ( projectiveDimension
  , castelnuovoMumfordRegularity
  , moduleDepth
  , isCohenMacaulay
  , isGorenstein
  , isLinearResolution
  , isPureResolution
  , HomologicalInvariants(..)
  , computeHomologicalInvariants
  ) where

import GHC.Generics (Generic)
import Syzygy.Resolution.FreeResolution (ChainComplex(..))
import Syzygy.Resolution.Betti
  ( BettiTable(..)
  , computeBettiTable
  )

-- | Summary of homological invariants extracted from the minimal free resolution.
data HomologicalInvariants = HomologicalInvariants
  { invProjDim       :: !Int   -- ^ Projective dimension pd_R(M)
  , invDepth         :: !Int   -- ^ Depth via Auslander-Buchsbaum formula: depth(M) = n - pd(M)
  , invRegularity    :: !Int   -- ^ Castelnuovo-Mumford regularity reg(M)
  , invIsCohenMacaulay :: !Bool -- ^ Whether M is Cohen-Macaulay
  , invIsGorenstein  :: !Bool  -- ^ Whether M is Gorenstein (Cohen-Macaulay with beta_p = 1)
  , invIsLinear      :: !Bool  -- ^ Whether the resolution is linear (differentials linear)
  , invIsPure        :: !Bool  -- ^ Whether the resolution is pure (each F_i has single degree shift)
  } deriving (Eq, Show, Generic)

-- | Returns the projective dimension pd_R(M) = max { i | F_i /= 0 in minimal resolution }.
projectiveDimension :: BettiTable -> Int
projectiveDimension = bettiProjDim

-- | Returns the Castelnuovo-Mumford regularity reg(M) = max { j - i | beta_{i, j} /= 0 }.
castelnuovoMumfordRegularity :: BettiTable -> Int
castelnuovoMumfordRegularity = bettiRegularity

-- | Computes depth(M) via the Auslander-Buchsbaum theorem:
-- pd_R(M) + depth(M) = depth(R) = n.
moduleDepth :: Int -> BettiTable -> Int
moduleDepth numVars bt = max 0 (numVars - bettiProjDim bt)

-- | Checks whether M is Cohen-Macaulay given its Krull dimension:
-- M is Cohen-Macaulay iff depth(M) == dim(M).
isCohenMacaulay :: Int -> Int -> BettiTable -> Bool
isCohenMacaulay numVars krullDim bt =
  moduleDepth numVars bt == krullDim

-- | Checks whether M is Gorenstein:
-- Cohen-Macaulay with top Betti number beta_p == 1.
isGorenstein :: Int -> Int -> BettiTable -> Bool
isGorenstein numVars krullDim bt =
  isCohenMacaulay numVars krullDim bt &&
  case reverse (filter (/= 0) (bettiTotal bt)) of
    (topBetti:_) -> topBetti == 1
    []           -> False

-- | Checks whether the resolution is pure:
-- each free module F_i is generated in a single degree (all shifts in F_i are identical).
isPureResolution :: ChainComplex c -> Bool
isPureResolution cc =
  all singleDegree (ccShifts cc)
  where
    singleDegree [] = True
    singleDegree (d:ds) = all (== d) ds

-- | Checks whether the resolution is linear:
-- pure resolution where shifts increase by exactly 1 at each homological step (deg(d_k) = 1).
isLinearResolution :: ChainComplex c -> Bool
isLinearResolution cc =
  isPureResolution cc &&
  let shifts = [ case sList of (d:_) -> d; [] -> 0 | sList <- ccShifts cc ]
  in case shifts of
       [] -> True
       [_] -> True
       (firstShift:restShifts) ->
         and [ nextShift == prevShift + 1 | (prevShift, nextShift) <- zip (firstShift : restShifts) restShifts ]

-- | Computes all primary homological invariants from a minimal chain complex and ring dimension.
computeHomologicalInvariants :: Int -> Int -> ChainComplex c -> HomologicalInvariants
computeHomologicalInvariants numVars krullDim cc =
  let bt = computeBettiTable cc
      pDim = bettiProjDim bt
      reg = bettiRegularity bt
      depth = moduleDepth numVars bt
      cm = isCohenMacaulay numVars krullDim bt
      gor = isGorenstein numVars krullDim bt
      pureRes = isPureResolution cc
      linRes = isLinearResolution cc
  in HomologicalInvariants
       { invProjDim         = pDim
       , invDepth           = depth
       , invRegularity      = reg
       , invIsCohenMacaulay = cm
       , invIsGorenstein    = gor
       , invIsLinear        = linRes
       , invIsPure          = pureRes
       }
