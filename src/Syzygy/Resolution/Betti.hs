{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.Betti
  ( BettiTable(..)
  , computeBettiTable
  , formatBettiTable
  , totalBettiNumbers
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Resolution.FreeResolution (ChainComplex(..))

-- | Graded Betti table of a module M:
-- Encodes the invariant dimensions beta_{i, j} = dim_k Tor_i^R(M, k)_j
-- indexed by homological degree i (columns) and internal degree j,
-- displayed using the standard Macaulay2 row index (j - i).
data BettiTable = BettiTable
  { bettiTotal      :: ![Int]    -- ^ Total Betti numbers [beta_0, beta_1, ..., beta_p]
  , bettiMatrix     :: ![[Int]]  -- ^ 2D table entries (row = j - i, col = i)
  , bettiMinRow     :: !Int      -- ^ Minimum degree shift offset (min j - i)
  , bettiMaxRow     :: !Int      -- ^ Maximum degree shift offset (max j - i)
  , bettiProjDim    :: !Int      -- ^ Projective dimension pd(M)
  , bettiRegularity :: !Int      -- ^ Castelnuovo-Mumford regularity reg(M)
  } deriving (Eq, Show, Generic)

-- | Computes the total Betti numbers [beta_0, ..., beta_p] directly from a chain complex.
totalBettiNumbers :: ChainComplex c -> [Int]
totalBettiNumbers = ccRanks

-- | Computes the graded Betti table from a graded chain complex.
-- For each free module F_i with generator degree shifts [d_{i, 1}, ..., d_{i, r_i}],
-- the pair (i, d_{i, k}) contributes 1 to beta_{i, d_{i, k}}, which appears at row (d_{i, k} - i).
computeBettiTable :: ChainComplex c -> BettiTable
computeBettiTable cc =
  let shifts = ccShifts cc
      ranks = ccRanks cc
      -- Collect all (col, deg) contributions
      contributions = [ (colIdx, deg)
                      | (colIdx, degList) <- zip [0 :: Int ..] shifts
                      , deg <- degList
                      ]
      -- Map (col, row = deg - col) -> count
      counts = foldl (\acc (col, deg) ->
                       let row = deg - col
                       in Map.insertWith (+) (row, col) 1 acc
                     ) Map.empty contributions
      pDim = if null ranks then 0 else length ranks - 1
      allRows = [ deg - col | (col, deg) <- contributions ]
      minR = if null allRows then 0 else minimum allRows
      maxR = if null allRows then 0 else maximum allRows
      tableRows = [ [ Map.findWithDefault 0 (r, c) counts | c <- [0 .. pDim] ]
                  | r <- [minR .. maxR]
                  ]
      reg = maxR
  in BettiTable
       { bettiTotal      = ranks
       , bettiMatrix     = tableRows
       , bettiMinRow     = minR
       , bettiMaxRow     = maxR
       , bettiProjDim    = pDim
       , bettiRegularity = reg
       }

-- | Formats the Betti table in standard algebraic geometry notation (Macaulay2 / Singular layout).
formatBettiTable :: BettiTable -> String
formatBettiTable bt =
  let totals = bettiTotal bt
      pDim = bettiProjDim bt
      colWidth = 5
      padLeft n s = replicate (max 0 (n - length s)) ' ' ++ s

      header = "       total: " ++ concatMap (padLeft colWidth . show) totals
      sep = "       " ++ replicate (7 + (pDim + 1) * colWidth) '-'

      formatRow (rowIdx, entries) =
        let rowLabel = padLeft 5 (show rowIdx) ++ ": "
            rowEntries = concatMap (\val ->
                                      padLeft colWidth (if val == 0 then "-" else show val)
                                   ) entries
        in rowLabel ++ rowEntries

      rows = zip [bettiMinRow bt .. bettiMaxRow bt] (bettiMatrix bt)
      rowStrs = map formatRow rows
  in unlines (header : sep : rowStrs)
