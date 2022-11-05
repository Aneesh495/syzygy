{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.Minimization
  ( minimizeComplex
  , isMinimalComplex
  , extractConstant
  ) where

import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (sparseDegree, sparseOne)
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  )
import Syzygy.Polynomial.Arithmetic (subPoly, mulPoly)
import Syzygy.Resolution.FreeResolution (ChainComplex(..))

-- | Extracts the constant term c in k of a polynomial if deg(f) == 0.
-- Returns Nothing if the polynomial is non-constant or zero.
extractConstant :: (Eq c, Num c) => Poly c -> Maybe c
extractConstant (Poly terms _) =
  case terms of
    [Term c m] | sparseDegree m == 0 -> Just c
    _ -> Nothing

-- | Checks whether a chain complex is minimal:
-- A complex is minimal if and only if no differential matrix d_k contains
-- an invertible constant entry (c in k, c /= 0).
isMinimalComplex :: (Eq c, Num c) => ChainComplex c -> Bool
isMinimalComplex cc =
  all (not . hasConstantEntry) (ccDifferentials cc)
  where
    hasConstantEntry mat =
      any (any hasConst) mat

    hasConst p =
      case extractConstant p of
        Just c  -> c /= 0
        Nothing -> False

-- | Minimizes a chain complex of free modules by eliminating trivial summands
-- 0 -> R -> R -> 0 induced by non-zero constant entries in differential matrices d_k.
-- By Nakayama's Lemma, this yields the unique minimal graded free resolution
-- whose module ranks are the invariant graded Betti numbers.
minimizeComplex :: forall c. (Eq c, Fractional c)
                => MonomialOrder
                -> ChainComplex c
                -> ChainComplex c
minimizeComplex polyOrder cc =
  let diffs = ccDifferentials cc
      ranks = ccRanks cc
      shifts = ccShifts cc
  in go 0 diffs ranks shifts
  where
    go !diffIdx currentDiffs currentRanks currentShifts
      | diffIdx >= length currentDiffs =
          ChainComplex currentRanks currentShifts currentDiffs
      | otherwise =
          let mat = currentDiffs !! diffIdx
          in case findScalarPivot mat of
               Nothing ->
                 -- No scalar pivot in d_{diffIdx+1}; advance to next differential
                 go (diffIdx + 1) currentDiffs currentRanks currentShifts
               Just (pivotR, pivotC, pivotVal) ->
                 -- Eliminate pivot row and column in d_{diffIdx+1} and update adjacent matrices
                 let invP = recip pivotVal
                     -- In d_{diffIdx+1}: row pivotR, col pivotC
                     -- Clean other entries in row pivotR: subtract multiple of col pivotC
                     (matCleaned, adjNext) = eliminateColumn pivotR pivotC invP mat
                                               (if diffIdx + 1 < length currentDiffs
                                                then Just (currentDiffs !! (diffIdx + 1))
                                                else Nothing)
                     -- Clean other entries in col pivotC: subtract multiple of row pivotR
                     (matFinal, adjPrev) = eliminateRow pivotR pivotC matCleaned
                                             (if diffIdx > 0
                                              then Just (currentDiffs !! (diffIdx - 1))
                                              else Nothing)
                     -- Delete row pivotR and col pivotC from matFinal
                     matPruned = deleteRowCol pivotR pivotC matFinal
                     -- Update diffs list
                     diffsUpdated = updateList diffIdx matPruned currentDiffs
                     diffsUpdatedPrev = case adjPrev of
                                          Just pMat -> updateList (diffIdx - 1) (deleteCol pivotR pMat) diffsUpdated
                                          Nothing   -> diffsUpdated
                     diffsUpdatedNext = case adjNext of
                                          Just nMat -> updateList (diffIdx + 1) (deleteRow pivotC nMat) diffsUpdatedPrev
                                          Nothing   -> diffsUpdatedPrev
                     -- Update ranks and shifts
                     ranksUpdated = updateRanks diffIdx pivotR pivotC currentRanks
                     shiftsUpdated = updateShifts diffIdx pivotR pivotC currentShifts
                 in go diffIdx diffsUpdatedNext ranksUpdated shiftsUpdated

    findScalarPivot mat =
      let indexed = [ (r, c, val)
                    | (r, row) <- zip [0 :: Int ..] mat
                    , (c, p) <- zip [0 :: Int ..] row
                    , Just val <- [extractConstant p]
                    , val /= 0
                    ]
      in case indexed of
           (p:_) -> Just p
           []    -> Nothing

    -- Eliminates other entries in row pivotR using col pivotC
    eliminateColumn pivotR pivotC invP mat maybeNext =
      let numCols = case mat of (r:_) -> length r; [] -> 0
          pivotCol = [ row !! pivotC | row <- mat ]
          scale = invP
          -- For each other column j /= pivotC:
          -- factor = (mat !! pivotR) !! j * scale
          -- col_j' = col_j - factor * pivotCol
          mat' = [ [ if j == pivotC
                     then (mat !! i) !! j
                     else let factor = mulPoly ((mat !! pivotR) !! j) (Poly [Term scale sparseOne] polyOrder)
                          in subPoly ((mat !! i) !! j) (mulPoly factor (pivotCol !! i))
                   | j <- [0 .. numCols - 1]
                   ]
                 | i <- [0 .. length mat - 1]
                 ]
          next' = maybeNext
      in (mat', next')

    -- Eliminates other entries in col pivotC using row pivotR
    eliminateRow pivotR pivotC mat maybePrev =
      let numRows = length mat
          pivotRow = mat !! pivotR
          mat' = [ if i == pivotR
                   then mat !! i
                   else let factor = (mat !! i) !! pivotC
                        in zipWith (\x y -> subPoly x (mulPoly factor y)) (mat !! i) pivotRow
                 | i <- [0 .. numRows - 1]
                 ]
          prev' = maybePrev
      in (mat', prev')

    deleteRow r mat = [ row | (i, row) <- zip [0 :: Int ..] mat, i /= r ]
    deleteCol c mat = [ [ val | (j, val) <- zip [0 :: Int ..] row, j /= c ] | row <- mat ]
    deleteRowCol r c mat = deleteCol c (deleteRow r mat)

    updateList idx val xs =
      [ if i == idx then val else x | (i, x) <- zip [0 :: Int ..] xs ]

    updateRanks diffIdx _r _c ranks =
      [ if i == diffIdx
        then rank - 1
        else if i == diffIdx + 1
             then rank - 1
             else rank
      | (i, rank) <- zip [0 :: Int ..] ranks
      ]

    updateShifts diffIdx r c shifts =
      [ if i == diffIdx
        then deleteAt r shiftList
        else if i == diffIdx + 1
             then deleteAt c shiftList
             else shiftList
      | (i, shiftList) <- zip [0 :: Int ..] shifts
      ]

    deleteAt idx xs = [ x | (i, x) <- zip [0 :: Int ..] xs, i /= idx ]
