{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.FreeResolution
  ( ChainComplex(..)
  , freeResolution
  , verifyComplex
  , resolutionLength
  , composeDifferentials
  ) where

import GHC.Generics (Generic)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Polynomial.Base
  ( Poly(..)
  , zeroPoly
  , isZeroPoly
  , polyDegree
  )
import Syzygy.Polynomial.Arithmetic (addPoly, mulPoly)
import Syzygy.Resolution.Module
  ( isZeroModPoly
  , modPolyFromComponents
  , modPolyToComponents
  )
import Syzygy.Resolution.ModuleOrder (ModuleOrder(..))
import Syzygy.Resolution.Syzygy (syzygyModule)

-- | Graded chain complex of free modules representing a free resolution:
-- 0 -> F_p -> ... -> F_1 -> F_0 -> M -> 0
-- where each differential d_k : F_k -> F_{k-1} is represented as a matrix of polynomials.
data ChainComplex c = ChainComplex
  { ccRanks         :: ![Int]         -- ^ Ranks [r_0, r_1, ..., r_p] of the free modules
  , ccShifts        :: ![[Int]]       -- ^ Degree shifts for basis elements of each F_k
  , ccDifferentials :: ![[[Poly c]]]  -- ^ Differential matrices [d_1, d_2, ..., d_p]
  } deriving (Eq, Show, Generic)

-- | Returns the projective length p of the free resolution.
-- By Hilbert's Syzygy Theorem, p <= n where n is the number of variables.
resolutionLength :: ChainComplex c -> Int
resolutionLength cc =
  case ccRanks cc of
    [] -> 0
    ranks -> length ranks - 1

-- | Computes a graded free resolution of a polynomial ideal I = <f_1, ..., f_s> in R = k[x_0, ..., x_{n-1}].
-- By Hilbert's Syzygy Theorem (1890), the resolution terminates in at most n steps.
freeResolution :: forall c. (Eq c, Fractional c)
               => MonomialOrder
               -> Int
               -> [Poly c]
               -> ChainComplex c
freeResolution polyOrder numVars polys =
  let nonZeroPolys = filter (not . isZeroPoly) polys
  in if null nonZeroPolys
     then ChainComplex [1] [[0]] []
     else
       let -- F_0 = R^1, F_1 = R^s with generators f_1, ..., f_s
           d1Matrix = [ nonZeroPolys ] -- 1 x s matrix
           f0Shifts = [0]
           f1Shifts = [ polyDegree p | p <- nonZeroPolys ]
           f1Gens = [ modPolyFromComponents [p] | p <- nonZeroPolys ]
           mOrder = TOP polyOrder []
       in resolveNext 1 [1, length nonZeroPolys] [f0Shifts, f1Shifts] [d1Matrix] f1Gens mOrder
  where
    resolveNext !step !ranks !shifts !diffs !currentGens mOrder
      | step > numVars =
          -- Hilbert Syzygy Theorem bound reached
          ChainComplex (reverse ranks) (reverse shifts) (reverse diffs)
      | otherwise =
          let syzs = syzygyModule mOrder numVars currentGens
              nonZeroSyzs = filter (not . isZeroModPoly) syzs
          in if null nonZeroSyzs
             then ChainComplex (reverse ranks) (reverse shifts) (reverse diffs)
             else
               let prevRank = case ranks of (r:_) -> r; [] -> 1
                   nextRank = length nonZeroSyzs
                   -- Convert each syzygy (column of d_{step+1}) into a list of polynomials of length prevRank
                   dNextColumns = [ modPolyToComponents prevRank polyOrder s | s <- nonZeroSyzs ]
                   -- Transpose columns into rows: d_{step+1} is prevRank x nextRank matrix
                   dNextMatrix = [ [ dNextColumns !! colIdx !! rowIdx | colIdx <- [0 .. nextRank - 1] ]
                                 | rowIdx <- [0 .. prevRank - 1]
                                 ]
                   nextShifts = replicate nextRank (step + 1)
                   ranks' = nextRank : ranks
                   shifts' = nextShifts : shifts
                   diffs' = dNextMatrix : diffs
               in resolveNext (step + 1) ranks' shifts' diffs' nonZeroSyzs mOrder

-- | Multiplies two polynomial matrices: C = A * B.
-- A has size (rA x cA), B has size (rB x cB) where cA == rB.
composeDifferentials :: (Eq c, Num c) => MonomialOrder -> [[Poly c]] -> [[Poly c]] -> [[Poly c]]
composeDifferentials polyOrder matA matB =
  let numRowsA = length matA
      numColsA = case matA of (r:_) -> length r; [] -> 0
      numColsB = case matB of (r:_) -> length r; [] -> 0
  in [ [ foldl addPoly (zeroPoly polyOrder)
           [ mulPoly ((matA !! i) !! k) ((matB !! k) !! j) | k <- [0 .. numColsA - 1] ]
       | j <- [0 .. numColsB - 1]
       ]
     | i <- [0 .. numRowsA - 1]
     ]

-- | Verifies the exactness of the chain complex: d_k * d_{k+1} == 0 for all consecutive differentials.
verifyComplex :: (Eq c, Num c) => MonomialOrder -> ChainComplex c -> Bool
verifyComplex polyOrder cc =
  let diffs = ccDifferentials cc
  in case diffs of
       [] -> True
       [_] -> True
       _ -> all isZeroMatrix [ composeDifferentials polyOrder d1 d2 | (d1, d2) <- zip diffs (drop 1 diffs) ]
  where
    isZeroMatrix mat = all (all isZeroPoly) mat
