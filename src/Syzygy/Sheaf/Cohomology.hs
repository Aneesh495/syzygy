{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Sheaf.Cohomology
Description : Sheaf cochains, coboundary operator delta, and sheaf cohomology H^k(X; F)
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Cellular sheaf cohomology computes the global obstructions and consistent
assignments of local data constrained by restriction maps (Curry 2014).

The space of k-cochains is the direct sum of stalks over all k-cells:
  C^k(X; F) = direct_sum_{sigma in X_k} F(sigma)

The coboundary operator delta^k : C^k(X; F) -> C^{k+1}(X; F) is defined on
a cochain x in C^k by its evaluation on a (k+1)-cell tau:
  (delta^k x)_tau = sum_{sigma <| tau} [sigma : tau] * rho_{sigma <| tau}(x_sigma)

The k-th sheaf cohomology vector space is:
  H^k(X; F) = ker(delta^k) / im(delta^{k-1})

with Betti number b_k(X; F) = dim H^k(X; F) = dim C^k - rank(delta^k) - rank(delta^{k-1}).
For the rank-1 constant sheaf F = Q, this recovers classical cellular cohomology.
-}
module Syzygy.Sheaf.Cohomology
  ( cochainDimension
  , coboundaryMatrix
  , verifyCoboundarySquareZero
  , sheafBettiNumber
  , allSheafBettiNumbers
  , sheafGlobalSections
  , sheafEulerCharacteristic
  ) where

import qualified Data.Map.Strict as Map
import Syzygy.Linear.Vector (Vector(..))
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matMul
  , matZero
  )
import Syzygy.Linear.Gaussian
  ( matRank
  , nullSpace
  )
import Syzygy.Sheaf.Complex
  ( Cell(..)
  , CellComplex(..)
  , complexDimension
  , cellsOfDim
  )
import Syzygy.Sheaf.CellularSheaf
  ( CellularSheaf(..)
  , stalkDimension
  , restrictionMap
  )

-- | Computes the total dimension of the k-cochain vector space:
-- dim C^k(X; F) = sum_{sigma in X_k} dim F(sigma).
cochainDimension :: CellularSheaf a -> Int -> Int
cochainDimension sh k =
  let cells = cellsOfDim (csComplex sh) k
  in sum [ stalkDimension sh (cellId c) | c <- cells ]

-- | Computes the sheaf coboundary matrix delta^k : C^k(X; F) -> C^{k+1}(X; F).
-- The block at row tau in X_{k+1} and column sigma in X_k is:
--   [sigma : tau] * rho_{sigma <| tau}
coboundaryMatrix :: forall a. (Eq a, Fractional a) => CellularSheaf a -> Int -> Matrix a
coboundaryMatrix sh k =
  let cc = csComplex sh
      kCells = cellsOfDim cc k
      kp1Cells = cellsOfDim cc (k + 1)
      dimCk = cochainDimension sh k
      dimCkp1 = cochainDimension sh (k + 1)

      -- Compute column offset for each k-cell
      colOffsets = computeOffsets [ (cellId c, stalkDimension sh (cellId c)) | c <- kCells ]
      -- Compute row offset for each (k+1)-cell
      rowOffsets = computeOffsets [ (cellId c, stalkDimension sh (cellId c)) | c <- kp1Cells ]

      -- Populate matrix block by block
      fillEntries =
        [ let signVal = fromIntegral (Map.findWithDefault 0 (sId, tId) (ccIncidences cc))
              Matrix rDims cDims rEntries = restrictionMap sh sId tId
              rBase = rowOffsets Map.! tId
              cBase = colOffsets Map.! sId
          in [ ((rBase + r, cBase + c), signVal * ((rEntries !! r) !! c))
             | r <- [0 .. rDims - 1]
             , c <- [0 .. cDims - 1]
             ]
        | tau <- kp1Cells
        , let tId = cellId tau
        , sigma <- kCells
        , let sId = cellId sigma
        , Map.member (sId, tId) (ccIncidences cc)
        ]
      entryMap = Map.fromList (concat fillEntries)
      fullGrid =
        [ [ Map.findWithDefault 0 (r, c) entryMap | c <- [0 .. dimCk - 1] ]
        | r <- [0 .. dimCkp1 - 1]
        ]
  in Matrix dimCkp1 dimCk fullGrid
  where
    computeOffsets pairs =
      let go _ [] = []
          go !acc ((cId, d) : rest) = (cId, acc) : go (acc + d) rest
      in Map.fromList (go 0 pairs)

-- | Verifies that delta^{k+1} . delta^k = 0 as matrices.
verifyCoboundarySquareZero :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Bool
verifyCoboundarySquareZero sh k =
  let dk = coboundaryMatrix sh k
      dkp1 = coboundaryMatrix sh (k + 1)
      numR = case dkp1 of Matrix r _ _ -> r
      numC = case dk of Matrix _ c _ -> c
  in if numR == 0 || numC == 0
     then True
     else
       let prod = matMul dkp1 dk
           zeroM = matZero numR numC
       in prod == zeroM

-- | Computes the k-th sheaf Betti number b_k(X; F) = dim H^k(X; F).
-- Uses rank-nullity: b_k = dim C^k - rank(delta^k) - rank(delta^{k-1}).
sheafBettiNumber :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Int
sheafBettiNumber sh k =
  let dimCk = cochainDimension sh k
  in if dimCk == 0
     then 0
     else
       let rankDk = if k >= complexDimension (csComplex sh)
                    then 0
                    else matRank (coboundaryMatrix sh k)
           rankDkm1 = if k <= 0
                      then 0
                      else matRank (coboundaryMatrix sh (k - 1))
       in max 0 (dimCk - rankDk - rankDkm1)

-- | Computes all sheaf Betti numbers [b_0, b_1, ..., b_d] up to the dimension of X.
allSheafBettiNumbers :: (Eq a, Fractional a) => CellularSheaf a -> [Int]
allSheafBettiNumbers sh =
  let maxDim = complexDimension (csComplex sh)
  in [ sheafBettiNumber sh k | k <- [0 .. maxDim] ]

-- | Computes a basis for the space of global sections H^0(X; F) = ker(delta^0).
sheafGlobalSections :: (Eq a, Fractional a) => CellularSheaf a -> [Vector a]
sheafGlobalSections sh =
  let d0 = coboundaryMatrix sh 0
  in nullSpace d0

-- | Computes the Euler-Poincare characteristic of the sheaf:
-- chi(X; F) = sum_{k=0}^d (-1)^k * dim C^k(X; F) = sum_{k=0}^d (-1)^k * b_k(X; F).
sheafEulerCharacteristic :: CellularSheaf a -> Int
sheafEulerCharacteristic sh =
  let maxDim = complexDimension (csComplex sh)
      terms = [ (if even k then 1 else -1) * cochainDimension sh k | k <- [0 .. maxDim] ]
  in sum terms
