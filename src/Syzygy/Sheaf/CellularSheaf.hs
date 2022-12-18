{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Sheaf.CellularSheaf
Description : Cellular sheaves, stalks, restriction maps, and local systems
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Cellular sheaves on regular cell complexes represent data over topological spaces
with linear constraints between incident cells (Curry 2014).

A cellular sheaf F on a cell complex X assigns:
1. To each cell sigma in X, a stalk vector space F(sigma) of dimension d_sigma.
2. To each incidence relation sigma <| tau (where sigma is a facet of tau),
   a linear restriction map:
     rho_{sigma <| tau} : F(sigma) -> F(tau)
   represented as a matrix of size d_tau x d_sigma.

Functoriality requires that for any codimension-2 incidence sigma <| gamma with
intermediate facets tau_1, tau_2, the diagram commutes:
  rho_{tau_1 <| gamma} . rho_{sigma <| tau_1} = rho_{tau_2 <| gamma} . rho_{sigma <| tau_2}
-}
module Syzygy.Sheaf.CellularSheaf
  ( CellularSheaf(..)
  , emptySheaf
  , stalkDimension
  , restrictionMap
  , constantSheaf
  , constantSheaf1D
  , circleMonodromySheaf
  , verifySheafCommutativity
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matIdentity
  , matZero
  , matMul
  )
import Syzygy.Sheaf.Complex
  ( Cell(..)
  , CellComplex(..)
  , cellsOfDim
  )

-- | A cellular sheaf over a cell complex X.
data CellularSheaf a = CellularSheaf
  { csComplex      :: !CellComplex                       -- ^ Underlying cell complex
  , csStalkDims    :: !(Map.Map Int Int)                 -- ^ Map from cellId to stalk dimension
  , csRestrictions :: !(Map.Map (Int, Int) (Matrix a))   -- ^ Map from (facetId, cellId) to restriction matrix
  } deriving (Eq, Show, Generic)

-- | Constructs an empty sheaf on an empty complex.
emptySheaf :: CellComplex -> CellularSheaf a
emptySheaf cc = CellularSheaf cc Map.empty Map.empty

-- | Retrieves the stalk dimension d_sigma of a cell. Defaults to 0 if unspecified.
stalkDimension :: CellularSheaf a -> Int -> Int
stalkDimension sh cId =
  Map.findWithDefault 0 cId (csStalkDims sh)

-- | Retrieves the restriction map rho_{sigma <| tau} from cell sigma to cell tau.
-- Returns a zero matrix of appropriate size if no map is specified.
restrictionMap :: Num a => CellularSheaf a -> Int -> Int -> Matrix a
restrictionMap sh facetId cId =
  case Map.lookup (facetId, cId) (csRestrictions sh) of
    Just m  -> m
    Nothing ->
      let dFacet = stalkDimension sh facetId
          dCell  = stalkDimension sh cId
      in matZero dCell dFacet

-- | Constructs the constant sheaf F = Q^d on a cell complex X.
-- Stalk dimension is d on every cell, and all restriction maps are identity matrices.
constantSheaf :: Num a => CellComplex -> Int -> CellularSheaf a
constantSheaf cc dimVal =
  let stalkMap = Map.fromList [ (cellId c, dimVal) | c <- ccCells cc ]
      ident = matIdentity dimVal
      restrMap = Map.fromList [ ((fId, cId), ident) | (fId, cId) <- Map.keys (ccIncidences cc) ]
  in CellularSheaf cc stalkMap restrMap

-- | Constructs the standard rank-1 constant sheaf F = Q.
constantSheaf1D :: Num a => CellComplex -> CellularSheaf a
constantSheaf1D cc = constantSheaf cc 1

-- | Constructs a rank-2 local system (flat bundle) over the circle S^1
-- with holonomy rotation angle theta in radians.
-- For theta = 0, this is the trivial bundle; for theta = pi, this is the Mobius bundle.
circleMonodromySheaf :: (Floating a, Eq a) => CellComplex -> a -> CellularSheaf a
circleMonodromySheaf cc theta =
  let stalkMap = Map.fromList [ (cellId c, 2) | c <- ccCells cc ]
      rotMat = Matrix 2 2 [[cos theta, -sin theta], [sin theta, cos theta]]
      ident = matIdentity 2
      -- Apply rotation to the final edge, identity to all other edges
      incKeys = Map.keys (ccIncidences cc)
      makeRestr (fId, cId) =
        if cId == maximum (map cellId (cellsOfDim cc 1)) && fId == 0
        then ((fId, cId), rotMat)
        else ((fId, cId), ident)
      restrMap = Map.fromList (map makeRestr incKeys)
  in CellularSheaf cc stalkMap restrMap

-- | Verifies that all restriction maps satisfy the commutativity relation
-- across all codimension-2 inclusions.
verifySheafCommutativity :: (Eq a, Num a) => CellularSheaf a -> Int -> Bool
verifySheafCommutativity sh k =
  let cc = csComplex sh
      gammaCells = cellsOfDim cc (k + 1)
      sigmaCells = cellsOfDim cc (k - 1)
      midCells   = cellsOfDim cc k

      checkPair sigma gamma =
        let sId = cellId sigma
            gId = cellId gamma
            -- Find all intermediate k-cells tau incident to both sigma and gamma
            intermediates =
              [ cellId tau
              | tau <- midCells
              , let tId = cellId tau
              , Map.member (sId, tId) (ccIncidences cc)
              , Map.member (tId, gId) (ccIncidences cc)
              ]
        in case intermediates of
             (t1 : t2 : _) ->
               let rho1 = matMul (restrictionMap sh t1 gId) (restrictionMap sh sId t1)
                   rho2 = matMul (restrictionMap sh t2 gId) (restrictionMap sh sId t2)
               in rho1 == rho2
             _ -> True
  in all (\g -> all (\s -> checkPair s g) sigmaCells) gammaCells
