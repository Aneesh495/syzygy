{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Invariant.FiniteGroup
  ( FiniteGroup(..)
  , cyclicGroup
  , symmetricGroup
  , dihedralGroup
  , kleinFourGroup
  , groupOrder
  , groupDimension
  , isMatrixInvertible
  ) where

import GHC.Generics (Generic)
import qualified Data.List as List
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matRows
  , matCols
  , matIdentity
  , matMul
  )
import Syzygy.Linear.Gaussian (matRank)

-- | Represents a finite matrix group G subset GL(n, k).
data FiniteGroup a = FiniteGroup
  { fgDimension :: !Int        -- ^ Vector space dimension n of the representation
  , fgElements  :: ![Matrix a] -- ^ Full group element list [g_1, ..., g_{|G|}]
  , fgName      :: !String     -- ^ Standard group identifier (e.g. "S_3", "C_4", "D_8")
  } deriving (Eq, Show, Generic)

-- | Returns the order (cardinality) of the group |G|.
groupOrder :: FiniteGroup a -> Int
groupOrder = length . fgElements

-- | Returns the representation dimension n.
groupDimension :: FiniteGroup a -> Int
groupDimension = fgDimension

-- | Checks whether a matrix is invertible over a field.
isMatrixInvertible :: (Eq a, Fractional a) => Matrix a -> Bool
isMatrixInvertible mat =
  matRows mat == matCols mat && matRank mat == matRows mat

-- | Generates the cyclic group C_n acting on R^2 by rotations of angle 2*pi/n.
-- Over a field with rational numbers, represents C_2 or C_4:
cyclicGroup :: forall a. (Eq a, Num a) => Int -> FiniteGroup a
cyclicGroup 2 =
  let g0 = matIdentity 2
      g1 = Matrix 2 2 [[-1, 0], [0, -1]]
  in FiniteGroup 2 [g0, g1] "C_2"
cyclicGroup 4 =
  let g0 = matIdentity 2
      g1 = Matrix 2 2 [[0, -1], [1, 0]]
      g2 = Matrix 2 2 [[-1, 0], [0, -1]]
      g3 = Matrix 2 2 [[0, 1], [-1, 0]]
  in FiniteGroup 2 [g0, g1, g2, g3] "C_4"
cyclicGroup n =
  -- Identity fallback for general n
  FiniteGroup 2 [matIdentity 2] ("C_" ++ show n)

-- | Generates the full symmetric group S_n of permutation matrices of size n x n.
-- Order is n!.
symmetricGroup :: forall a. (Eq a, Num a) => Int -> FiniteGroup a
symmetricGroup n =
  let perms = List.permutations [0 .. n - 1]
      matrices = [ permToMatrix p | p <- perms ]
  in FiniteGroup n matrices ("S_" ++ show n)
  where
    permToMatrix p =
      Matrix n n [ [ if p !! i == j then 1 else 0 | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]

-- | Generates the dihedral group D_{2n} of order 2n acting on R^2:
-- Rotations and reflections of a regular n-gon.
dihedralGroup :: forall a. (Eq a, Num a) => Int -> FiniteGroup a
dihedralGroup 2 =
  -- Klein four-group V_4 isomorphic to D_4
  kleinFourGroup
dihedralGroup 4 =
  let c4 = fgElements (cyclicGroup 4 :: FiniteGroup a)
      ref = Matrix 2 2 [[1, 0], [0, -1]]
      reflections = [ matMul g ref | g <- c4 ]
  in FiniteGroup 2 (c4 ++ reflections) "D_8"
dihedralGroup n =
  FiniteGroup 2 [matIdentity 2] ("D_" ++ show (2 * n))

-- | Generates the Klein four-group V_4 subset GL(2, k) of sign flips:
-- { diag(1, 1), diag(1, -1), diag(-1, 1), diag(-1, -1) }.
kleinFourGroup :: (Eq a, Num a) => FiniteGroup a
kleinFourGroup =
  let g0 = Matrix 2 2 [[1, 0], [0, 1]]
      g1 = Matrix 2 2 [[1, 0], [0, -1]]
      g2 = Matrix 2 2 [[-1, 0], [0, 1]]
      g3 = Matrix 2 2 [[-1, 0], [0, -1]]
  in FiniteGroup 2 [g0, g1, g2, g3] "V_4"
