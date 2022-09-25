{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Linear.Vector
  ( Vector(..)
  , SparseVector(..)
  , vecFromList
  , vecToList
  , vecDim
  , vecIndex
  , vecZero
  , vecAdd
  , vecSub
  , vecScale
  , vecDot
  , svecFromList
  , svecToList
  , svecDim
  , svecIndex
  , svecZero
  , svecAdd
  , svecScale
  , svecDot
  , svecToDense
  , denseToSparse
  ) where

import qualified Data.Map.Strict as Map
import GHC.Generics (Generic)

-- | Dense vector of elements with fixed dimension.
data Vector a = Vector
  { vDim      :: !Int
  , vElements :: ![a]
  } deriving (Eq, Generic)

-- | Constructs a dense vector from a list of elements.
vecFromList :: [a] -> Vector a
vecFromList xs = Vector (length xs) xs

-- | Converts a dense vector to a list.
vecToList :: Vector a -> [a]
vecToList (Vector _ xs) = xs

-- | Returns the dimension of the dense vector.
vecDim :: Vector a -> Int
vecDim (Vector d _) = d

-- | 0-indexed element lookup in dense vector.
vecIndex :: Vector a -> Int -> a
vecIndex (Vector d xs) i
  | i < 0 || i >= d = error "vecIndex: index out of bounds"
  | otherwise       = xs !! i

-- | Zero dense vector of dimension d.
vecZero :: Num a => Int -> Vector a
vecZero d = Vector d (replicate d 0)

-- | Vector addition: u + v.
vecAdd :: Num a => Vector a -> Vector a -> Vector a
vecAdd (Vector d1 xs) (Vector d2 ys)
  | d1 /= d2  = error "vecAdd: dimensions do not match"
  | otherwise = Vector d1 (zipWith (+) xs ys)

-- | Vector subtraction: u - v.
vecSub :: Num a => Vector a -> Vector a -> Vector a
vecSub (Vector d1 xs) (Vector d2 ys)
  | d1 /= d2  = error "vecSub: dimensions do not match"
  | otherwise = Vector d1 (zipWith (-) xs ys)

-- | Scalar multiplication: c * v.
vecScale :: Num a => a -> Vector a -> Vector a
vecScale c (Vector d xs) = Vector d (map (c *) xs)

-- | Dot product: u . v.
vecDot :: Num a => Vector a -> Vector a -> a
vecDot (Vector d1 xs) (Vector d2 ys)
  | d1 /= d2  = error "vecDot: dimensions do not match"
  | otherwise = sum (zipWith (*) xs ys)

-- | Sparse vector of fixed dimension storing only non-zero entries in a Map.
data SparseVector a = SparseVector
  { svDim     :: !Int
  , svEntries :: !(Map.Map Int a)
  } deriving (Eq, Generic)

-- | Constructs a sparse vector from an association list [(index, value)].
svecFromList :: (Eq a, Num a) => Int -> [(Int, a)] -> SparseVector a
svecFromList d entries =
  let nonZero = filter (\(idx, val) -> idx >= 0 && idx < d && val /= 0) entries
  in SparseVector d (Map.fromList nonZero)

-- | Converts a sparse vector to an association list of non-zero entries.
svecToList :: SparseVector a -> [(Int, a)]
svecToList (SparseVector _ m) = Map.toList m

-- | Returns the dimension of the sparse vector.
svecDim :: SparseVector a -> Int
svecDim (SparseVector d _) = d

-- | 0-indexed element lookup in sparse vector. Returns 0 for missing indices.
svecIndex :: Num a => SparseVector a -> Int -> a
svecIndex (SparseVector d m) i
  | i < 0 || i >= d = error "svecIndex: index out of bounds"
  | otherwise       = Map.findWithDefault 0 i m

-- | Zero sparse vector of dimension d.
svecZero :: Int -> SparseVector a
svecZero d = SparseVector d Map.empty

-- | Sparse vector addition.
svecAdd :: (Eq a, Num a) => SparseVector a -> SparseVector a -> SparseVector a
svecAdd (SparseVector d1 m1) (SparseVector d2 m2)
  | d1 /= d2  = error "svecAdd: dimensions do not match"
  | otherwise =
      let sumMap = Map.filter (/= 0) (Map.unionWith (+) m1 m2)
      in SparseVector d1 sumMap

-- | Sparse vector scalar scaling.
svecScale :: (Eq a, Num a) => a -> SparseVector a -> SparseVector a
svecScale 0 (SparseVector d _) = SparseVector d Map.empty
svecScale c (SparseVector d m) =
  let scaled = Map.map (c *) m
  in SparseVector d scaled

-- | Sparse vector dot product.
svecDot :: Num a => SparseVector a -> SparseVector a -> a
svecDot (SparseVector d1 m1) (SparseVector d2 m2)
  | d1 /= d2  = error "svecDot: dimensions do not match"
  | otherwise =
      let isect = Map.intersectionWith (*) m1 m2
      in Map.foldl' (+) 0 isect

-- | Converts a sparse vector to a dense vector.
svecToDense :: Num a => SparseVector a -> Vector a
svecToDense sv@(SparseVector d _) =
  Vector d [ svecIndex sv i | i <- [0 .. d - 1] ]

-- | Converts a dense vector to a sparse vector.
denseToSparse :: (Eq a, Num a) => Vector a -> SparseVector a
denseToSparse (Vector d xs) =
  svecFromList d (zip [0 .. d - 1] xs)

instance Show a => Show (Vector a) where
  show (Vector _ xs) = "[" ++ unwords (map show xs) ++ "]"

instance Show a => Show (SparseVector a) where
  show (SparseVector d m) =
    "SparseVector(dim=" ++ show d ++ ", " ++ show (Map.toList m) ++ ")"
