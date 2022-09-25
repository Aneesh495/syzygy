{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Linear.Matrix
  ( Matrix(..)
  , matFromRows
  , matFromCols
  , matZero
  , matIdentity
  , matDiag
  , matRows
  , matCols
  , matGet
  , matRow
  , matCol
  , matToList
  , matAdd
  , matSub
  , matScale
  , matMul
  , matVecMul
  , matTranspose
  , matTrace
  , subMatrix
  , minorMatrix
  , swapRows
  , scaleRow
  , addRowMultiple
  ) where

import Data.List (transpose)
import GHC.Generics (Generic)
import Syzygy.Linear.Vector (Vector(..), vecFromList, vecToList)

-- | Exact m x n matrix represented as row-major lists of elements.
data Matrix a = Matrix
  { mNumRows :: !Int
  , mNumCols :: !Int
  , mData    :: ![[a]]
  } deriving (Eq, Generic)

-- | Constructs a matrix from a list of rows.
matFromRows :: [[a]] -> Matrix a
matFromRows [] = Matrix 0 0 []
matFromRows rows@(r0:_) =
  let numR = length rows
      numC = length r0
  in if all (\r -> length r == numC) rows
     then Matrix numR numC rows
     else error "matFromRows: inconsistent row lengths"

-- | Constructs a matrix from a list of columns.
matFromCols :: [[a]] -> Matrix a
matFromCols cols = matTranspose (matFromRows cols)

-- | Zero matrix of size m x n.
matZero :: Num a => Int -> Int -> Matrix a
matZero m n = Matrix m n (replicate m (replicate n 0))

-- | n x n Identity matrix.
matIdentity :: Num a => Int -> Matrix a
matIdentity n =
  Matrix n n [ [ if i == j then 1 else 0 | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]

-- | Diagonal matrix with specified diagonal entries.
matDiag :: Num a => [a] -> Matrix a
matDiag ds =
  let n = length ds
  in Matrix n n [ [ if i == j then ds !! i else 0 | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]

-- | Number of rows.
matRows :: Matrix a -> Int
matRows = mNumRows

-- | Number of columns.
matCols :: Matrix a -> Int
matCols = mNumCols

-- | 0-indexed element lookup at (row i, col j).
matGet :: Matrix a -> Int -> Int -> a
matGet (Matrix numR numC rows) i j
  | i < 0 || i >= numR = error "matGet: row index out of bounds"
  | j < 0 || j >= numC = error "matGet: col index out of bounds"
  | otherwise          = (rows !! i) !! j

-- | Extracts row i as a dense vector.
matRow :: Matrix a -> Int -> Vector a
matRow (Matrix numR _ rows) i
  | i < 0 || i >= numR = error "matRow: index out of bounds"
  | otherwise          = vecFromList (rows !! i)

-- | Extracts column j as a dense vector.
matCol :: Matrix a -> Int -> Vector a
matCol mat j
  | j < 0 || j >= matCols mat = error "matCol: index out of bounds"
  | otherwise                 = vecFromList [ matGet mat i j | i <- [0 .. matRows mat - 1] ]

-- | Converts matrix to list of rows.
matToList :: Matrix a -> [[a]]
matToList = mData

-- | Matrix addition: A + B.
matAdd :: Num a => Matrix a -> Matrix a -> Matrix a
matAdd (Matrix r1 c1 rows1) (Matrix r2 c2 rows2)
  | r1 /= r2 || c1 /= c2 = error "matAdd: dimensions do not match"
  | otherwise = Matrix r1 c1 (zipWith (zipWith (+)) rows1 rows2)

-- | Matrix subtraction: A - B.
matSub :: Num a => Matrix a -> Matrix a -> Matrix a
matSub (Matrix r1 c1 rows1) (Matrix r2 c2 rows2)
  | r1 /= r2 || c1 /= c2 = error "matSub: dimensions do not match"
  | otherwise = Matrix r1 c1 (zipWith (zipWith (-)) rows1 rows2)

-- | Scalar multiplication: c * A.
matScale :: Num a => a -> Matrix a -> Matrix a
matScale c (Matrix r cols rows) =
  Matrix r cols (map (map (c *)) rows)

-- | Matrix multiplication: A * B.
matMul :: forall a. Num a => Matrix a -> Matrix a -> Matrix a
matMul (Matrix r1 c1 rows1) (Matrix r2 c2 rows2)
  | c1 /= r2 = error "matMul: inner dimensions do not match"
  | otherwise =
      let cols2 = transpose rows2
          newRows = [ [ sum (zipWith (*) row col) | col <- cols2 ] | row <- rows1 ]
      in Matrix r1 c2 newRows

-- | Matrix-vector multiplication: A * v.
matVecMul :: Num a => Matrix a -> Vector a -> Vector a
matVecMul mat v
  | matCols mat /= length (vecToList v) = error "matVecMul: dimensions do not match"
  | otherwise =
      let xs = vecToList v
          newXs = [ sum (zipWith (*) row xs) | row <- mData mat ]
      in vecFromList newXs

-- | Transpose of a matrix: A^T.
matTranspose :: Matrix a -> Matrix a
matTranspose (Matrix r cols rows) =
  Matrix cols r (transpose rows)

-- | Trace of a square matrix: sum of diagonal entries.
matTrace :: Num a => Matrix a -> a
matTrace mat
  | matRows mat /= matCols mat = error "matTrace: matrix must be square"
  | otherwise = sum [ matGet mat i i | i <- [0 .. matRows mat - 1] ]

-- | Submatrix formed by rows in [rStart .. rEnd] and cols in [cStart .. cEnd].
subMatrix :: Int -> Int -> Int -> Int -> Matrix a -> Matrix a
subMatrix rStart rEnd cStart cEnd mat =
  let rows = [ [ matGet mat i j | j <- [cStart .. cEnd] ] | i <- [rStart .. rEnd] ]
  in matFromRows rows

-- | Minor matrix formed by deleting row delR and column delC.
minorMatrix :: Int -> Int -> Matrix a -> Matrix a
minorMatrix delR delC (Matrix r cols rows) =
  let filteredRows = [ [ val | (j, val) <- zip [0 ..] row, j /= delC ]
                     | (i, row) <- zip [0 ..] rows
                     , i /= delR
                     ]
  in Matrix (r - 1) (cols - 1) filteredRows

-- | Elementary row operation: swaps row i and row j.
swapRows :: Int -> Int -> Matrix a -> Matrix a
swapRows i j mat@(Matrix r cols rows)
  | i == j = mat
  | otherwise =
      let rowI = rows !! i
          rowJ = rows !! j
          newRows = [ if k == i then rowJ else if k == j then rowI else row
                    | (k, row) <- zip [0 .. r - 1] rows
                    ]
      in Matrix r cols newRows

-- | Elementary row operation: scales row i by scalar factor c.
scaleRow :: Num a => Int -> a -> Matrix a -> Matrix a
scaleRow i c (Matrix r cols rows) =
  let newRows = [ if k == i then map (c *) row else row | (k, row) <- zip [0 .. r - 1] rows ]
  in Matrix r cols newRows

-- | Elementary row operation: targetRow <- targetRow + lambda * sourceRow.
addRowMultiple :: Num a => Int -> Int -> a -> Matrix a -> Matrix a
addRowMultiple target source lambda (Matrix r cols rows) =
  let sourceRow = rows !! source
      newRows = [ if k == target
                  then zipWith (\t s -> t + lambda * s) row sourceRow
                  else row
                | (k, row) <- zip [0 .. r - 1] rows
                ]
  in Matrix r cols newRows

instance Show a => Show (Matrix a) where
  show (Matrix r cols rows) =
    "Matrix " ++ show r ++ "x" ++ show cols ++ ":\n" ++
    unlines [ "  [" ++ unwords (map show row) ++ "]" | row <- rows ]
