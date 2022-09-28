{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Linear.SparseMatrix
  ( SparseMat(..)
  , COO(..)
  , CSR(..)
  , smatFromCOO
  , smatToCOO
  , smatFromRows
  , smatToRows
  , smatFromCSR
  , smatToCSR
  , smatZero
  , smatRows
  , smatCols
  , smatGet
  , smatSet
  , smatAdd
  , smatScale
  , smatMulVec
  , smatMul
  , smatRow
  , smatEliminate
  , sparseGaussianElimination
  , extractReducedBasisRows
  ) where

import qualified Data.Map.Strict as Map
import GHC.Generics (Generic)
import Syzygy.Linear.Vector
  ( SparseVector(..)
  , svecFromList
  , svecToList
  , svecDim
  , svecIndex
  , svecZero
  , svecAdd
  , svecScale
  , svecDot
  )

-- | Coordinate format (COO) for sparse matrices: (row, col, value) triples.
data COO a = COO
  { cooRows    :: !Int
  , cooCols    :: !Int
  , cooEntries :: ![(Int, Int, a)]
  } deriving (Eq, Show, Generic)

-- | Compressed Sparse Row (CSR) format: values, column indices, and row pointers.
data CSR a = CSR
  { csrRows    :: !Int
  , csrCols    :: !Int
  , csrValues  :: ![a]
  , csrColInd  :: ![Int]
  , csrRowPtr  :: ![Int]
  } deriving (Eq, Show, Generic)

-- | Sparse matrix represented as a Map from row index to SparseVector.
-- Rows with all zero entries are omitted from the map.
data SparseMat a = SparseMat
  { smRows :: !Int
  , smCols :: !Int
  , smData :: !(Map.Map Int (SparseVector a))
  } deriving (Eq, Generic)

-- | Number of rows.
smatRows :: SparseMat a -> Int
smatRows = smRows

-- | Number of columns.
smatCols :: SparseMat a -> Int
smatCols = smCols

-- | Empty zero sparse matrix.
smatZero :: Int -> Int -> SparseMat a
smatZero r c = SparseMat r c Map.empty

-- | Constructs a SparseMat from COO triples.
smatFromCOO :: (Eq a, Num a) => COO a -> SparseMat a
smatFromCOO (COO r c entries) =
  let grouped = foldl insertEntry Map.empty entries
      insertEntry acc (row, col, val)
        | val == 0  = acc
        | otherwise = Map.insertWith (++) row [(col, val)] acc
      svecMap = Map.map (svecFromList c) grouped
  in SparseMat r c svecMap

-- | Converts SparseMat to COO triples.
smatToCOO :: SparseMat a -> COO a
smatToCOO (SparseMat r c rowMap) =
  let triples = [ (row, col, val)
                | (row, svec) <- Map.toList rowMap
                , (col, val) <- svecToList svec
                ]
  in COO r c triples

-- | Constructs SparseMat from a list of SparseVectors representing rows.
smatFromRows :: (Eq a, Num a) => Int -> [SparseVector a] -> SparseMat a
smatFromRows c rows =
  let r = length rows
      nonZeroRows = [ (idx, sv)
                    | (idx, sv) <- zip [0 ..] rows
                    , not (Map.null (svEntries sv))
                    ]
  in SparseMat r c (Map.fromList nonZeroRows)

-- | Converts SparseMat to a list of SparseVectors for each row in [0 .. r-1].
smatToRows :: SparseMat a -> [SparseVector a]
smatToRows (SparseMat r c rowMap) =
  [ Map.findWithDefault (svecZero c) i rowMap | i <- [0 .. r - 1] ]

-- | Converts SparseMat to CSR format.
smatToCSR :: SparseMat a -> CSR a
smatToCSR sm@(SparseMat r c _) =
  let rows = smatToRows sm
      entriesPerRow = map svecToList rows
      vals = concatMap (map snd) entriesPerRow
      colIndices = concatMap (map fst) entriesPerRow
      rowCounts = map length entriesPerRow
      rowPointers = scanl (+) 0 rowCounts
  in CSR r c vals colIndices rowPointers

-- | Converts CSR to SparseMat format.
smatFromCSR :: (Eq a, Num a) => CSR a -> SparseMat a
smatFromCSR (CSR _ c vals colInd rowPtr) =
  let rowSlices = zipWith (\start end -> (start, end)) rowPtr (drop 1 rowPtr)
      buildRow (start, end) =
        let rowVals = take (end - start) (drop start vals)
            rowCols = take (end - start) (drop start colInd)
        in svecFromList c (zip rowCols rowVals)
      rows = map buildRow rowSlices
  in smatFromRows c rows

-- | 0-indexed element lookup at (row i, col j).
smatGet :: Num a => SparseMat a -> Int -> Int -> a
smatGet (SparseMat r c rowMap) i j
  | i < 0 || i >= r = error "smatGet: row out of bounds"
  | j < 0 || j >= c = error "smatGet: col out of bounds"
  | otherwise =
      case Map.lookup i rowMap of
        Just svec -> svecIndex svec j
        Nothing   -> 0

-- | Sets element at (row i, col j).
smatSet :: (Eq a, Num a) => Int -> Int -> a -> SparseMat a -> SparseMat a
smatSet i j val (SparseMat r c rowMap)
  | i < 0 || i >= r = error "smatSet: row out of bounds"
  | j < 0 || j >= c = error "smatSet: col out of bounds"
  | otherwise =
      let currentSvec = Map.findWithDefault (svecZero c) i rowMap
          newEntries = if val == 0
                       then Map.delete j (svEntries currentSvec)
                       else Map.insert j val (svEntries currentSvec)
          newSvec = SparseVector c newEntries
          newMap = if Map.null newEntries
                   then Map.delete i rowMap
                   else Map.insert i newSvec rowMap
      in SparseMat r c newMap

-- | Extracts row i as a SparseVector.
smatRow :: SparseMat a -> Int -> SparseVector a
smatRow (SparseMat r c rowMap) i
  | i < 0 || i >= r = error "smatRow: row out of bounds"
  | otherwise       = Map.findWithDefault (svecZero c) i rowMap

-- | Matrix addition: A + B.
smatAdd :: (Eq a, Num a) => SparseMat a -> SparseMat a -> SparseMat a
smatAdd (SparseMat r1 c1 m1) (SparseMat r2 c2 m2)
  | r1 /= r2 || c1 /= c2 = error "smatAdd: dimension mismatch"
  | otherwise =
      let unionMap = Map.filter (\sv -> not (Map.null (svEntries sv)))
                                (Map.unionWith svecAdd m1 m2)
      in SparseMat r1 c1 unionMap

-- | Scalar multiplication: c * A.
smatScale :: (Eq a, Num a) => a -> SparseMat a -> SparseMat a
smatScale 0 (SparseMat r c _) = smatZero r c
smatScale factor (SparseMat r c rowMap) =
  let scaledMap = Map.map (svecScale factor) rowMap
  in SparseMat r c scaledMap

-- | Matrix-vector product: A * v.
smatMulVec :: (Eq a, Num a) => SparseMat a -> SparseVector a -> SparseVector a
smatMulVec (SparseMat r c rowMap) v
  | c /= svecDim v = error "smatMulVec: dimension mismatch"
  | otherwise =
      let nonZeroProducts = [ (i, svecDot svec v)
                            | (i, svec) <- Map.toList rowMap
                            ]
          entries = filter (\(_, val) -> val /= 0) nonZeroProducts
      in SparseVector r (Map.fromList entries)

-- | Sparse matrix multiplication: A * B.
smatMul :: forall a. (Eq a, Num a) => SparseMat a -> SparseMat a -> SparseMat a
smatMul (SparseMat r1 c1 m1) b@(SparseMat r2 c2 _)
  | c1 /= r2 = error "smatMul: inner dimensions mismatch"
  | otherwise =
      let bCols = [ SparseVector r2 (extractCol j) | j <- [0 .. c2 - 1] ]
          extractCol j = Map.fromList [ (i, smatGet b i j)
                                      | i <- [0 .. r2 - 1]
                                      , smatGet b i j /= 0
                                      ]
          multiplyRow svec =
            let rowEntries = [ (j, svecDot svec (bCols !! j))
                             | j <- [0 .. c2 - 1]
                             , svecDot svec (bCols !! j) /= 0
                             ]
            in svecFromList c2 rowEntries
          newRowMap = Map.filter (\sv -> not (Map.null (svEntries sv)))
                                 (Map.map multiplyRow m1)
      in SparseMat r1 c2 newRowMap

-- | Elementary row elimination: targetRow <- targetRow - lambda * pivotRow.
smatEliminate :: (Eq a, Num a) => a -> SparseVector a -> SparseVector a -> SparseVector a
smatEliminate lambda pivotRow targetRow =
  svecAdd targetRow (svecScale (-lambda) pivotRow)

-- | High-performance Sparse Gaussian Elimination on Macaulay matrices.
-- Given a sparse matrix of candidate polynomial multiples,
-- reduces rows against known pivots and each other using Markowitz-inspired pivoting.
-- Returns (pivots, newlyReducedRows).
sparseGaussianElimination :: forall a. (Eq a, Fractional a)
                          => SparseMat a
                          -> (Map.Map Int (SparseVector a), [SparseVector a])
sparseGaussianElimination sm =
  let rows = smatToRows sm
      nonZeroRows = filter (\sv -> not (Map.null (svEntries sv))) rows
  in go nonZeroRows Map.empty []
  where
    go [] !pivots !newRows = (pivots, reverse newRows)
    go (r:rs) !pivots !newRows =
      -- Reduce row r against existing pivots
      let reducedR = reduceAgainstPivots r pivots
      in if Map.null (svEntries reducedR)
         then go rs pivots newRows  -- Reduces to zero: redundant
         else
           let leadCol = fst (Map.findMin (svEntries reducedR))
               leadVal = svecIndex reducedR leadCol
               normalizedR = svecScale (recip leadVal) reducedR
               pivots' = Map.insert leadCol normalizedR pivots
           in go rs pivots' (normalizedR : newRows)

    reduceAgainstPivots sv pivots
      | Map.null (svEntries sv) = sv
      | otherwise =
          let leadCol = fst (Map.findMin (svEntries sv))
          in case Map.lookup leadCol pivots of
               Just pivotSv ->
                 let leadVal = svecIndex sv leadCol
                     elimSv = smatEliminate leadVal pivotSv sv
                 in reduceAgainstPivots elimSv pivots
               Nothing -> sv

-- | Extracts newly reduced basis rows that cannot be reduced further by input pivots.
extractReducedBasisRows :: (Eq a, Fractional a) => SparseMat a -> [SparseVector a]
extractReducedBasisRows sm =
  let (_, newRows) = sparseGaussianElimination sm
  in newRows

instance Show a => Show (SparseMat a) where
  show sm@(SparseMat r c rowMap) =
    "SparseMat " ++ show r ++ "x" ++ show c ++ " (non-zero rows=" ++ show (Map.size rowMap) ++ "):\n" ++
    unlines [ "  row " ++ show i ++ ": " ++ show (svecToList sv)
            | (i, sv) <- Map.toList (smData sm)
            ]
