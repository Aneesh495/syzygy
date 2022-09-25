{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Linear.Gaussian
  ( rref
  , ref
  , matRank
  , matInverse
  , solveLinearSystem
  , nullSpace
  , columnSpace
  , bareissDeterminant
  ) where

import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matRows
  , matCols
  , matGet
  , matIdentity
  , swapRows
  , scaleRow
  , addRowMultiple
  )
import Syzygy.Linear.Vector (Vector(..), vecFromList)

-- | Computes the row echelon form (REF) of an m x n matrix over a field.
ref :: forall a. (Eq a, Fractional a) => Matrix a -> Matrix a
ref initialMat = go 0 0 initialMat
  where
    numR = matRows initialMat
    numC = matCols initialMat

    go !r !c mat
      | r >= numR || c >= numC = mat
      | otherwise =
          -- Find pivot row in column c at or below row r
          case findPivot r c mat of
            Nothing -> go r (c + 1) mat
            Just pivotR ->
              let matSwapped = if pivotR /= r then swapRows r pivotR mat else mat
                  pivotVal = matGet matSwapped r c
                  matScaled = scaleRow r (recip pivotVal) matSwapped
                  matElim = eliminateBelow r c matScaled
              in go (r + 1) (c + 1) matElim

    findPivot r c mat =
      let candidates = [ (i, matGet mat i c) | i <- [r .. numR - 1] ]
      in case filter (\(_, val) -> val /= 0) candidates of
           ((pIdx, _):_) -> Just pIdx
           []            -> Nothing

    eliminateBelow r c mat =
      foldl (\acc i ->
               let lead = matGet acc i c
               in if lead == 0
                  then acc
                  else addRowMultiple i r (-lead) acc
            ) mat [r + 1 .. numR - 1]

-- | Computes the reduced row echelon form (RREF) of an m x n matrix over a field
-- using Gauss-Jordan elimination with row pivoting.
rref :: forall a. (Eq a, Fractional a) => Matrix a -> Matrix a
rref initialMat =
  let refMat = ref initialMat
  in backSubstitute (matRows refMat - 1) refMat
  where
    numC = matCols initialMat

    backSubstitute !r mat
      | r < 0 = mat
      | otherwise =
          case findLeadingCol r mat of
            Nothing -> backSubstitute (r - 1) mat
            Just c ->
              let matElimAbove = foldl (\acc i ->
                                         let lead = matGet acc i c
                                         in if lead == 0
                                            then acc
                                            else addRowMultiple i r (-lead) acc
                                       ) mat [0 .. r - 1]
              in backSubstitute (r - 1) matElimAbove

    findLeadingCol r mat =
      let cols = [ (c, matGet mat r c) | c <- [0 .. numC - 1] ]
      in case filter (\(_, val) -> val /= 0) cols of
           ((c, _):_) -> Just c
           []         -> Nothing

-- | Computes the rank of a matrix (number of non-zero rows in RREF).
matRank :: (Eq a, Fractional a) => Matrix a -> Int
matRank mat =
  let rrefMat = rref mat
      nonZeroRows = filter (any (/= 0)) (mData rrefMat)
  in length nonZeroRows

-- | Computes the inverse of a square matrix A via Gauss-Jordan elimination on [A | I].
-- Returns Just A^{-1} if invertible, or Nothing if singular.
matInverse :: forall a. (Eq a, Fractional a) => Matrix a -> Maybe (Matrix a)
matInverse mat
  | numR /= numC = Nothing
  | otherwise =
      let ident = matIdentity numR
          augmented = Matrix numR (2 * numC) [ rA ++ rI | (rA, rI) <- zip (mData mat) (mData ident) ]
          reduced = rref augmented
          leftPart = [ take numC row | row <- mData reduced ]
          rightPart = [ drop numC row | row <- mData reduced ]
      in if leftPart == mData ident
         then Just (Matrix numR numC rightPart)
         else Nothing
  where
    numR = matRows mat
    numC = matCols mat

-- | Solves the linear system A * x = b over a field.
-- Returns Just x if a unique solution exists, or Nothing if inconsistent / underdetermined.
solveLinearSystem :: (Eq a, Fractional a) => Matrix a -> Vector a -> Maybe (Vector a)
solveLinearSystem a b
  | matRows a /= length (vElements b) = Nothing
  | otherwise =
      let aug = Matrix (matRows a) (matCols a + 1)
                  [ row ++ [val] | (row, val) <- zip (mData a) (vElements b) ]
          red = rref aug
          r = matRank a
          augR = length (filter (any (/= 0)) (mData red))
      in if r /= augR
         then Nothing  -- Inconsistent system
         else if r < matCols a
         then Nothing  -- Infinitely many solutions (not unique)
         else
           let sol = [ last row | row <- take (matCols a) (mData red) ]
           in Just (vecFromList sol)

-- | Computes a basis for the nullspace (kernel) of matrix A:
-- { x | A * x = 0 }.
nullSpace :: forall a. (Eq a, Fractional a) => Matrix a -> [Vector a]
nullSpace mat =
  let red = rref mat
      numC = matCols mat
      pivotCols = [ findPivotCol r | r <- mData red, any (/= 0) r ]
      freeCols = [ c | c <- [0 .. numC - 1], c `notElem` pivotCols ]
  in [ makeKernelVector freeC red pivotCols numC | freeC <- freeCols ]
  where
    findPivotCol row =
      let nonZeros = [ idx | (idx, val) <- zip [0 ..] row, val /= 0 ]
      in case nonZeros of
           (nz:_) -> nz
           []     -> -1

    makeKernelVector freeC red pivotCols numC =
      let vecEntries = [ computeEntry i | i <- [0 .. numC - 1] ]
          computeEntry i
            | i == freeC = 1
            | i `elem` pivotCols =
                case [ - (row !! freeC) | row <- mData red, findPivotCol row == i ] of
                  (val:_) -> val
                  []      -> 0
            | otherwise = 0
      in vecFromList vecEntries

-- | Computes a basis for the column space (image) of matrix A:
-- extracts the columns corresponding to the pivot columns in RREF.
columnSpace :: (Eq a, Fractional a) => Matrix a -> [Vector a]
columnSpace mat =
  let red = rref mat
      pivotCols = [ findPivotCol row | row <- mData red, any (/= 0) row ]
  in [ vecFromList [ matGet mat i c | i <- [0 .. matRows mat - 1] ] | c <- pivotCols ]
  where
    findPivotCol row =
      let nonZeros = [ idx | (idx, val) <- zip [0 ..] row, val /= 0 ]
      in case nonZeros of
           (nz:_) -> nz
           []     -> -1

-- | Fraction-free Bareiss algorithm for computing the determinant of an n x n matrix
-- over any integral domain without fraction swell in O(n^3) time.
-- All divisions in the Bareiss recursion are guaranteed to have remainder zero!
bareissDeterminant :: forall a. (Eq a, Integral a) => Matrix a -> a
bareissDeterminant mat
  | numR /= numC = error "bareissDeterminant: matrix must be square"
  | numR == 0 = 1
  | numR == 1 = matGet mat 0 0
  | otherwise = go 0 1 (mData mat) 1
  where
    numR = matRows mat
    numC = matCols mat

    go !k !signVal rows !prevPivot
      | k >= numR - 1 =
          signVal * ((rows !! (numR - 1)) !! (numC - 1))
      | otherwise =
          -- Find pivot in column k at or below row k
          case [ i | i <- [k .. numR - 1], (rows !! i) !! k /= 0 ] of
            [] -> 0  -- Singular matrix: determinant is 0
            (pivotI:_) ->
              let (newRows, newSign) = if pivotI /= k
                                       then (swapList k pivotI rows, -signVal)
                                       else (rows, signVal)
                  p = (newRows !! k) !! k
                  elimRow i =
                    let rowI = newRows !! i
                        rowK = newRows !! k
                        valIK = rowI !! k
                    in [ if j < k + 1
                         then 0
                         else ((rowI !! j) * p - (rowK !! j) * valIK) `quot` prevPivot
                       | j <- [0 .. numC - 1]
                       ]
                  nextRows = take (k + 1) newRows ++ [ elimRow i | i <- [k + 1 .. numR - 1] ]
              in go (k + 1) newSign nextRows p

    swapList i j xs
      | i == j = xs
      | otherwise =
          let elemI = xs !! i
              elemJ = xs !! j
          in [ if k == i then elemJ else if k == j then elemI else x | (k, x) <- zip [0 ..] xs ]
