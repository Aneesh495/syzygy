{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Linear.SmithNormalForm
  ( SNF(..)
  , smithNormalForm
  , invariantFactors
  , abelianGroupStructure
  , verifySNF
  , swapCols
  , scaleCol
  , addColMultiple
  ) where

import GHC.Generics (Generic)
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matRows
  , matCols
  , matGet
  , matIdentity
  , matMul
  , swapRows
  , scaleRow
  , addRowMultiple
  )

-- | Result of Smith Normal Form decomposition: U * A * V = D
-- where D is diagonal with d_1 | d_2 | ... | d_r, and U, V are unimodular.
data SNF a = SNF
  { snfD :: !(Matrix a)  -- ^ Diagonal matrix of invariant factors
  , snfU :: !(Matrix a)  -- ^ Unimodular row transformation matrix
  , snfV :: !(Matrix a)  -- ^ Unimodular column transformation matrix
  } deriving (Eq, Show, Generic)

-- | Extracts the non-zero invariant factors d_1, ..., d_r from the SNF diagonal.
invariantFactors :: SNF Integer -> [Integer]
invariantFactors (SNF d _ _) =
  let r = min (matRows d) (matCols d)
      diags = [ abs (matGet d i i) | i <- [0 .. r - 1] ]
  in filter (/= 0) diags

-- | Computes the abstract finitely generated abelian group decomposition:
-- coker(A) = Z^betti + Z_{t_1} + ... + Z_{t_k}.
-- Returns (freeRank, torsionCoefficients).
abelianGroupStructure :: SNF Integer -> (Int, [Integer])
abelianGroupStructure snf =
  let invs = invariantFactors snf
      torsion = filter (> 1) invs
      rank = length invs
      freeRank = max 0 (matRows (snfD snf) - rank)
  in (freeRank, torsion)

-- | Elementary column operations on a matrix (acting by post-multiplication).
swapCols :: Int -> Int -> Matrix a -> Matrix a
swapCols j1 j2 (Matrix r cols rows)
  | j1 == j2 = Matrix r cols rows
  | otherwise =
      let swapRow rList =
            let v1 = rList !! j1
                v2 = rList !! j2
            in [ if k == j1 then v2 else if k == j2 then v1 else v
               | (k, v) <- zip [0 .. cols - 1] rList
               ]
      in Matrix r cols (map swapRow rows)

scaleCol :: Num a => Int -> a -> Matrix a -> Matrix a
scaleCol j c (Matrix r cols rows) =
  let scaleR rList =
        [ if k == j then c * v else v | (k, v) <- zip [0 .. cols - 1] rList ]
  in Matrix r cols (map scaleR rows)

addColMultiple :: Num a => Int -> Int -> a -> Matrix a -> Matrix a
addColMultiple target source lambda (Matrix r cols rows) =
  let addR rList =
        let sVal = rList !! source
        in [ if k == target then v + lambda * sVal else v
           | (k, v) <- zip [0 .. cols - 1] rList
           ]
  in Matrix r cols (map addR rows)

-- | Computes the exact Smith Normal Form of an integer matrix A: U * A * V = D
-- where D is diagonal, U and V are unimodular, and d_1 | d_2 | ... | d_r.
smithNormalForm :: forall a. (Eq a, Integral a) => Matrix a -> SNF a
smithNormalForm initialA =
  let uInit = matIdentity m
      vInit = matIdentity n
  in go 0 initialA uInit vInit
  where
    m = matRows initialA
    n = matCols initialA

    go !step mat u v
      | step >= min m n =
          let (dFinal, uFinal) = makeDiagonalPositive mat u
          in SNF dFinal uFinal v
      | otherwise =
          -- Find non-zero pivot entry in submatrix [step .. m-1] x [step .. n-1]
          case findPivot step mat of
            Nothing ->
              let (dFinal, uFinal) = makeDiagonalPositive mat u
              in SNF dFinal uFinal v
            Just (pRow, pCol) ->
              let -- Move pivot to position (step, step)
                  mat1 = if pRow /= step then swapRows step pRow mat else mat
                  u1   = if pRow /= step then swapRows step pRow u else u
                  mat2 = if pCol /= step then swapCols step pCol mat1 else mat1
                  v2   = if pCol /= step then swapCols step pCol v else v
                  -- Clear column and row at step
                  (mat3, u3, v3) = clearStep step mat2 u1 v2
              in go (step + 1) mat3 u3 v3

    findPivot step mat =
      let nonZeros = [ (i, j, abs (matGet mat i j))
                     | i <- [step .. m - 1]
                     , j <- [step .. n - 1]
                     , matGet mat i j /= 0
                     ]
      in case nonZeros of
           [] -> Nothing
           _  ->
             -- Choose pivot with minimal absolute value to accelerate Euclidean steps
             let (minI, minJ, _) = foldl1' (\acc@(_, _, v1) cand@(_, _, v2) ->
                                              if v2 < v1 then cand else acc
                                           ) nonZeros
             in Just (minI, minJ)

    foldl1' _ [] = error "foldl1': empty list"
    foldl1' f (x:xs) = foldl f x xs

    clearStep step mat u v =
      -- Iteratively clear entries in row step and col step until only diagonal entry remains
      let (matR, uR) = eliminateCol step mat u
          (matC, vC) = eliminateRow step matR v
      in if any (\i -> matGet matC i step /= 0) [step + 1 .. m - 1] ||
            any (\j -> matGet matC step j /= 0) [step + 1 .. n - 1]
         then clearStep step matC uR vC
         else (matC, uR, vC)

    eliminateCol step mat u =
      foldl (\(curM, curU) i ->
               let pivotVal = matGet curM step step
                   targetVal = matGet curM i step
               in if targetVal == 0
                  then (curM, curU)
                  else
                    let q = targetVal `div` pivotVal
                        nextM = addRowMultiple i step (-q) curM
                        nextU = addRowMultiple i step (-q) curU
                    in (nextM, nextU)
            ) (mat, u) [step + 1 .. m - 1]

    eliminateRow step mat v =
      foldl (\(curM, curV) j ->
               let pivotVal = matGet curM step step
                   targetVal = matGet curM step j
               in if targetVal == 0
                  then (curM, curV)
                  else
                    let q = targetVal `div` pivotVal
                        nextM = addColMultiple j step (-q) curM
                        nextV = addColMultiple j step (-q) curV
                    in (nextM, nextV)
            ) (mat, v) [step + 1 .. n - 1]

    makeDiagonalPositive mat u =
      let r = min m n
      in foldl (\(curM, curU) i ->
                  let diagVal = matGet curM i i
                  in if diagVal < 0
                     then (scaleRow i (-1) curM, scaleRow i (-1) curU)
                     else (curM, curU)
               ) (mat, u) [0 .. r - 1]

-- | Verifies that U * A * V == D and all invariant factors satisfy d_i | d_{i+1}.
verifySNF :: forall a. (Eq a, Integral a) => Matrix a -> SNF a -> Bool
verifySNF origA (SNF d u v) =
  let prod = matMul (matMul u origA) v
      isEqual = prod == d
      diags = [ matGet d i i | i <- [0 .. min (matRows d) (matCols d) - 1], matGet d i i /= 0 ]
      dividesConsecutive [] = True
      dividesConsecutive [_] = True
      dividesConsecutive (d1:d2:rest) = (d2 `rem` d1 == 0) && dividesConsecutive (d2:rest)
  in isEqual && dividesConsecutive diags
