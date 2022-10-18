{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.FGLM
  ( isZeroDimensional
  , quotientStandardBasis
  , quotientDimension
  , fglm
  , fglmToLex
  ) where

import qualified Data.List as List
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Data.Maybe (mapMaybe)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  , sparseMul
  , sparseDivides
  , sparseMaxVar
  , compareSparse
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , isZeroPoly
  , leadingMonomial
  , polyOrder
  )
import Syzygy.Polynomial.Division (normalForm)

-- | Checks whether an ideal represented by a Groebner basis G is zero-dimensional.
-- Over an algebraically closed field, an ideal is zero-dimensional (has finitely many solutions)
-- if and only if for every variable x_i in the polynomial ring, there exists an element
-- in the Groebner basis whose leading monomial is a pure power of x_i.
isZeroDimensional :: (Eq a, Fractional a) => Int -> [Poly a] -> Bool
isZeroDimensional numVars gb =
  let nonZeroGB = filter (not . isZeroPoly) gb
      lms = mapMaybe leadingMonomial nonZeroGB
  in all (hasPurePower lms) [0 .. numVars - 1]
  where
    hasPurePower lms varId =
      any (\(SparseMonomial m) ->
            case Map.toList m of
              [(v, expVal)] -> v == varId && expVal > 0
              _             -> False
          ) lms

-- | Computes the standard monomial basis of the quotient ring A = k[x_0, ..., x_{n-1}] / I.
-- These are all monomials that are not divisible by any leading monomial of the Groebner basis G.
-- For a zero-dimensional ideal, this basis is finite and its cardinality equals
-- the vector space dimension of A.
quotientStandardBasis :: (Eq a, Fractional a) => Int -> [Poly a] -> [SparseMonomial]
quotientStandardBasis numVars gb =
  let nonZeroGB = filter (not . isZeroPoly) gb
      lms = mapMaybe leadingMonomial nonZeroGB
  in bfs lms Set.empty [sparseOne] []
  where
    bfs _ _ [] acc = reverse acc
    bfs lms !visited (m:q) acc
      | Set.member m visited = bfs lms visited q acc
      | any (`sparseDivides` m) lms = bfs lms (Set.insert m visited) q acc
      | otherwise =
          let visited' = Set.insert m visited
              children = [ sparseMul m (sparseVar v 1) | v <- [0 .. numVars - 1] ]
              q' = q ++ filter (\c -> not (Set.member c visited')) children
          in bfs lms visited' q' (m : acc)

-- | Returns the vector space dimension of the quotient ring k[x_0, ..., x_{n-1}] / I.
-- Returns Nothing if the ideal is positive-dimensional.
quotientDimension :: (Eq a, Fractional a) => Int -> [Poly a] -> Maybe Int
quotientDimension numVars gb
  | isZeroDimensional numVars gb = Just (length (quotientStandardBasis numVars gb))
  | otherwise = Nothing

-- | Given linearly independent column vectors v_1, ..., v_k in k^D and a vector v in k^D,
-- determines if v lies in their linear span.
-- If so, returns the unique coordinates [c_1, ..., c_k] such that sum_{i=1}^k c_i * v_i = v.
-- If v is linearly independent from {v_1, ..., v_k}, returns Nothing.
findCoordinates :: forall a. (Eq a, Fractional a) => [[a]] -> [a] -> Maybe [a]
findCoordinates cols targetVec
  | null cols =
      if all (== 0) targetVec
      then Just []
      else Nothing
  | otherwise =
      let numRows = length targetVec
          numCols = length cols
          -- Construct augmented matrix of size numRows x (numCols + 1)
          initialMat = [ [ col !! r | col <- cols ] ++ [targetVec !! r] | r <- [0 .. numRows - 1] ]
      in gaussianElim 0 0 initialMat numRows numCols
  where
    gaussianElim !r !c mat numRows numCols
      | c >= numCols =
          -- All basis columns eliminated; check if remaining rows are consistent (all zeros in target column)
          let remainingInconsistent = any (\row -> (row !! numCols) /= 0) (drop numCols mat)
          in if remainingInconsistent
             then Nothing
             else Just [ (mat !! i) !! numCols | i <- [0 .. numCols - 1] ]
      | r >= numRows = Nothing
      | otherwise =
          case findPivot r c mat numRows of
            Nothing ->
              -- Should not happen since cols are linearly independent, but handle safely
              Nothing
            Just pivotR ->
              let matSwapped = if pivotR /= r then swapRows r pivotR mat else mat
                  pivotVal = (matSwapped !! r) !! c
                  matScaled = scaleRow r (recip pivotVal) matSwapped
                  matElim = elimRows r c matScaled numRows
              in gaussianElim (r + 1) (c + 1) matElim numRows numCols

    findPivot r c mat numRows =
      let candidates = [ (i, (mat !! i) !! c) | i <- [r .. numRows - 1] ]
      in case filter (\(_, val) -> val /= 0) candidates of
           ((idx, _):_) -> Just idx
           []           -> Nothing

    swapRows r1 r2 mat =
      [ if i == r1 then mat !! r2
        else if i == r2 then mat !! r1
        else mat !! i
      | i <- [0 .. length mat - 1]
      ]

    scaleRow r factor mat =
      [ if i == r then map (* factor) (mat !! i) else mat !! i
      | i <- [0 .. length mat - 1]
      ]

    elimRows r c mat numRows =
      [ if i == r
        then mat !! i
        else let factor = (mat !! i) !! c
             in if factor == 0
                then mat !! i
                else zipWith (\x y -> x - factor * y) (mat !! i) (mat !! r)
      | i <- [0 .. numRows - 1]
      ]

-- | Converts a polynomial to its coordinate vector with respect to the quotient basis B.
polyToVector :: (Eq a, Fractional a) => [SparseMonomial] -> Poly a -> [a]
polyToVector basis (Poly terms _) =
  let termMap = Map.fromList [ (m, c) | Term c m <- terms ]
  in [ Map.findWithDefault 0 b termMap | b <- basis ]

-- | The FGLM algorithm (Faugere, Gianni, Lazard, Mora 1993).
-- Converts a Groebner basis G with respect to an initial term order (typically DegRevLex)
-- into the unique reduced Groebner basis with respect to a target term order (typically Lex).
-- Runs in O(n * D^3) time where D is the vector space dimension of the quotient ring,
-- bypassing the exponential complexity of Buchberger's algorithm for lexicographic bases.
fglm :: forall a. (Ord a, Fractional a)
     => MonomialOrder
     -> [Poly a]
     -> Maybe [Poly a]
fglm targetOrder gb =
  let nonZeroGB = filter (not . isZeroPoly) gb
      maxVarInGB = maximum (-1 : [ sparseMaxVar m | g <- nonZeroGB, Just m <- [leadingMonomial g] ])
      numVars = maxVarInGB + 1
      origOrder = case nonZeroGB of
                    (g:_) -> polyOrder g
                    []    -> DegRevLex
  in if numVars <= 0
     then Just nonZeroGB
     else if not (isZeroDimensional numVars nonZeroGB)
          then Nothing
          else
            let stdBasis = quotientStandardBasis numVars nonZeroGB
                dimD = length stdBasis
            in if dimD == 0
               then Just [Poly [Term 1 sparseOne] targetOrder]
               else Just (fglmLoop numVars origOrder stdBasis nonZeroGB)
  where
    fglmLoop numVars origOrder stdBasis nonZeroGB =
      go [] [] [] (Set.singleton sparseOne) [sparseOne]
      where
        -- cmpTarget orders monomials ascending (smallest first)
        cmpTarget m1 m2 = compareSparse targetOrder numVars m1 m2

        insertCandidate m [] = [m]
        insertCandidate m (x:xs) =
          case cmpTarget m x of
            LT -> m : x : xs
            EQ -> x : xs
            GT -> x : insertCandidate m xs

        go !targetGB !targetLMs !linIndepCols !visitedQueue candidateQueue =
          case candidateQueue of
            [] -> reverse targetGB
            (m:remCandidates)
              -- If m is divisible by any already-found leading monomial in the target order, skip it
              | any (`sparseDivides` m) targetLMs ->
                  go targetGB targetLMs linIndepCols visitedQueue remCandidates
              | otherwise ->
                  let -- Compute normal form of monomial m with respect to the original Groebner basis
                      pM = Poly [Term 1 m] origOrder
                      nfM = normalForm pM nonZeroGB
                      vecM = polyToVector stdBasis nfM
                      existingCols = map snd linIndepCols
                  in case findCoordinates existingCols vecM of
                       Nothing ->
                         -- m is linearly independent! Add (m, vecM) to the linear basis
                         let linIndepCols' = linIndepCols ++ [(m, vecM)]
                             -- Enqueue new boundary candidates m * x_i
                             newCandidates = [ sparseMul m (sparseVar v 1) | v <- [0 .. numVars - 1] ]
                             (visited', q') = foldl (\(!vAcc, !qAcc) cand ->
                                                      if Set.member cand vAcc || any (`sparseDivides` cand) targetLMs
                                                      then (vAcc, qAcc)
                                                      else (Set.insert cand vAcc, insertCandidate cand qAcc)
                                                    ) (visitedQueue, remCandidates) newCandidates
                         in go targetGB targetLMs linIndepCols' visited' q'

                       Just coords ->
                         -- m is linearly dependent on earlier monomials!
                         -- Produces a new Groebner basis element: m - sum(c_i * m_i) = 0
                         let terms = Term 1 m : [ Term (-c) m_i | (c, (m_i, _)) <- zip coords linIndepCols, c /= 0 ]
                             newPoly = Poly (List.sortBy (\(Term _ ma) (Term _ mb) -> cmpTarget mb ma) terms) targetOrder
                             targetGB' = newPoly : targetGB
                             targetLMs' = m : targetLMs
                         in go targetGB' targetLMs' linIndepCols visitedQueue remCandidates

-- | Convenience helper to convert a Groebner basis to pure Lexicographical order.
fglmToLex :: (Ord a, Fractional a) => [Poly a] -> Maybe [Poly a]
fglmToLex = fglm Lex
