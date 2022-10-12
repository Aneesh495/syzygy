{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.F4
  ( f4
  , f4WithStats
  , symbolicPreprocessing
  , buildMacaulayMatrix
  , extractF4NewGenerators
  , F4Stats(..)
  , emptyF4Stats
  ) where

import Data.List (sortBy)
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import GHC.Generics (Generic)
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseDivides
  , sparseDiv
  , sparseMaxVar
  , compareSparse
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , isZeroPoly
  , polyTerms
  , polyOrder
  , leadingMonomial
  , polyDegree
  , fromTerms
  )
import Syzygy.Polynomial.Arithmetic
  ( mulTermPoly
  , monicPoly
  )
import Syzygy.Linear.Vector (svecFromList, svecToList)
import Syzygy.Linear.SparseMatrix
  ( SparseMat(..)
  , smatFromRows
  , sparseGaussianElimination
  )
import Syzygy.Grobner.SPolynomial
  ( CriticalPair(..)
  , sPolynomialWithMultipliers
  )
import Syzygy.Grobner.Criteria
  ( updatePairSetGebauerMoller
  , sortPairsBySugar
  )

-- | Performance statistics collected during Faugere F4 execution.
data F4Stats = F4Stats
  { f4Rounds            :: !Int  -- ^ Number of linear algebra matrix elimination rounds
  , f4PairsProcessed    :: !Int  -- ^ Total critical pairs selected
  , f4MacaulayRows      :: !Int  -- ^ Peak rows in a Macaulay matrix
  , f4MacaulayCols      :: !Int  -- ^ Peak columns (monomials) in a Macaulay matrix
  , f4ZeroReductions    :: !Int  -- ^ Rows reducing to 0 during matrix elimination
  , f4NewGenerators     :: !Int  -- ^ New monic generators extracted from reduced matrices
  } deriving (Eq, Show, Generic)

-- | Initial empty statistics for F4.
emptyF4Stats :: F4Stats
emptyF4Stats = F4Stats 0 0 0 0 0 0

-- | Symbolic Preprocessing phase of Faugere F4.
-- Given a set of selected critical pairs and the current basis,
-- generates a set of shifted polynomials L = { u_i * f_i } such that
-- every monomial appearing in L is either irreducible or has a divisor in L.
symbolicPreprocessing :: forall c. (Eq c, Fractional c)
                      => [CriticalPair c] -> [Poly c] -> [Poly c]
symbolicPreprocessing selectedPairs basis =
  let initialL = concatMap getPairGenerators selectedPairs
      initialM = foldl (\s p -> Set.union s (Set.fromList (map getM (polyTerms p))))
                       Set.empty initialL
  in loop (Map.fromList [ (leadingM p, p) | p <- initialL ]) initialM Set.empty
  where
    getM (Term _ m) = m
    leadingM p = case leadingMonomial p of
      Just m  -> m
      Nothing -> SparseMonomial Map.empty

    getPairGenerators cp =
      let p1 = cpPoly1 cp
          p2 = cpPoly2 cp
          (_, t1, t2) = sPolynomialWithMultipliers p1 p2
      in [ mulTermPoly t1 p1, mulTermPoly t2 p2 ]

    loop !polyMap !unprocessedM !processedM
      | Set.null unprocessedM = Map.elems polyMap
      | otherwise =
          let (m, restM) = Set.deleteFindMin unprocessedM
              processedM' = Set.insert m processedM
          in if Map.member m polyMap
             then loop polyMap restM processedM'
             else
               case findReductor m basis of
                 Just reductor ->
                   let lmR = case leadingMonomial reductor of
                               Just lm -> lm
                               Nothing -> error "symbolicPreprocessing: reductor has no leading monomial"
                       quotM = sparseDiv lmR m
                       shiftedReductor = mulTermPoly (Term 1 quotM) reductor
                       newMonomials = Set.fromList [ getM t | t <- polyTerms shiftedReductor ]
                       newUnprocessed = Set.difference (Set.union restM newMonomials) processedM'
                       polyMap' = Map.insert m shiftedReductor polyMap
                   in loop polyMap' newUnprocessed processedM'
                 Nothing ->
                   loop polyMap restM processedM'

    findReductor m bs =
      case [ g | g <- bs, not (isZeroPoly g),
                 case leadingMonomial g of
                   Just lmG -> sparseDivides lmG m
                   Nothing  -> False
           ] of
        (g:_) -> Just g
        []    -> Nothing

-- | Builds the sparse Macaulay matrix for a list of polynomials.
-- Columns correspond to all distinct monomials in descending MonomialOrder.
-- Rows correspond to the polynomials.
buildMacaulayMatrix :: forall c. (Eq c, Fractional c)
                    => [Poly c] -> (SparseMat c, [SparseMonomial])
buildMacaulayMatrix polys =
  let ord = case polys of
              (p:_) -> polyOrder p
              []    -> error "buildMacaulayMatrix: empty polynomial list"
      allMonomials = Set.toList (foldl (\acc p ->
                                          foldl (\s (Term _ m) -> Set.insert m s)
                                                acc (polyTerms p)
                                       ) Set.empty polys)
      numVars = maximum (0 : map sparseMaxVar allMonomials)
      sortedMonomials = sortBy (\m1 m2 ->
                                  case compareSparse ord numVars m1 m2 of
                                    GT -> LT
                                    LT -> GT
                                    EQ -> EQ
                               ) allMonomials
      colIndexMap = Map.fromList (zip sortedMonomials [0 :: Int ..])
      numCols = length sortedMonomials

      makeRow p =
        let entries = [ (Map.findWithDefault 0 m colIndexMap, c)
                      | Term c m <- polyTerms p
                      ]
        in svecFromList numCols entries
      rows = map makeRow polys
  in (smatFromRows numCols rows, sortedMonomials)

-- | Extracts newly reduced basis generators from the reduced Macaulay matrix.
-- A row corresponds to a new generator if its leading monomial was NOT a leading monomial
-- of any row before Gaussian elimination.
extractF4NewGenerators :: forall c. (Eq c, Fractional c)
                       => [Poly c]               -- ^ Polynomials before matrix reduction
                       -> SparseMat c            -- ^ Reduced sparse matrix
                       -> [SparseMonomial]       -- ^ Monomial column index mapping
                       -> [Poly c]               -- ^ Newly discovered polynomials
extractF4NewGenerators inputPolys reducedMat colMonomials =
  let ord = case inputPolys of
              (p:_) -> polyOrder p
              []    -> error "extractF4NewGenerators: empty inputs"
      inputLeadingMonomials = Set.fromList [ lm | p <- inputPolys
                                                , Just lm <- [leadingMonomial p]
                                                ]
      reducedRows = extractBasisRows reducedMat
  in [ p | sv <- reducedRows
         , let p = rowToPoly ord sv
         , not (isZeroPoly p)
         , case leadingMonomial p of
             Just lm -> not (Set.member lm inputLeadingMonomials)
             Nothing -> False
     ]
  where
    extractBasisRows sm =
      let (pivots, _) = sparseGaussianElimination sm
      in Map.elems pivots

    rowToPoly ord sv =
      let entries = svecToList sv
          terms = [ Term val (colMonomials !! col)
                  | (col, val) <- entries
                  , col < length colMonomials
                  , val /= 0
                  ]
      in fromTerms ord terms

-- | Faugere F4 algorithm for computing Groebner bases via linear algebra.
f4 :: (Eq c, Fractional c) => [Poly c] -> [Poly c]
f4 inputPolys = fst (f4WithStats inputPolys)

-- | Faugere F4 algorithm with comprehensive execution statistics.
f4WithStats :: forall c. (Eq c, Fractional c)
            => [Poly c] -> ([Poly c], F4Stats)
f4WithStats inputPolys =
  let nonZeroInputs = [ monicPoly p | p <- inputPolys, not (isZeroPoly p) ]
  in case nonZeroInputs of
       [] -> ([], emptyF4Stats)
       [p] -> ([p], emptyF4Stats)
       (p0:rest) ->
         let initBasis = [p0]
             initSugars = [polyDegree p0]
             (fullBasis, fullSugars, initialQueue) =
               foldl installInitial (initBasis, initSugars, []) (zip [1 :: Int ..] rest)
             sortedQueue = sortPairsBySugar initialQueue
         in loop fullBasis fullSugars sortedQueue emptyF4Stats
  where
    installInitial (currB, currS, currQ) (idx, p) =
      let sugarP = polyDegree p
          newQ = updatePairSetGebauerMoller currB currS p idx sugarP currQ
      in (currB ++ [p], currS ++ [sugarP], newQ)

    loop !basis _ [] !stats = (basis, stats)
    loop !basis !sugars queue !stats =
      -- Selection step: select all pairs of minimal sugar degree d
      let minDeg = minimum (map cpSugar queue)
          (selectedPairs, remainingQueue) = span (\cp -> cpSugar cp == minDeg) queue
      in if null selectedPairs
         then (basis, stats)
         else
           let -- Step 1: Symbolic Preprocessing
               lPolys = symbolicPreprocessing selectedPairs basis
               -- Step 2: Build Macaulay Matrix
               (macaulayMat, colMonomials) = buildMacaulayMatrix lPolys
               -- Step 3: Sparse Gaussian Elimination
               (reducedPivots, _) = sparseGaussianElimination macaulayMat
               reducedMat = smatFromRows (length colMonomials) (Map.elems reducedPivots)
               -- Step 4: Extract New Generators
               newGenerators = extractF4NewGenerators lPolys reducedMat colMonomials
               monicNew = [ monicPoly g | g <- newGenerators, not (isZeroPoly g) ]
           in if null monicNew
              then
                let stats' = stats
                      { f4Rounds = f4Rounds stats + 1
                      , f4PairsProcessed = f4PairsProcessed stats + length selectedPairs
                      , f4ZeroReductions = f4ZeroReductions stats + length selectedPairs
                      }
                in loop basis sugars remainingQueue stats'
              else
                let updateWithGen (currB, currS, currQ) g =
                      let gIdx = length currB
                          gSugar = polyDegree g
                          q' = updatePairSetGebauerMoller currB currS g gIdx gSugar currQ
                      in (currB ++ [g], currS ++ [gSugar], q')
                    (newBasis, newSugars, newQueue) =
                      foldl updateWithGen (basis, sugars, remainingQueue) monicNew
                    stats' = stats
                      { f4Rounds = f4Rounds stats + 1
                      , f4PairsProcessed = f4PairsProcessed stats + length selectedPairs
                      , f4MacaulayRows = max (f4MacaulayRows stats) (length lPolys)
                      , f4MacaulayCols = max (f4MacaulayCols stats) (length colMonomials)
                      , f4NewGenerators = f4NewGenerators stats + length monicNew
                      }
                in loop newBasis newSugars newQueue stats'
