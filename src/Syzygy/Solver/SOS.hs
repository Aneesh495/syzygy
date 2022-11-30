{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.SOS
  ( SOSDecomposition(..)
  , isSOSPolynomial
  , decomposeSOS
  , monomialBasisUpToDegree
  , motzkinPolynomial
  , gramMatrixForPoly
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseVar
  , sparseMul
  , sparseDegree
  , sparseMaxVar
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , polyDegree
  )
import Syzygy.Linear.Matrix (Matrix(..))

-- | Represents an exact Sum of Squares decomposition:
-- p(x) = sum_{k=1}^r c_k * q_k(x)^2 where c_k > 0.
data SOSDecomposition a = SOSDecomposition
  { sosComponents :: ![(a, Poly a)] -- ^ Pairs (c_k, q_k(x)) such that sum c_k * q_k^2 == p
  , sosGramMatrix :: !(Matrix a)     -- ^ Positive semidefinite Gram matrix Q
  , sosBasis      :: ![SparseMonomial] -- ^ Monomial basis m_d(x)
  } deriving (Eq, Show, Generic)

-- | Generates all sparse monomials in numVars variables of degree at most maxDeg.
monomialBasisUpToDegree :: Int -> Int -> [SparseMonomial]
monomialBasisUpToDegree numVars maxDeg =
  filter (\m -> sparseDegree m <= maxDeg) (generate 0)
  where
    generate varId
      | varId >= numVars = [sparseOne]
      | otherwise =
          [ sparseMul (sparseVar varId e) rest
          | e <- [0 .. maxDeg]
          , rest <- generate (varId + 1)
          , sparseDegree rest + e <= maxDeg
          ]

-- | Computes the candidate Gram matrix Q for polynomial p with respect to monomial basis m_d(x):
-- Equates p(x) = sum_{i,j} Q_{i,j} * m_i * m_j.
gramMatrixForPoly :: forall a. (Ord a, Fractional a)
                  => Poly a
                  -> [SparseMonomial]
                  -> Matrix a
gramMatrixForPoly (Poly terms _) basis =
  let termMap = Map.fromList [ (m, c) | Term c m <- terms ]
      n = length basis
      -- Collect all pairs (i, j) contributing to each monomial
      productMap = foldl (\acc (i, bi) ->
                           foldl (\accInner (j, bj) ->
                                   let prod = sparseMul bi bj
                                   in Map.insertWith (++) prod [(i, j)] accInner
                                 ) acc (zip [0 :: Int ..] basis)
                         ) Map.empty (zip [0 :: Int ..] basis)
      -- Construct diagonal or symmetric split of coefficients
      entries = [ [ entryVal i j productMap termMap | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]
  in Matrix n n entries
  where
    entryVal i j productMap termMap =
      let bi = basis !! i
          bj = basis !! j
          prod = sparseMul bi bj
          totalCoeff = Map.findWithDefault 0 prod termMap
          numPairs = case Map.lookup prod productMap of
                       Just pairs -> length pairs
                       Nothing    -> 1
      in totalCoeff / fromIntegral numPairs

-- | Checks whether a matrix is positive semidefinite using LDL^T Cholesky-like factorization:
-- Q = L * D * L^T where all diagonal entries D_{i,i} >= 0.
choleskyLDLT :: forall a. (Ord a, Fractional a) => Matrix a -> Maybe [(a, [a])]
choleskyLDLT (Matrix n _ d) =
  go 0 d []
  where
    go !i mat acc
      | i >= n = Just (reverse acc)
      | otherwise =
          let pivot = (mat !! i) !! i
          in if pivot < 0
             then Nothing
             else if pivot == 0
                  then
                    -- Degenerate zero pivot: check if whole row is zero
                    let rowNonZero = any (/= 0) (drop i (mat !! i))
                    in if rowNonZero
                       then Nothing
                       else go (i + 1) mat ((0, replicate n 0) : acc)
                  else
                    let invPivot = recip pivot
                        lCol = [ if r < i then 0 else if r == i then 1 else (mat !! r) !! i * invPivot
                               | r <- [0 .. n - 1]
                               ]
                        -- Schur complement: A' = A - pivot * lCol * lCol^T
                        mat' = [ [ (mat !! r) !! c - pivot * (lCol !! r) * (lCol !! c)
                                 | c <- [0 .. n - 1]
                                 ]
                               | r <- [0 .. n - 1]
                               ]
                    in go (i + 1) mat' ((pivot, lCol) : acc)

-- | Checks whether a polynomial is a Sum of Squares (SOS).
isSOSPolynomial :: (Ord a, Fractional a) => Poly a -> Bool
isSOSPolynomial p =
  case decomposeSOS p of
    Just _  -> True
    Nothing -> False

-- | Computes an exact Sum of Squares decomposition of polynomial p:
-- p = sum c_k * q_k(x)^2 with c_k >= 0.
-- Returns Nothing if p is not a sum of squares.
decomposeSOS :: forall a. (Ord a, Fractional a) => Poly a -> Maybe (SOSDecomposition a)
decomposeSOS p@(Poly terms ord) =
  let deg = polyDegree p
  in if odd deg
     then Nothing -- odd degree polynomials cannot be non-negative everywhere
     else
       let halfDeg = deg `div` 2
           numVars = maximum (-1 : [ sparseMaxVar m | Term _ m <- terms ]) + 1
           basis = monomialBasisUpToDegree numVars halfDeg
           gram = gramMatrixForPoly p basis
       in case choleskyLDLT gram of
            Nothing -> Nothing
            Just ldltFactors ->
              let components = [ (c, makePoly lVec)
                               | (c, lVec) <- ldltFactors
                               , c > 0
                               ]
              in Just SOSDecomposition
                   { sosComponents = components
                   , sosGramMatrix = gram
                   , sosBasis      = basis
                   }
  where
    makePoly lVec =
      let basis = monomialBasisUpToDegree (maximum (-1 : [ sparseMaxVar m | Term _ m <- terms ]) + 1) (polyDegree p `div` 2)
          termList = [ Term coeff m | (coeff, m) <- zip lVec basis, coeff /= 0 ]
      in Poly termList ord

-- | The famous Motzkin polynomial (Theodore Motzkin 1967):
-- M(x, y) = x^4 * y^2 + x^2 * y^4 - 3 * x^2 * y^2 + 1.
-- Non-negative everywhere on R^2 (by AM-GM inequality), but NOT a sum of squares of polynomials!
motzkinPolynomial :: (Ord a, Fractional a) => Poly a
motzkinPolynomial =
  let x4y2 = sparseMul (sparseVar 0 4) (sparseVar 1 2)
      x2y4 = sparseMul (sparseVar 0 2) (sparseVar 1 4)
      x2y2 = sparseMul (sparseVar 0 2) (sparseVar 1 2)
      oneM = sparseOne
      terms = [ Term 1 x4y2, Term 1 x2y4, Term (-3) x2y2, Term 1 oneM ]
  in Poly terms DegRevLex
