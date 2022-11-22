{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.Multiplication
  ( coordinateMultiplicationMatrices
  , verifyCommutingOperators
  , genericLinearCombination
  , companionMatrix
  , characteristicPolynomial
  ) where

import Syzygy.Monomial.Sparse (sparseVar)
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  )
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matRows
  , matMul
  , matScale
  , matAdd
  )
import Syzygy.Solver.Quotient
  ( QuotientRing(..)
  , multiplicationMatrix
  )

-- | Computes the multiplication matrices M_{x_0}, M_{x_1}, ..., M_{x_{n-1}}
-- for all coordinate variables in the quotient ring A = k[x_0, ..., x_{n-1}] / I.
coordinateMultiplicationMatrices :: (Ord a, Fractional a) => QuotientRing a -> [Matrix a]
coordinateMultiplicationMatrices qr =
  let n = qrNumVars qr
      ord = qrOrder qr
  in [ let xVar = Poly [Term 1 (sparseVar v 1)] ord
       in multiplicationMatrix qr xVar
     | v <- [0 .. n - 1]
     ]

-- | Verifies that all coordinate multiplication operators mutually commute:
-- M_{x_i} * M_{x_j} == M_{x_j} * M_{x_i} for all i < j.
-- By algebra theory, in the quotient ring A = k[x]/I, multiplication is commutative,
-- so this identity is an invariant test for correct Groebner basis and quotient arithmetic.
verifyCommutingOperators :: (Eq a, Num a) => [Matrix a] -> Bool
verifyCommutingOperators mats =
  and [ let prod1 = matMul (mats !! i) (mats !! j)
            prod2 = matMul (mats !! j) (mats !! i)
        in prod1 == prod2
      | i <- [0 .. length mats - 1]
      , j <- [i + 1 .. length mats - 1]
      ]

-- | Computes a generic linear combination of multiplication operators:
-- M_l = sum_{i=0}^{n-1} c_i * M_{x_i}.
-- For generic coefficients c, the eigenvalues of M_l are all distinct,
-- separating the roots of the polynomial system (Shape Lemma).
genericLinearCombination :: (Eq a, Num a) => [a] -> [Matrix a] -> Matrix a
genericLinearCombination weights mats =
  case (weights, mats) of
    ([], []) -> Matrix 0 0 []
    (w:ws, m:ms) ->
      foldl matAdd (matScale w m) [ matScale wi mi | (wi, mi) <- zip ws ms ]
    _ -> Matrix 0 0 []

-- | Constructs the classic companion matrix of a monic univariate polynomial
-- p(x) = x^n + c_{n-1} x^{n-1} + ... + c_1 x + c_0:
-- [ 0  0  ...  -c_0     ]
-- [ 1  0  ...  -c_1     ]
-- [ 0  1  ...  -c_2     ]
-- [ ...        ...      ]
-- [ 0  0  ... 1 -c_{n-1} ]
companionMatrix :: forall a. (Eq a, Num a) => [a] -> Matrix a
companionMatrix coeffs =
  -- coeffs = [c_0, c_1, ..., c_{n-1}] for monic degree n
  let n = length coeffs
  in if n == 0
     then Matrix 0 0 []
     else
       let row i = [ if j == n - 1
                     then -(coeffs !! i)
                     else if i == j + 1
                          then 1
                          else 0
                   | j <- [0 .. n - 1]
                   ]
       in Matrix n n [ row i | i <- [0 .. n - 1] ]

-- | Computes the characteristic polynomial det(lambda * I - M) of a matrix
-- using the Faddeev-LeVerrier algorithm over any field of characteristic zero.
-- Returns the coefficients [c_0, c_1, ..., c_n] where c_n = 1.
characteristicPolynomial :: forall a. (Eq a, Fractional a) => Matrix a -> [a]
characteristicPolynomial mat =
  let n = matRows mat
  in if n == 0
     then [1]
     else
       let ident = matIdentity n
           -- Faddeev-LeVerrier iteration
           -- M_k = M * (M_{k-1} - c_{n-k+1} * I)
           -- c_{n-k} = - (1 / k) * Trace(M_k)
           go !k !accMat coeffs
             | k > n = reverse coeffs
             | otherwise =
                 let prevC = case coeffs of (c:_) -> c; [] -> 1
                     mAdj = if k == 1 then mat else matMul mat (matAdd accMat (matScale (-prevC) ident))
                     tr = matTrace mAdj
                     ck = -tr / fromIntegral k
                 in go (k + 1) mAdj (ck : coeffs)
       in go 1 mat [1]
  where
    matIdentity dim =
      Matrix dim dim [ [ if i == j then 1 else 0 | j <- [0 .. dim - 1] ] | i <- [0 .. dim - 1] ]
    matTrace (Matrix r _ d) = sum [ (d !! i) !! i | i <- [0 .. r - 1] ]
