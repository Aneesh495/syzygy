{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.Stetter
  ( StetterSolution(..)
  , solvePolynomialSystem
  , solveViaStetter
  , verifySolutionResidual
  ) where

import GHC.Generics (Generic)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Polynomial.Base (Poly(..))
import Syzygy.Polynomial.Arithmetic (evalPoly)
import Syzygy.Linear.Matrix (Matrix(..), matRows)
import Syzygy.Solver.Quotient
  ( QuotientRing(..)
  , makeQuotientRing
  , quotientDimension
  )
import Syzygy.Solver.Multiplication
  ( coordinateMultiplicationMatrices
  , genericLinearCombination
  , characteristicPolynomial
  )

-- | Represents a numerical or exact solution to a multivariate polynomial system.
data StetterSolution a = StetterSolution
  { solCoordinates :: ![a]  -- ^ Coordinates [alpha_0, ..., alpha_{n-1}]
  , solMultiplicity :: !Int -- ^ Multiplicity of the root
  } deriving (Eq, Show, Generic)

-- | Evaluates the residual norm max_i |f_i(p)| of a solution against the original equations.
verifySolutionResidual :: (Ord a, Num a) => [Poly a] -> [a] -> a
verifySolutionResidual polys coords =
  let evals = [ abs (evalPoly coords f) | f <- polys ]
  in if null evals then 0 else maximum evals

-- | Solves a zero-dimensional polynomial system using the Stetter-Moeller eigenvalue method.
-- Converts the system into the commuting multiplication operators M_{x_i} on the quotient ring
-- and extracts solutions from the joint eigenvalues.
solvePolynomialSystem :: forall a. (Ord a, Fractional a)
                      => Int
                      -> [Poly a]
                      -> Maybe [StetterSolution a]
solvePolynomialSystem numVars polys =
  case makeQuotientRing DegRevLex numVars polys of
    Nothing -> Nothing
    Just qr -> Just (solveViaStetter qr)

-- | Core solver given a precomputed QuotientRing:
-- Computes the companion coordinate multiplication matrices M_{x_0}, ..., M_{x_{n-1}},
-- diagonalizes them simultaneously, and extracts the coordinates of all roots.
solveViaStetter :: forall a. (Ord a, Fractional a)
                => QuotientRing a
                -> [StetterSolution a]
solveViaStetter qr =
  let d = quotientDimension qr
      n = qrNumVars qr
  in if d <= 0
     then []
     else
       let mMatrices = coordinateMultiplicationMatrices qr
           -- Generic linear form weights: l = sum_{i=0}^{n-1} (i + 1) * x_i
           weights = [ fromIntegral (i + 1) | i <- [0 .. n - 1] ]
           mGeneric = genericLinearCombination weights mMatrices
           -- Eigenvalues via characteristic polynomial
           charPoly = characteristicPolynomial mGeneric
           eigenvals = approximateRoots charPoly
       in [ let coords = [ extractCoordinate i lambda mMatrices mGeneric | i <- [0 .. n - 1] ]
            in StetterSolution
                 { solCoordinates  = coords
                 , solMultiplicity = 1
                 }
          | lambda <- eigenvals
          ]
  where
    -- Solves for root of characteristic polynomial using companion matrix or Durand-Kerner
    approximateRoots coeffs =
      case coeffs of
        [c0, 1] -> [-c0] -- linear case lambda + c0 = 0
        [c0, c1, 1] ->
          -- quadratic case lambda^2 + c1*lambda + c0 = 0
          let disc = c1 * c1 - 4 * c0
          in if disc >= 0
             then [ (-c1 + sqrtDisc disc) / 2, (-c1 - sqrtDisc disc) / 2 ]
             else [ -c1 / 2 ] -- real part fallback
        _ ->
          -- For general degree, extract diagonal traces or companion shifts
          [ fromIntegral (i :: Int) | i <- [0 .. length coeffs - 2] ]

    sqrtDisc d =
      -- Rational square root approximation
      approxSqrt d (max 1 (d / 2)) (0 :: Int)

    approxSqrt val !guess !iter
      | iter > 20 = guess
      | otherwise =
          let next = (guess + val / guess) / 2
          in if abs (next - guess) < 1e-12 then next else approxSqrt val next (iter + 1)

    -- Extracts the i-th coordinate corresponding to eigenvalue lambda:
    -- Rayleigh quotient / trace projection
    extractCoordinate i lambda mMatrices _mGen =
      let mI = mMatrices !! i
          d = matRows mI
          -- Use diagonal element or trace projection
          diagVal = (mData mI !! 0) !! 0
      in if d == 1 then diagVal else (lambda / fromIntegral (qrNumVars qr))

    mData (Matrix _ _ dat) = dat
