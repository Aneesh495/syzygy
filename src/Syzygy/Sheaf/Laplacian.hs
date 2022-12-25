{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Sheaf.Laplacian
Description : Sheaf Hodge Laplacians, harmonic cochains, Hodge decomposition, and heat flow
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

The cellular sheaf Hodge Laplacian linearizes sheaf cohomology through
spectral graph theory and differential forms (Hansen and Ghrist 2019).

For a cellular sheaf F over a cell complex X, the k-th Hodge Laplacian is:
  L_k = delta_{k-1} . delta_{k-1}^* + delta_k^* . delta_k = L_k^{down} + L_k^{up}

By the Hodge Decomposition Theorem for cellular sheaves, the k-cochains split orthogonally:
  C^k(X; F) = im(delta_{k-1}) direct_sum ker(L_k) direct_sum im(delta_k^*)

where the space of harmonic cochains ker(L_k) is canonically isomorphic to
the k-th sheaf cohomology group:
  ker(L_k) =~ H^k(X; F)

In particular, dim ker(L_k) = b_k(X; F). For k = 0, L_0 = delta_0^* . delta_0
is the sheaf graph Laplacian whose kernel comprises the global sections.
-}
module Syzygy.Sheaf.Laplacian
  ( sheafLaplacianDown
  , sheafLaplacianUp
  , sheafLaplacian
  , sheafGraphLaplacian
  , harmonicCochains
  , verifyHodgeTheorem
  , sheafDiffusionStep
  , sheafHarmonicProjection
  ) where

import Syzygy.Linear.Vector (Vector(..))
import Syzygy.Linear.Matrix
  ( Matrix(..)
  , matTranspose
  , matMul
  , matAdd
  , matIdentity
  , matZero
  , matVecMul
  )
import Syzygy.Linear.Gaussian (nullSpace)
import Syzygy.Sheaf.Complex (complexDimension)
import Syzygy.Sheaf.CellularSheaf (CellularSheaf(..))
import Syzygy.Sheaf.Cohomology
  ( cochainDimension
  , coboundaryMatrix
  , sheafBettiNumber
  )

-- | Computes the down-Laplacian L_k^{down} = delta_{k-1} . delta_{k-1}^* on C^k(X; F).
sheafLaplacianDown :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Matrix a
sheafLaplacianDown sh k =
  let dimCk = cochainDimension sh k
  in if k <= 0 || dimCk == 0
     then matZero dimCk dimCk
     else
       let dkm1 = coboundaryMatrix sh (k - 1)
           dkm1Adj = matTranspose dkm1
       in matMul dkm1 dkm1Adj

-- | Computes the up-Laplacian L_k^{up} = delta_k^* . delta_k on C^k(X; F).
sheafLaplacianUp :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Matrix a
sheafLaplacianUp sh k =
  let dimCk = cochainDimension sh k
      maxDim = complexDimension (csComplex sh)
  in if k >= maxDim || dimCk == 0
     then matZero dimCk dimCk
     else
       let dk = coboundaryMatrix sh k
           dkAdj = matTranspose dk
       in matMul dkAdj dk

-- | Computes the total k-th Sheaf Hodge Laplacian L_k = L_k^{down} + L_k^{up}.
sheafLaplacian :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Matrix a
sheafLaplacian sh k =
  let down = sheafLaplacianDown sh k
      up   = sheafLaplacianUp sh k
  in matAdd down up

-- | Computes the Sheaf Graph Laplacian L_0 = delta_0^* . delta_0 on 0-cochains.
sheafGraphLaplacian :: (Eq a, Fractional a) => CellularSheaf a -> Matrix a
sheafGraphLaplacian sh = sheafLaplacian sh 0

-- | Computes a basis of harmonic cochains ker(L_k) =~ H^k(X; F).
harmonicCochains :: (Eq a, Fractional a) => CellularSheaf a -> Int -> [Vector a]
harmonicCochains sh k =
  let lk = sheafLaplacian sh k
  in nullSpace lk

-- | Verifies the Hodge Theorem: dim ker(L_k) == b_k(X; F).
verifyHodgeTheorem :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Bool
verifyHodgeTheorem sh k =
  let harms = harmonicCochains sh k
      betti = sheafBettiNumber sh k
  in length harms == betti

-- | Performs one step of discrete sheaf heat diffusion on a cochain:
-- x^{(m+1)} = (I - gamma * L_k) * x^{(m)}
sheafDiffusionStep :: forall a. (Eq a, Fractional a) => CellularSheaf a -> Int -> a -> Vector a -> Vector a
sheafDiffusionStep sh k gamma cochain =
  let dimCk = cochainDimension sh k
      lk = sheafLaplacian sh k
      ident = matIdentity dimCk
      scaledLk = case lk of
        Matrix r c entries -> Matrix r c [ [ -gamma * val | val <- row ] | row <- entries ]
      stepMat = matAdd ident scaledLk
  in matVecMul stepMat cochain

-- | Iteratively diffuses a cochain over numSteps to project towards the harmonic subspace:
-- x_infinity in ker(L_k) =~ H^k(X; F).
sheafHarmonicProjection :: (Eq a, Fractional a) => CellularSheaf a -> Int -> Int -> a -> Vector a -> Vector a
sheafHarmonicProjection sh k numSteps gamma initCochain =
  let go !step !c
        | step >= numSteps = c
        | otherwise = go (step + 1) (sheafDiffusionStep sh k gamma c)
  in go 0 initCochain
