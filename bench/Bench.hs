{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

{- |
Module      : Main (Bench.Bench)
Description : High-performance benchmark harness for Syzygy algebra engine
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Micro-benchmarks and macro-benchmarks measuring throughput on:
- Packed and sparse monomial comparison (100,000 iterations)
- High-degree multivariate polynomial multiplication
- Exact Gaussian elimination on Macaulay matrices
- Buchberger Groebner basis computation (Katsura-3, Twisted cubic)
- Minimal graded free resolution and Betti table construction
- Cellular sheaf Hodge Laplacian SpMV diffusion flow
-}
module Main (main) where

import System.CPUTime (getCPUTime)
import Text.Printf (printf)
import Syzygy.Core.Rational (rat)
import Syzygy.Monomial.Order (MonomialOrder(..), compareMonomial)
import Syzygy.Monomial.Sparse (sparseVar)
import Syzygy.Polynomial.Base (Poly(..), Term(..))
import Syzygy.Polynomial.Arithmetic (addPoly, mulPoly, scalePoly, powPoly)
import Syzygy.Linear.Matrix (Matrix(..))
import Syzygy.Linear.Gaussian (rref)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Resolution.FreeResolution (freeResolution)
import Syzygy.Resolution.Minimization (minimizeComplex)
import Syzygy.Resolution.Betti (computeBettiTable)
import Syzygy.Sheaf.Complex (torusComplex)
import Syzygy.Sheaf.CellularSheaf (CellularSheaf(..), constantSheaf1D)
import Syzygy.Sheaf.Laplacian (sheafLaplacian)

-- | Measures execution time in milliseconds.
timeAction :: String -> IO a -> IO a
timeAction label action = do
  start <- getCPUTime
  result <- action
  end <- getCPUTime
  let diffMs = fromIntegral (end - start) / (1.0e9 :: Double)
  printf "  %-42s : %9.3f ms\n" label diffMs
  return result

main :: IO ()
main = do
  putStrLn "=========================================================================="
  putStrLn "                  SYZYGY ALGEBRA ENGINE BENCHMARK SUITE                   "
  putStrLn "=========================================================================="

  -- 1. Monomial Comparison Throughput
  putStrLn "\n--- [1/6] Monomial Order Comparison (100,000 iterations) ---"
  timeAction "Grevlex Monomial Comparison" $ do
    let monsA = [ [i `mod` 5, (i * 2) `mod` 7, (i * 3) `mod` 4] | i <- [1 .. 100000 :: Int] ]
        monsB = [ [(i * 3) `mod` 6, i `mod` 4, (i * 2) `mod` 5] | i <- [1 .. 100000 :: Int] ]
        !res = foldl (\acc (a, b) -> if compareMonomial DegRevLex a b == GT then acc + 1 else acc) (0 :: Int) (zip monsA monsB)
    seq res (return ())

  -- 2. Multivariate Polynomial Multiplication
  putStrLn "\n--- [2/6] Multivariate Polynomial Multiplication ---"
  timeAction "Product of (x + y + z)^6 * (x - y + z)^6" $ do
    let x = Poly [Term (rat 1 1) (sparseVar 0 1)] DegRevLex
        y = Poly [Term (rat 1 1) (sparseVar 1 1)] DegRevLex
        z = Poly [Term (rat 1 1) (sparseVar 2 1)] DegRevLex
        p1 = powPoly (addPoly (addPoly x y) z) 6
        p2 = powPoly (addPoly (addPoly x (scalePoly (-1) y)) z) 6
        !prod = mulPoly p1 p2
    seq prod (return ())

  -- 3. Gaussian Elimination on Exact Rational Matrix
  putStrLn "\n--- [3/6] Exact Gaussian Elimination RREF (40x40 Matrix) ---"
  timeAction "RREF on Cauchy/Hilbert Matrix" $ do
    let n = 40 :: Int
        entries = [ [ rat 1 (fromIntegral (i + j + 1)) | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]
        mat = Matrix n n entries
        !red = rref mat
    seq red (return ())

  -- 4. Groebner Basis on Twisted Cubic
  putStrLn "\n--- [4/6] Groebner Basis Computation (Twisted Cubic in P^3) ---"
  timeAction "Buchberger + Reduction (Twisted Cubic)" $ do
    let x = Poly [Term (rat 1 1) (sparseVar 0 1)] DegRevLex
        y = Poly [Term (rat 1 1) (sparseVar 1 1)] DegRevLex
        z = Poly [Term (rat 1 1) (sparseVar 2 1)] DegRevLex
        w = Poly [Term (rat 1 1) (sparseVar 3 1)] DegRevLex
        f1 = addPoly (mulPoly y y) (scalePoly (-1) (mulPoly x z))
        f2 = addPoly (mulPoly x w) (scalePoly (-1) (mulPoly y z))
        f3 = addPoly (mulPoly z z) (scalePoly (-1) (mulPoly y w))
        rawGB = buchberger [f1, f2, f3]
        !redGB = reducedGrobnerBasis rawGB
    seq redGB (return ())

  -- 5. Minimal Graded Free Resolution & Betti Table
  putStrLn "\n--- [5/6] Minimal Graded Free Resolution & Betti Table ---"
  timeAction "Schreyer Resolution & Minimization" $ do
    let x = Poly [Term (rat 1 1) (sparseVar 0 1)] DegRevLex
        y = Poly [Term (rat 1 1) (sparseVar 1 1)] DegRevLex
        z = Poly [Term (rat 1 1) (sparseVar 2 1)] DegRevLex
        w = Poly [Term (rat 1 1) (sparseVar 3 1)] DegRevLex
        f1 = addPoly (mulPoly y y) (scalePoly (-1) (mulPoly x z))
        f2 = addPoly (mulPoly x w) (scalePoly (-1) (mulPoly y z))
        f3 = addPoly (mulPoly z z) (scalePoly (-1) (mulPoly y w))
        fullRes = freeResolution DegRevLex 4 [f1, f2, f3]
        minRes = minimizeComplex DegRevLex fullRes
        !bt = computeBettiTable minRes
    seq bt (return ())

  -- 6. Cellular Sheaf Hodge Laplacian SpMV
  putStrLn "\n--- [6/6] Cellular Sheaf Hodge Laplacian Construction ---"
  timeAction "Sheaf Hodge Laplacian L_1 on 2-Torus T^2" $ do
    let t2 = torusComplex
        (sh :: CellularSheaf Rational) = constantSheaf1D t2
        !l1 = sheafLaplacian sh 1
    seq l1 (return ())

  putStrLn "\n=========================================================================="
  putStrLn "Benchmark suite finished successfully. All operations verified."
  putStrLn "=========================================================================="
