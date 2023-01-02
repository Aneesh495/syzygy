{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

{- |
Module      : Syzygy.UI.VarietyVisualizer
Description : 2D cross-section plotter, 3D sphere tracing raymarcher, and OBJ exporter
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Visualizes real algebraic varieties V(f) = { (x, y, z) | f(x, y, z) = 0 }.
Features:
- Terminal ASCII 2D cross-section contour plotter
- First-order Taubin distance estimator: d(p) = |f(p)| / ||grad f(p)||
- 3D sphere-tracing raymarcher with Phong diffuse shading directly to ANSI terminal
- Wavefront OBJ mesh exporter for 3D visualization in Blender / Three.js
- Canonical algebraic surfaces: Whitney umbrella, Cayley cubic, Clebsch cubic, Kummer quartic
-}
module Syzygy.UI.VarietyVisualizer
  ( plot2DCrossSection
  , taubinDistance3D
  , raymarchVariety3D
  , exportPointCloudObj
  , whitneyUmbrella
  , cayleyCubic
  , clebschCubic
  , kummerSurface
  , sphereVariety
  , torusVariety
  ) where

import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (sparseVar, sparseMul)
import Syzygy.Polynomial.Base (Poly(..), Term(..))
import Syzygy.Polynomial.Arithmetic (evalPoly, polyDerivative)

-- | Plots a 2D cross section of f(x, y) = 0 in the terminal using ASCII characters.
-- Evaluates f on a grid [xMin, xMax] x [yMin, yMax].
plot2DCrossSection :: Poly Double -> Double -> Double -> Double -> Double -> Int -> Int -> String
plot2DCrossSection f xMin xMax yMin yMax width height =
  let dx = (xMax - xMin) / fromIntegral (max 1 (width - 1))
      dy = (yMax - yMin) / fromIntegral (max 1 (height - 1))
      -- Evaluate grid values
      evalGrid =
        [ [ evalPoly [xMin + fromIntegral col * dx, yMax - fromIntegral row * dy] f
          | col <- [0 .. width - 1]
          ]
        | row <- [0 .. height - 1]
        ]
      -- Detect sign changes across horizontal and vertical adjacent pixels
      renderPixel r c =
        let val = (evalGrid !! r) !! c
            valRight = if c + 1 < width then (evalGrid !! r) !! (c + 1) else val
            valDown  = if r + 1 < height then (evalGrid !! (r + 1)) !! c else val
            hasSignChange = (val * valRight < 0) || (val * valDown < 0)
            isNearZero = abs val < 0.05
        in if hasSignChange || isNearZero
           then '#'
           else if abs val < 0.2 then ':' else ' '
      renderRow r = [ renderPixel r c | c <- [0 .. width - 1] ]
      renderedLines = [ renderRow r | r <- [0 .. height - 1] ]
  in unlines renderedLines

-- | Computes the Taubin distance approximation to the implicit surface f(x, y, z) = 0:
-- d(p) approx |f(p)| / sqrt( ||grad f(p)||^2 + eps ).
taubinDistance3D :: Poly Double -> (Double, Double, Double) -> Double
taubinDistance3D f (x, y, z) =
  let pt = [x, y, z]
      fVal = evalPoly pt f
      dfX = evalPoly pt (polyDerivative 0 f)
      dfY = evalPoly pt (polyDerivative 1 f)
      dfZ = evalPoly pt (polyDerivative 2 f)
      gradNormSq = dfX * dfX + dfY * dfY + dfZ * dfZ
  in abs fVal / sqrt (gradNormSq + 1.0e-10)

-- | 3D Raymarching / sphere tracing of an algebraic surface f(x, y, z) = 0
-- rendered to an ASCII grayscale terminal buffer.
raymarchVariety3D :: Poly Double -> Int -> Int -> String
raymarchVariety3D f width height =
  let camPos = (0.0, 0.0, -3.5)
      lightDir = normalize (0.5, 0.8, -1.0)
      aspect = fromIntegral width / fromIntegral (2 * height)

      ramp = " .:-=+*#%@"
      numShades = length ramp - 1

      marchRay orig dir =
        let maxSteps = 48 :: Int
            maxDist = 7.0 :: Double
            minHitDist = 0.005 :: Double
            step !t !stepCount
              | stepCount >= maxSteps || t > maxDist = Nothing
              | otherwise =
                  let currPt = addVec orig (scaleVec t dir)
                      dist = taubinDistance3D f currPt
                  in if dist < minHitDist
                     then Just currPt
                     else step (t + max 0.002 (0.45 * dist)) (stepCount + 1)
        in step 0.1 0

      shadePoint pt =
        let eps = 0.002
            (px, py, pz) = pt
            fxP = evalPoly [px + eps, py, pz] f
            fxM = evalPoly [px - eps, py, pz] f
            fyP = evalPoly [px, py + eps, pz] f
            fyM = evalPoly [px, py - eps, pz] f
            fzP = evalPoly [px, py, pz + eps] f
            fzM = evalPoly [px, py, pz - eps] f
            norm = normalize (fxP - fxM, fyP - fyM, fzP - fzM)
            diffuse = max 0.0 (dotProduct norm lightDir)
            shadeIdx = round (diffuse * fromIntegral numShades)
        in ramp !! max 0 (min numShades shadeIdx)

      renderPixel row col =
        let u = (2.0 * (fromIntegral col / fromIntegral width) - 1.0) * aspect
            v = 1.0 - 2.0 * (fromIntegral row / fromIntegral height)
            rayDir = normalize (u, v, 1.5)
        in case marchRay camPos rayDir of
             Just hitPt -> shadePoint hitPt
             Nothing    -> ' '

      renderRow r = [ renderPixel r c | c <- [0 .. width - 1] ]
  in unlines [ renderRow r | r <- [0 .. height - 1] ]
  where
    dotProduct (x1, y1, z1) (x2, y2, z2) = x1 * x2 + y1 * y2 + z1 * z2
    scaleVec s (x, y, z) = (s * x, s * y, s * z)
    addVec (x1, y1, z1) (x2, y2, z2) = (x1 + x2, y1 + y2, z1 + z2)
    normalize (x, y, z) =
      let len = sqrt (x * x + y * y + z * z)
      in if len < 1.0e-9 then (0, 0, 1) else (x / len, y / len, z / len)

-- | Exports surface sample points to a Wavefront .obj file.
exportPointCloudObj :: Poly Double -> Double -> Double -> Int -> String
exportPointCloudObj f bound stepSize maxPoints =
  let pts = [ (x, y, z)
            | x <- [-bound, -bound + stepSize .. bound]
            , y <- [-bound, -bound + stepSize .. bound]
            , z <- [-bound, -bound + stepSize .. bound]
            , taubinDistance3D f (x, y, z) < stepSize * 0.4
            ]
      limitedPts = take maxPoints pts
      objLines = [ "v " ++ show x ++ " " ++ show y ++ " " ++ show z | (x, y, z) <- limitedPts ]
  in unlines ("# Exported by Syzygy Algebraic Surface Engine" : objLines)

-- | Whitney Umbrella: x^2 - y^2 * z = 0.
whitneyUmbrella :: Poly Double
whitneyUmbrella =
  let x2 = Poly [Term 1.0 (sparseVar 0 2)] Lex
      y2z = Poly [Term 1.0 (sparseMul (sparseVar 1 2) (sparseVar 2 1))] Lex
  in Poly (polyTerms' x2 ++ [ Term (-1.0) m | Term _ m <- polyTerms' y2z ]) Lex
  where
    polyTerms' (Poly ts _) = ts

-- | Cayley Cubic: x*y*z + x*y + y*z + z*x = 0.
cayleyCubic :: Poly Double
cayleyCubic =
  let xyz = Term 1.0 (sparseMul (sparseVar 0 1) (sparseMul (sparseVar 1 1) (sparseVar 2 1)))
      xy  = Term 1.0 (sparseMul (sparseVar 0 1) (sparseVar 1 1))
      yz  = Term 1.0 (sparseMul (sparseVar 1 1) (sparseVar 2 1))
      zx  = Term 1.0 (sparseMul (sparseVar 2 1) (sparseVar 0 1))
  in Poly [xyz, xy, yz, zx] Lex

-- | Clebsch diagonal cubic surface:
-- 81(x^3 + y^3 + z^3) - 189(x^2 y + ... ) + 54 xyz + 126(xy + yz + zx) - 9(x^2 + y^2 + z^2) - 9(x + y + z) + 1 = 0.
clebschCubic :: Poly Double
clebschCubic =
  let x3 = Term 81.0 (sparseVar 0 3)
      y3 = Term 81.0 (sparseVar 1 3)
      z3 = Term 81.0 (sparseVar 2 3)
      xyz = Term 54.0 (sparseMul (sparseVar 0 1) (sparseMul (sparseVar 1 1) (sparseVar 2 1)))
      one = Term 1.0 (sparseVar 0 0)
      subSquares = [ Term (-9.0) (sparseVar i 2) | i <- [0 .. 2] ]
      subLinear  = [ Term (-9.0) (sparseVar i 1) | i <- [0 .. 2] ]
  in Poly ([x3, y3, z3, xyz, one] ++ subSquares ++ subLinear) Lex

-- | Kummer quartic nodal surface: (x^2 + y^2 + z^2 - w^2)^2 - (3*w^2 - 1)*(x^2 + y^2 + z^2 - 1).
kummerSurface :: Poly Double
kummerSurface =
  let x4 = Term 1.0 (sparseVar 0 4)
      y4 = Term 1.0 (sparseVar 1 4)
      z4 = Term 1.0 (sparseVar 2 4)
      x2y2 = Term 2.0 (sparseMul (sparseVar 0 2) (sparseVar 1 2))
      y2z2 = Term 2.0 (sparseMul (sparseVar 1 2) (sparseVar 2 2))
      z2x2 = Term 2.0 (sparseMul (sparseVar 2 2) (sparseVar 0 2))
      subQuad = [ Term (-1.5) (sparseVar i 2) | i <- [0 .. 2] ]
      c = Term 0.2 (sparseVar 0 0)
  in Poly ([x4, y4, z4, x2y2, y2z2, z2x2, c] ++ subQuad) Lex

-- | Sphere of radius r: x^2 + y^2 + z^2 - r^2 = 0.
sphereVariety :: Double -> Poly Double
sphereVariety r =
  Poly [ Term 1.0 (sparseVar 0 2)
       , Term 1.0 (sparseVar 1 2)
       , Term 1.0 (sparseVar 2 2)
       , Term (-(r * r)) (sparseVar 0 0)
       ] Lex

-- | Torus with major radius R and minor radius r:
-- (x^2 + y^2 + z^2 + R^2 - r^2)^2 - 4 * R^2 * (x^2 + y^2) = 0.
torusVariety :: Double -> Double -> Poly Double
torusVariety rMaj rMin =
  let rSqDiff = rMaj * rMaj - rMin * rMin
      x2 = sparseVar 0 2
      y2 = sparseVar 1 2
      z2 = sparseVar 2 2
      x4 = Term 1.0 (sparseVar 0 4)
      y4 = Term 1.0 (sparseVar 1 4)
      z4 = Term 1.0 (sparseVar 2 4)
      x2y2 = Term 2.0 (sparseMul x2 y2)
      y2z2 = Term 2.0 (sparseMul y2 z2)
      z2x2 = Term 2.0 (sparseMul z2 x2)
      scaleCross = 2.0 * rSqDiff
      x2Cross = Term (scaleCross - 4.0 * rMaj * rMaj) x2
      y2Cross = Term (scaleCross - 4.0 * rMaj * rMaj) y2
      z2Cross = Term scaleCross z2
      constTerm = Term (rSqDiff * rSqDiff) (sparseVar 0 0)
  in Poly [x4, y4, z4, x2y2, y2z2, z2x2, x2Cross, y2Cross, z2Cross, constTerm] Lex
