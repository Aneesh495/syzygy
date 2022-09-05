{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Syzygy.Core.Dual
  ( Dual(..)
  , HyperDual(..)
  , dual
  , hyperDual
  , primal
  , tangent
  , diff
  , diff2
  , grad
  , hessian
  , directionalDeriv
  ) where

import GHC.Generics (Generic)

-- | First-order dual number: a + b * epsilon where epsilon^2 = 0.
-- Enables exact evaluation of first-order derivatives without discretization errors.
data Dual a = Dual !a !a
  deriving (Eq, Ord, Generic)

-- | Smart constructor for a dual number.
dual :: a -> a -> Dual a
dual p t = Dual p t

-- | Extracts the primal value f(x).
primal :: Dual a -> a
primal (Dual p _) = p

-- | Extracts the tangent derivative value f'(x).
tangent :: Dual a -> a
tangent (Dual _ t) = t

-- | Computes the exact first derivative of a function at a point x.
diff :: Num a => (Dual a -> Dual a) -> a -> a
diff f x = tangent (f (Dual x 1))

-- | Computes the gradient vector of a multivariate function at a point xs.
grad :: Num a => ([Dual a] -> Dual a) -> [a] -> [a]
grad f xs =
  [ tangent (f (perturb i xs)) | i <- [0 .. length xs - 1] ]
  where
    perturb targetIdx list =
      [ if idx == targetIdx then Dual val 1 else Dual val 0
      | (idx, val) <- zip [0 :: Int ..] list
      ]

-- | Computes the directional derivative of a multivariate function along direction v at point xs.
directionalDeriv :: Num a => ([Dual a] -> Dual a) -> [a] -> [a] -> a
directionalDeriv f xs vs =
  tangent (f (zipWith Dual xs vs))

-- | Second-order hyper-dual number:
-- f + f1*eps1 + f2*eps2 + f12*eps1*eps2 where eps1^2 = eps2^2 = 0.
-- Computes exact first and second derivatives simultaneously.
data HyperDual a = HyperDual !a !a !a !a
  deriving (Eq, Ord, Generic)

-- | Smart constructor for hyper-dual number.
hyperDual :: a -> a -> a -> a -> HyperDual a
hyperDual f0 f1 f2 f12 = HyperDual f0 f1 f2 f12

-- | Computes value, first derivative, and second derivative at point x.
-- Returns (f(x), f'(x), f''(x)).
diff2 :: Num a => (HyperDual a -> HyperDual a) -> a -> (a, a, a)
diff2 f x =
  let HyperDual f0 f1 _ f12 = f (HyperDual x 1 1 0)
  in (f0, f1, f12)

-- | Computes the exact Hessian matrix of a multivariate function at point xs.
hessian :: Num a => ([HyperDual a] -> HyperDual a) -> [a] -> [[a]]
hessian f xs =
  let n = length xs
  in [ [ evaluateSecondCross i j | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]
  where
    evaluateSecondCross i j =
      let inputs =
            [ if idx == i && idx == j
              then HyperDual x 1 1 0
              else if idx == i
              then HyperDual x 1 0 0
              else if idx == j
              then HyperDual x 0 1 0
              else HyperDual x 0 0 0
            | (idx, x) <- zip [0 :: Int ..] xs
            ]
          HyperDual _ _ _ f12 = f inputs
      in f12

instance Show a => Show (Dual a) where
  show (Dual p t) = show p ++ " + " ++ show t ++ "*eps"

instance Show a => Show (HyperDual a) where
  show (HyperDual f0 f1 f2 f12) =
    show f0 ++ " + " ++ show f1 ++ "*eps1 + " ++ show f2 ++ "*eps2 + " ++ show f12 ++ "*eps1*eps2"

instance Num a => Num (Dual a) where
  (Dual p1 t1) + (Dual p2 t2) = Dual (p1 + p2) (t1 + t2)
  (Dual p1 t1) - (Dual p2 t2) = Dual (p1 - p2) (t1 - t2)
  (Dual p1 t1) * (Dual p2 t2) = Dual (p1 * p2) (p1 * t2 + t1 * p2)
  negate (Dual p t) = Dual (-p) (-t)
  abs (Dual p t) = Dual (abs p) (t * signum p)
  signum (Dual p _) = Dual (signum p) 0
  fromInteger n = Dual (fromInteger n) 0

instance Fractional a => Fractional (Dual a) where
  fromRational r = Dual (fromRational r) 0
  (Dual p1 t1) / (Dual p2 t2) =
    Dual (p1 / p2) ((t1 * p2 - p1 * t2) / (p2 * p2))
  recip (Dual p t) =
    Dual (recip p) (-t / (p * p))

instance Floating a => Floating (Dual a) where
  pi = Dual pi 0
  exp (Dual p t) = let e = exp p in Dual e (t * e)
  log (Dual p t) = Dual (log p) (t / p)
  sin (Dual p t) = Dual (sin p) (t * cos p)
  cos (Dual p t) = Dual (cos p) (-t * sin p)
  asin (Dual p t) = Dual (asin p) (t / sqrt (1 - p * p))
  acos (Dual p t) = Dual (acos p) (-t / sqrt (1 - p * p))
  atan (Dual p t) = Dual (atan p) (t / (1 + p * p))
  sinh (Dual p t) = Dual (sinh p) (t * cosh p)
  cosh (Dual p t) = Dual (cosh p) (t * sinh p)
  asinh (Dual p t) = Dual (asinh p) (t / sqrt (1 + p * p))
  acosh (Dual p t) = Dual (acosh p) (t / sqrt (p * p - 1))
  atanh (Dual p t) = Dual (atanh p) (t / (1 - p * p))

instance Num a => Num (HyperDual a) where
  (HyperDual a0 a1 a2 a12) + (HyperDual b0 b1 b2 b12) =
    HyperDual (a0 + b0) (a1 + b1) (a2 + b2) (a12 + b12)

  (HyperDual a0 a1 a2 a12) - (HyperDual b0 b1 b2 b12) =
    HyperDual (a0 - b0) (a1 - b1) (a2 - b2) (a12 - b12)

  (HyperDual a0 a1 a2 a12) * (HyperDual b0 b1 b2 b12) =
    let c0 = a0 * b0
        c1 = a0 * b1 + a1 * b0
        c2 = a0 * b2 + a2 * b0
        c12 = a0 * b12 + a1 * b2 + a2 * b1 + a12 * b0
    in HyperDual c0 c1 c2 c12

  negate (HyperDual a0 a1 a2 a12) =
    HyperDual (-a0) (-a1) (-a2) (-a12)

  abs (HyperDual a0 a1 a2 a12) =
    let s = signum a0
    in HyperDual (abs a0) (a1 * s) (a2 * s) (a12 * s)

  signum (HyperDual a0 _ _ _) =
    HyperDual (signum a0) 0 0 0

  fromInteger n = HyperDual (fromInteger n) 0 0 0

instance Fractional a => Fractional (HyperDual a) where
  fromRational r = HyperDual (fromRational r) 0 0 0

  (HyperDual a0 a1 a2 a12) / (HyperDual b0 b1 b2 b12) =
    let invB0 = 1 / b0
        invB02 = invB0 * invB0
        c0 = a0 * invB0
        c1 = (a1 * b0 - a0 * b1) * invB02
        c2 = (a2 * b0 - a0 * b2) * invB02
        c12 = (a12 * b0 - a0 * b12 - a1 * b2 - a2 * b1 + 2 * a0 * b1 * b2 * invB0) * invB02
    in HyperDual c0 c1 c2 c12

  recip (HyperDual b0 b1 b2 b12) =
    let invB0 = 1 / b0
        invB02 = invB0 * invB0
        c0 = invB0
        c1 = -b1 * invB02
        c2 = -b2 * invB02
        c12 = (-b12 + 2 * b1 * b2 * invB0) * invB02
    in HyperDual c0 c1 c2 c12
