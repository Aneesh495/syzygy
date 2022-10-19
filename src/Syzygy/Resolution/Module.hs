{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.Module
  ( ModMonomial(..)
  , ModTerm(..)
  , ModPoly(..)
  , zeroModPoly
  , isZeroModPoly
  , modPolyFromComponents
  , modPolyToComponents
  , addModPoly
  , subModPoly
  , scaleModPoly
  , mulPolyModPoly
  , modPolyDegree
  , modLeadingTerm
  , modLeadingMonomial
  , modLeadingCoeff
  ) where

import GHC.Generics (Generic)
import qualified Data.Map.Strict as Map
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseDegree
  , sparseMul
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , zeroPoly
  )

-- | A module monomial m * e_i representing monomial m in component e_i of a free module R^m.
data ModMonomial = ModMonomial
  { modMonoExp  :: !SparseMonomial
  , modMonoComp :: !Int
  } deriving (Eq, Ord, Show, Generic)

-- | A single module term c * m * e_i where c /= 0.
data ModTerm c = ModTerm
  { modTermCoeff :: !c
  , modTermExp   :: !SparseMonomial
  , modTermComp  :: !Int
  } deriving (Eq, Show, Generic)

-- | An element of a free module R^m over the polynomial ring R.
-- Represented canonically as a list of non-zero terms.
newtype ModPoly c = ModPoly
  { modPolyTerms :: [ModTerm c]
  } deriving (Eq, Show, Generic)

-- | The zero element in the module (empty list of terms).
zeroModPoly :: ModPoly c
zeroModPoly = ModPoly []

-- | Tests if a module element is zero.
isZeroModPoly :: ModPoly c -> Bool
isZeroModPoly (ModPoly terms) = null terms

-- | Converts a list of component polynomials [p_0, p_1, ..., p_{m-1}] into a module element
-- sum_{i=0}^{m-1} p_i * e_i.
modPolyFromComponents :: [Poly c] -> ModPoly c
modPolyFromComponents polys =
  let terms = [ ModTerm c m compIdx
              | (compIdx, Poly pTerms _) <- zip [0 :: Int ..] polys
              , Term c m <- pTerms
              ]
  in ModPoly terms

-- | Converts a module element into its component polynomials [p_0, ..., p_{m-1}]
-- of length rank (where rank >= maximum component index + 1).
modPolyToComponents :: (Eq c, Num c) => Int -> MonomialOrder -> ModPoly c -> [Poly c]
modPolyToComponents rank order (ModPoly terms) =
  let grouped = foldl (\acc (ModTerm c m compIdx) ->
                        Map.insertWith (++) compIdx [Term c m] acc
                      ) Map.empty terms
  in [ case Map.lookup compIdx grouped of
         Nothing     -> zeroPoly order
         Just cTerms -> Poly cTerms order
     | compIdx <- [0 .. rank - 1]
     ]

-- | Addition of two module elements.
addModPoly :: (Eq c, Num c) => ModPoly c -> ModPoly c -> ModPoly c
addModPoly (ModPoly terms1) (ModPoly terms2) =
  let combined = foldl (\acc (ModTerm c m compIdx) ->
                         Map.insertWith (+) (m, compIdx) c acc
                       ) Map.empty (terms1 ++ terms2)
      filtered = [ ModTerm c m compIdx
                 | ((m, compIdx), c) <- Map.toList combined
                 , c /= 0
                 ]
  in ModPoly filtered

-- | Subtraction of two module elements.
subModPoly :: (Eq c, Num c) => ModPoly c -> ModPoly c -> ModPoly c
subModPoly p1 p2 = addModPoly p1 (scaleModPoly (-1) p2)

-- | Scalar multiplication of a module element by a constant in k.
scaleModPoly :: (Eq c, Num c) => c -> ModPoly c -> ModPoly c
scaleModPoly scalar (ModPoly terms)
  | scalar == 0 = zeroModPoly
  | otherwise   = ModPoly [ ModTerm (scalar * c) m compIdx | ModTerm c m compIdx <- terms ]

-- | Module multiplication: multiplies a polynomial f in R with a module element v in R^m,
-- producing f * v in R^m.
mulPolyModPoly :: (Eq c, Num c) => Poly c -> ModPoly c -> ModPoly c
mulPolyModPoly (Poly polyTerms _) (ModPoly modTerms) =
  let products = [ ModTerm (cF * cM) (sparseMul mF mM) compIdx
                 | Term cF mF <- polyTerms
                 , ModTerm cM mM compIdx <- modTerms
                 ]
      combined = foldl (\acc (ModTerm c m compIdx) ->
                         Map.insertWith (+) (m, compIdx) c acc
                       ) Map.empty products
      filtered = [ ModTerm c m compIdx
                 | ((m, compIdx), c) <- Map.toList combined
                 , c /= 0
                 ]
  in ModPoly filtered

-- | Computes the degree of a module element given generator degree shifts [d_0, ..., d_{m-1}].
-- The degree of a module term m * e_i is deg(m) + d_i.
-- Returns Nothing for the zero module element.
modPolyDegree :: [Int] -> ModPoly c -> Maybe Int
modPolyDegree shifts (ModPoly terms)
  | null terms = Nothing
  | otherwise =
      let termDegrees = [ sparseDegree m + getShift compIdx
                        | ModTerm _ m compIdx <- terms
                        ]
      in Just (maximum termDegrees)
  where
    getShift compIdx
      | compIdx < length shifts = shifts !! compIdx
      | otherwise = 0

-- | Leading term of a module element (first term in canonical sorted order).
modLeadingTerm :: ModPoly c -> Maybe (ModTerm c)
modLeadingTerm (ModPoly terms) =
  case terms of
    (t:_) -> Just t
    []    -> Nothing

-- | Leading monomial of a module element.
modLeadingMonomial :: ModPoly c -> Maybe ModMonomial
modLeadingMonomial p = do
  ModTerm _ m compIdx <- modLeadingTerm p
  return (ModMonomial m compIdx)

-- | Leading coefficient of a module element.
modLeadingCoeff :: ModPoly c -> Maybe c
modLeadingCoeff p = do
  ModTerm c _ _ <- modLeadingTerm p
  return c
