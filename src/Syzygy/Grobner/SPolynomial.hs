{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.SPolynomial
  ( CriticalPair(..)
  , makeCriticalPair
  , sPolynomial
  , sPolynomialWithMultipliers
  , criticalPairLCM
  , criticalPairDegree
  , criticalPairSugar
  ) where

import GHC.Generics (Generic)
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseDegree
  , sparseLCM
  , sparseDiv
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , zeroPoly
  , isZeroPoly
  , leadingTerm
  , leadingMonomial
  , polyOrder
  )
import Syzygy.Polynomial.Arithmetic
  ( subPoly
  , mulTermPoly
  )

-- | Critical pair between polynomial f_i and f_j with precomputed LCM, degree, and sugar.
data CriticalPair c = CriticalPair
  { cpIndex1 :: !Int             -- ^ Index of first generator in current basis
  , cpIndex2 :: !Int             -- ^ Index of second generator in current basis
  , cpPoly1  :: !(Poly c)        -- ^ Polynomial f_1
  , cpPoly2  :: !(Poly c)        -- ^ Polynomial f_2
  , cpLCM    :: !SparseMonomial  -- ^ lcm(LM(f_1), LM(f_2))
  , cpDegree :: !Int             -- ^ deg(lcm(LM(f_1), LM(f_2)))
  , cpSugar  :: !Int             -- ^ Giovini sugar degree of the pair
  } deriving (Eq, Generic)

-- | Extracts the LCM monomial of the critical pair.
criticalPairLCM :: CriticalPair c -> SparseMonomial
criticalPairLCM = cpLCM

-- | Extracts the degree of the LCM of leading monomials.
criticalPairDegree :: CriticalPair c -> Int
criticalPairDegree = cpDegree

-- | Extracts the sugar degree of the critical pair.
criticalPairSugar :: CriticalPair c -> Int
criticalPairSugar = cpSugar

-- | Smart constructor for a critical pair with sugar calculation.
makeCriticalPair :: (Eq c, Fractional c)
                 => Int -> Int -> Poly c -> Poly c -> Int -> Int -> CriticalPair c
makeCriticalPair idx1 idx2 p1 p2 sugar1 sugar2 =
  case (leadingMonomial p1, leadingMonomial p2) of
    (Just m1, Just m2) ->
      let lcmM = sparseLCM m1 m2
          degLCM = sparseDegree lcmM
          u1 = sparseDiv m1 lcmM
          u2 = sparseDiv m2 lcmM
          degU1 = sparseDegree u1
          degU2 = sparseDegree u2
          sugarPair = max (sugar1 + degU1) (sugar2 + degU2)
      in CriticalPair idx1 idx2 p1 p2 lcmM degLCM sugarPair
    _ -> error "makeCriticalPair: cannot construct critical pair for zero polynomial"

-- | Computes the classical S-polynomial:
-- S(f, g) = (lcm / LT(f)) * f - (lcm / LT(g)) * g.
sPolynomial :: (Eq c, Fractional c) => Poly c -> Poly c -> Poly c
sPolynomial f g
  | isZeroPoly f || isZeroPoly g = zeroPoly (polyOrder f)
  | otherwise =
      case (leadingTerm f, leadingTerm g) of
        (Just (Term cF mF), Just (Term cG mG)) ->
          let lcmM = sparseLCM mF mG
              uF = sparseDiv mF lcmM
              uG = sparseDiv mG lcmM
              tF = Term (recip cF) uF
              tG = Term (recip cG) uG
              term1 = mulTermPoly tF f
              term2 = mulTermPoly tG g
          in subPoly term1 term2
        _ -> zeroPoly (polyOrder f)

-- | Computes the S-polynomial along with the monomial-multiplier terms:
-- returns (S(f, g), Term(1/cF, uF), Term(1/cG, uG)).
sPolynomialWithMultipliers :: (Eq c, Fractional c)
                           => Poly c -> Poly c -> (Poly c, Term c, Term c)
sPolynomialWithMultipliers f g =
  case (leadingTerm f, leadingTerm g) of
    (Just (Term cF mF), Just (Term cG mG)) ->
      let lcmM = sparseLCM mF mG
          uF = sparseDiv mF lcmM
          uG = sparseDiv mG lcmM
          tF = Term (recip cF) uF
          tG = Term (recip cG) uG
          term1 = mulTermPoly tF f
          term2 = mulTermPoly tG g
          sPoly = subPoly term1 term2
      in (sPoly, tF, tG)
    _ -> (zeroPoly (polyOrder f), Term 0 (SparseMonomial mempty), Term 0 (SparseMonomial mempty))

instance Show c => Show (CriticalPair c) where
  show cp =
    "CriticalPair(" ++ show (cpIndex1 cp) ++ ", " ++ show (cpIndex2 cp) ++
    ", deg=" ++ show (cpDegree cp) ++ ", sugar=" ++ show (cpSugar cp) ++ ")"
