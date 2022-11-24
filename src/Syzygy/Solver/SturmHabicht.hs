{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Solver.SturmHabicht
  ( sturmHabichtSequence
  , generalizedSturmHabicht
  , tarskiQuery
  , principalSturmHabichtCoeffs
  ) where

import Syzygy.Polynomial.Univariate
  ( UnivariatePoly(..)
  , uDegree
  , uLead
  , uDerivative
  , uSubresultantPRS
  )
import Syzygy.Solver.RealRoot
  ( countSignVariations
  , uEval
  , cauchyBound
  )

-- | Scales a univariate polynomial by a scalar.
uScale :: Num c => c -> UnivariatePoly c -> UnivariatePoly c
uScale c (UnivariatePoly cs) = UnivariatePoly (map (c *) cs)

-- | Computes the Sturm-Habicht sequence SH(P) of a univariate polynomial P
-- using signed subresultant polynomials.
-- Avoids coefficient explosion while preserving exact real root counts under specialization.
sturmHabichtSequence :: forall a. (Eq a, Fractional a)
                     => UnivariatePoly a
                     -> [UnivariatePoly a]
sturmHabichtSequence p =
  let pPrime = uDerivative p
  in generalizedSturmHabicht p pPrime

-- | Computes the generalized Sturm-Habicht sequence SH(P, Q):
-- SH_k(P, Q) = (-1)^{(p-k)(p-k-1)/2} * sRes_k(P, Q)
-- where sRes_k is the k-th polynomial subresultant.
generalizedSturmHabicht :: forall a. (Eq a, Fractional a)
                        => UnivariatePoly a
                        -> UnivariatePoly a
                        -> [UnivariatePoly a]
generalizedSturmHabicht p q =
  let prs = uSubresultantPRS p q
      degP = uDegree p
  in [ let k = uDegree sPoly
           sign = if even (((degP - k) * (degP - k - 1)) `div` 2) then 1 else -1
       in uScale sign sPoly
     | sPoly <- prs
     ]

-- | Computes the Tarski Query TaQ(Q, P) = #{x | P(x) = 0, Q(x) > 0} - #{x | P(x) = 0, Q(x) < 0}:
-- Evaluates the sign variations of the generalized Sturm-Habicht sequence SH(P, Q)
-- at -infty and +infty:
-- TaQ(Q, P) = V(SH(P, Q; -infty)) - V(SH(P, Q; +infty)).
tarskiQuery :: forall a. (Ord a, Fractional a)
            => UnivariatePoly a
            -> UnivariatePoly a
            -> Int
tarskiQuery p q =
  let shSeq = generalizedSturmHabicht p q
      boundP = cauchyBound p
      boundQ = cauchyBound q
      bigM = max boundP boundQ + 1
      valsNegInf = [ uEval s (-bigM) | s <- shSeq ]
      valsPosInf = [ uEval s bigM | s <- shSeq ]
      vNeg = countSignVariations valsNegInf
      vPos = countSignVariations valsPosInf
  in vNeg - vPos

-- | Extracts the principal subresultant coefficients (leading coefficients of the Sturm-Habicht sequence).
principalSturmHabichtCoeffs :: (Eq a, Fractional a) => UnivariatePoly a -> [a]
principalSturmHabichtCoeffs p =
  [ uLead s | s <- sturmHabichtSequence p ]
