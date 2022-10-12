{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.F5
  ( Signature(..)
  , SigPoly(..)
  , makeSigPoly
  , compareSignature
  , isRegularReduction
  , f5SyzygyCriterion
  , f5RewrittenCriterion
  , f5SignatureBuchberger
  , isRegularSequence
  ) where

import GHC.Generics (Generic)
import Syzygy.Monomial.Sparse
  ( SparseMonomial(..)
  , sparseOne
  , sparseMul
  , sparseDivides
  , sparseDiv
  , sparseMaxVar
  , compareSparse
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , isZeroPoly
  , leadingTerm
  , leadingMonomial
  , polyOrder
  )
import Syzygy.Polynomial.Arithmetic
  ( subPoly
  , mulTermPoly
  , monicPoly
  )
import Syzygy.Polynomial.Division (normalForm)
import Syzygy.Grobner.SPolynomial (sPolynomialWithMultipliers)

-- | Module signature (u, e_k) representing a monomial u times the k-th basis generator e_k.
-- Signatures are ordered using Position-Over-Term (POT) or Term-Over-Position (TOP).
data Signature = Signature
  { sigIndex    :: !Int             -- ^ Index k of the canonical basis generator e_k
  , sigMonomial :: !SparseMonomial  -- ^ Monomial multiplier u
  } deriving (Eq, Generic)

-- | A polynomial paired with its exact module signature.
data SigPoly c = SigPoly
  { spSignature  :: !Signature
  , spPolynomial :: !(Poly c)
  } deriving (Eq, Generic)

-- | Constructs a SigPoly for the k-th input generator with signature (1, e_k).
makeSigPoly :: Int -> Poly c -> SigPoly c
makeSigPoly k p = SigPoly (Signature k sparseOne) p

-- | Compares two signatures under Position-Over-Term (POT) order:
-- First compares the generator index k (higher index is greater).
-- Breaks ties by comparing monomials under the given MonomialOrder.
compareSignature :: Poly c -> Signature -> Signature -> Ordering
compareSignature pRef (Signature k1 m1) (Signature k2 m2) =
  case compare k1 k2 of
    EQ ->
      let ord = polyOrder pRef
          nVars = max (sparseMaxVar m1 + 1) (sparseMaxVar m2 + 1)
      in compareSparse ord nVars m1 m2
    other -> other

-- | Regular reduction criterion:
-- Polynomial f with signature sig(f) can be regularly reduced by g with signature sig(g)
-- via monomial term t = c * m only if sig(t * g) < sig(f).
-- This strictly preserves the signature of f throughout all reductions!
isRegularReduction :: Poly c -> Signature -> Signature -> SparseMonomial -> Bool
isRegularReduction pRef sigF (Signature kG mG) mQuot =
  let sigShiftedG = Signature kG (sparseMul mQuot mG)
  in compareSignature pRef sigShiftedG sigF == LT

-- | F5 Syzygy Criterion:
-- Checks whether signature (u, e_k) is divisible by the leading term of any known principal syzygy:
-- LM(f_j) * e_k - LM(f_k) * e_j = 0.
-- If u is divisible by LM(f_j) for some j < k, the signature will reduce to zero and must be rejected!
f5SyzygyCriterion :: [Poly c] -> Signature -> Bool
f5SyzygyCriterion basis (Signature k u) =
  let priorGenerators = take k basis
  in any (\pj -> case leadingMonomial pj of
                   Just lmJ -> sparseDivides lmJ u
                   Nothing  -> False
         ) priorGenerators

-- | F5 Rewritten Criterion (Gao-Volny-Wang):
-- Rejects a pair if its signature has already been produced or rewritten by a smaller generator.
f5RewrittenCriterion :: [SigPoly c] -> Signature -> Bool
f5RewrittenCriterion existingSigPolys targetSig =
  any (\sp ->
         let Signature kE mE = spSignature sp
         in kE == sigIndex targetSig &&
            sparseDivides mE (sigMonomial targetSig) &&
            mE /= sigMonomial targetSig
      ) existingSigPolys

-- | Regular reduction of a SigPoly against a list of SigPolys:
-- reduces f using only regular reductions (sig(t * g) < sig(f)).
regularNormalForm :: forall c. (Eq c, Fractional c)
                  => SigPoly c -> [SigPoly c] -> SigPoly c
regularNormalForm sigF [] = sigF
regularNormalForm sigF divisors = go sigF
  where
    pRef = spPolynomial sigF

    go (SigPoly sig currP)
      | isZeroPoly currP = SigPoly sig currP
      | otherwise =
          case findRegularReductor sig currP divisors of
            Just (termQuot, redPoly) ->
              let subtrahend = mulTermPoly termQuot (spPolynomial redPoly)
                  nextP = subPoly currP subtrahend
              in go (SigPoly sig nextP)
            Nothing ->
              SigPoly sig currP

    findRegularReductor sig currP divs =
      case leadingTerm currP of
        Just (Term cP mP) ->
          let check [] = Nothing
              check (d:ds) =
                case leadingTerm (spPolynomial d) of
                  Just (Term cD mD)
                    | sparseDivides mD mP ->
                        let mQuot = sparseDiv mD mP
                            cQuot = cP / cD
                            tQuot = Term cQuot mQuot
                        in if isRegularReduction pRef sig (spSignature d) mQuot
                           then Just (tQuot, d)
                           else check ds
                    | otherwise -> check ds
                  Nothing -> check ds
          in check divs
        Nothing -> Nothing

-- | Signature-based Groebner basis algorithm (F5/GVW framework).
-- Incremental computation over generators f_0, f_1, ..., f_{m-1}.
-- Rejects all critical pairs violating the F5 Syzygy Criterion or Rewritten Criterion,
-- completely eliminating redundant zero reductions on regular sequences!
f5SignatureBuchberger :: forall c. (Eq c, Fractional c) => [Poly c] -> [Poly c]
f5SignatureBuchberger [] = []
f5SignatureBuchberger inputPolys =
  let nonZeroInputs = [ monicPoly p | p <- inputPolys, not (isZeroPoly p) ]
      sigInputs = zipWith makeSigPoly [0 :: Int ..] nonZeroInputs
  in map spPolynomial (foldl addGenerator [] sigInputs)
  where
    addGenerator basis g =
      let gReduced = regularNormalForm g basis
      in if isZeroPoly (spPolynomial gReduced)
         then basis
         else
           let candidatePairs = [ (sp, gReduced) | sp <- basis ]
               validPairs = filter (isValidF5Pair (map spPolynomial (basis ++ [gReduced]))) candidatePairs
           in loop (basis ++ [gReduced]) validPairs

    isValidF5Pair currentBasis (sp1, sp2) =
      let p1 = spPolynomial sp1
          p2 = spPolynomial sp2
          (_, t1, _) = sPolynomialWithMultipliers p1 p2
          Term _ u1 = t1
          sigS = Signature (sigIndex (spSignature sp1)) (sparseMul u1 (sigMonomial (spSignature sp1)))
      in not (f5SyzygyCriterion currentBasis sigS)

    loop basis [] = basis
    loop basis ((sp1, sp2):pairs) =
      let p1 = spPolynomial sp1
          p2 = spPolynomial sp2
          (sPoly, t1, _) = sPolynomialWithMultipliers p1 p2
          Term _ u1 = t1
          sigSPoly = Signature (sigIndex (spSignature sp1)) (sparseMul u1 (sigMonomial (spSignature sp1)))
          sigS = SigPoly sigSPoly sPoly
          reduced = regularNormalForm sigS basis
      in if isZeroPoly (spPolynomial reduced)
         then loop basis pairs
         else
           let monicP = monicPoly (spPolynomial reduced)
               newSigPoly = SigPoly (spSignature reduced) monicP
               newPairs = [ (sp, newSigPoly) | sp <- basis ]
               validNewPairs = filter (isValidF5Pair (map spPolynomial (basis ++ [newSigPoly]))) newPairs
           in loop (basis ++ [newSigPoly]) (pairs ++ validNewPairs)

-- | Tests whether an ideal basis forms a regular sequence:
-- (f_0, ..., f_{m-1}) is a regular sequence if each f_k is a non-zero divisor in R / (f_0, ..., f_{k-1}).
-- In F5, a regular sequence produces ZERO redundant zero reductions!
isRegularSequence :: (Eq c, Fractional c) => [Poly c] -> Bool
isRegularSequence [] = True
isRegularSequence [_] = True
isRegularSequence (f0:fs) =
  let go _ [] = True
      go currentIdeal (g:rest) =
        let remG = normalForm g currentIdeal
        in if isZeroPoly remG
           then False
           else go (currentIdeal ++ [g]) rest
  in go [f0] fs

instance Show Signature where
  show (Signature k m) = "(" ++ show m ++ " * e_" ++ show k ++ ")"

instance (Show c, Eq c, Ord c, Num c) => Show (SigPoly c) where
  show (SigPoly sig p) = show sig ++ " | " ++ show p
