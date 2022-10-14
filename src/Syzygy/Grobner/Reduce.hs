{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Grobner.Reduce
  ( minimalGrobnerBasis
  , reducedGrobnerBasis
  , canonicalGroebnerBasis
  , idealEquality
  , isReducedGrobnerBasis
  ) where

import Syzygy.Monomial.Sparse
  ( sparseDivides
  )
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , isZeroPoly
  , polyTerms
  , leadingTerm
  , leadingMonomial
  )
import Syzygy.Polynomial.Arithmetic (monicPoly)
import Syzygy.Polynomial.Division (fullReduction)
import Syzygy.Grobner.F4 (f4)

-- | Computes a minimal Groebner basis from an arbitrary Groebner basis:
-- 1. Discards any generator g_i whose leading monomial LM(g_i) is divisible
--    by the leading monomial LM(g_j) of another generator j /= i.
-- 2. Normalizes each surviving polynomial to be monic (LC(g) = 1).
minimalGrobnerBasis :: forall c. (Eq c, Fractional c) => [Poly c] -> [Poly c]
minimalGrobnerBasis rawBasis =
  let nonZero = [ monicPoly p | p <- rawBasis, not (isZeroPoly p) ]
  in filterPruned nonZero []
  where
    filterPruned [] acc = reverse acc
    filterPruned (p:ps) acc =
      let otherGenerators = acc ++ ps
          isRedundant = case leadingMonomial p of
            Just lmP ->
              any (\other ->
                     case leadingMonomial other of
                       Just lmO -> sparseDivides lmO lmP
                       Nothing  -> False
                  ) otherGenerators
            Nothing -> True
      in if isRedundant
         then filterPruned ps acc
         else filterPruned ps (p : acc)

-- | Computes the UNIQUE Reduced Canonical Groebner Basis (Buchberger Uniqueness Theorem):
-- 1. First extracts a minimal Groebner basis.
-- 2. Fully tail-reduces each generator g_i against all other generators G \ {g_i},
--    ensuring no monomial in g_i lies in the leading ideal in(G \ {g_i}).
-- 3. Renormalizes all generators to be monic.
-- For a fixed MonomialOrder, the reduced Groebner basis is canonical and unique!
reducedGrobnerBasis :: forall c. (Eq c, Fractional c) => [Poly c] -> [Poly c]
reducedGrobnerBasis rawBasis =
  let minimal = minimalGrobnerBasis rawBasis
  in go minimal []
  where
    go [] processed = reverse processed
    go (g:gs) processed =
      let otherDivisors = processed ++ gs
          reducedG = fullReduction g otherDivisors
          monicG = monicPoly reducedG
      in if isZeroPoly monicG
         then go gs processed
         else go gs (monicG : processed)

-- | End-to-end solver pipeline: computes the canonical reduced Groebner basis
-- of an ideal from any arbitrary set of generating polynomials using Faugere F4
-- followed by complete inter-reduction.
canonicalGroebnerBasis :: (Eq c, Fractional c) => [Poly c] -> [Poly c]
canonicalGroebnerBasis inputPolys =
  let rawGB = f4 inputPolys
  in reducedGrobnerBasis rawGB

-- | Tests ideal equality: checks whether two polynomial systems generate the exact same ideal:
-- <f_1, ..., f_k> == <g_1, ..., g_m>.
-- Holds if and only if their canonical reduced Groebner bases are identical!
idealEquality :: (Eq c, Fractional c) => [Poly c] -> [Poly c] -> Bool
idealEquality sys1 sys2 =
  let gb1 = canonicalGroebnerBasis sys1
      gb2 = canonicalGroebnerBasis sys2
  in gb1 == gb2

-- | Tests if a given basis satisfies all conditions of a reduced Groebner basis.
isReducedGrobnerBasis :: forall c. (Eq c, Fractional c) => [Poly c] -> Bool
isReducedGrobnerBasis basis
  | null basis = True
  | otherwise =
      let nonZero = filter (not . isZeroPoly) basis
          allMonic = all (\p -> case leadingTerm p of
                                  Just (Term c _) -> c == 1
                                  Nothing         -> False
                         ) nonZero
          allTailReduced = all checkGenerator (zip [0 :: Int ..] nonZero)
          checkGenerator (i, p) =
            let others = [ g | (j, g) <- zip [0 ..] nonZero, j /= i ]
                terms = polyTerms p
            in all (\(Term _ m) ->
                      not (any (\g -> case leadingMonomial g of
                                       Just lmG -> sparseDivides lmG m
                                       Nothing  -> False
                               ) others)
                   ) terms
      in allMonic && allTailReduced
