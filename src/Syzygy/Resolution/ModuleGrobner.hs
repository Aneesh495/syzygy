{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Resolution.ModuleGrobner
  ( moduleDivision
  , moduleNormalForm
  , moduleSPolynomial
  , moduleBuchberger
  , moduleReduce
  , dividesModLeadingTerm
  ) where

import qualified Data.List as List
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse
  ( sparseDivides
  , sparseDiv
  , sparseLCM
  )
import Syzygy.Polynomial.Base
  ( Poly(..)
  , Term(..)
  , zeroPoly
  )
import Syzygy.Polynomial.Arithmetic (addPoly)
import Syzygy.Resolution.Module
  ( ModMonomial(..)
  , ModTerm(..)
  , ModPoly(..)
  , zeroModPoly
  , isZeroModPoly
  , subModPoly
  , scaleModPoly
  , mulPolyModPoly
  , modLeadingTerm
  , modLeadingMonomial
  , modLeadingCoeff
  )
import Syzygy.Resolution.ModuleOrder
  ( ModuleOrder(..)
  , canonicalizeModPoly
  )

-- | Checks if the leading term of g divides the leading term of f in the free module.
-- Division is only possible if both terms belong to the SAME module component
-- and the monomial of g divides the monomial of f.
-- If true, returns Just (scalarCoeff, monomialFactor) such that
-- (scalarCoeff * monomialFactor) * LT(g) = LT(f).
dividesModLeadingTerm :: (Eq c, Fractional c)
                      => ModPoly c
                      -> ModPoly c
                      -> Maybe (c, Poly c)
dividesModLeadingTerm f g =
  case (modLeadingTerm f, modLeadingTerm g) of
    (Just (ModTerm cF mF compF), Just (ModTerm cG mG compG))
      | compF == compG && sparseDivides mG mF ->
          let qCoeff = cF / cG
              qExp = sparseDiv mF mG
              polyTerm = Poly [Term qCoeff qExp] DegRevLex
          in Just (qCoeff, polyTerm)
      | otherwise -> Nothing
    _ -> Nothing

-- | Multivariate division algorithm for free module elements over R = k[x_0, ..., x_{n-1}].
-- Given f in R^m and divisors [g_1, ..., g_s] in R^m, computes quotients [q_1, ..., q_s] in R
-- and remainder r in R^m such that:
--   f = sum_{k=1}^s q_k * g_k + r
-- where no term of r is divisible by any leading term of g_k.
moduleDivision :: forall c. (Eq c, Fractional c)
               => ModuleOrder
               -> Int
               -> ModPoly c
               -> [ModPoly c]
               -> ([Poly c], ModPoly c)
moduleDivision order numVars f divisors =
  let nonZeroDivs = filter (not . isZeroModPoly) divisors
      numDivs = length divisors
      initQuotients = replicate numDivs (zeroPoly (extractPolyOrder order))
  in divLoop f initQuotients zeroModPoly nonZeroDivs
  where
    extractPolyOrder (TOP po _) = po
    extractPolyOrder (POT po) = po
    extractPolyOrder (SchreyerOrder bo _) = extractPolyOrder bo

    divLoop !p !quotients !remTerms divs
      | isZeroModPoly p =
          let finalRem = canonicalizeModPoly order numVars remTerms
          in (quotients, finalRem)
      | otherwise =
          case findReducer p divs of
            Just (divIdx, qPoly) ->
              let qDiv = mulPolyModPoly qPoly (divisors !! divIdx)
                  p' = canonicalizeModPoly order numVars (subModPoly p qDiv)
                  quotients' = updateAt divIdx (\q -> addPoly q qPoly) quotients
              in divLoop p' quotients' remTerms divs
            Nothing ->
              case modLeadingTerm p of
                Nothing -> (quotients, canonicalizeModPoly order numVars remTerms)
                Just lt ->
                  let pRest = ModPoly (drop 1 (modPolyTerms p))
                      remTerms' = ModPoly (lt : modPolyTerms remTerms)
                  in divLoop pRest quotients remTerms' divs

    findReducer p _divs =
      let candidates = [ (i, ltG)
                       | (i, g) <- zip [0 :: Int ..] divisors
                       , not (isZeroModPoly g)
                       , Just ltG <- [modLeadingTerm g]
                       ]
      in case modLeadingTerm p of
           Nothing -> Nothing
           Just (ModTerm cP mP compP) ->
             case filter (\(_, ModTerm _ mG compG) -> compG == compP && sparseDivides mG mP) candidates of
               ((i, ModTerm cG mG _):_) ->
                 let qCoeff = cP / cG
                     qExp = sparseDiv mP mG
                     qPoly = Poly [Term qCoeff qExp] (extractPolyOrder order)
                 in Just (i, qPoly)
               [] -> Nothing

    updateAt idx fn xs =
      [ if i == idx then fn x else x | (i, x) <- zip [0 :: Int ..] xs ]

-- | Computes the normal form (remainder) of a module element with respect to a list of generators.
moduleNormalForm :: (Eq c, Fractional c)
                 => ModuleOrder
                 -> Int
                 -> ModPoly c
                 -> [ModPoly c]
                 -> ModPoly c
moduleNormalForm order numVars f divs = snd (moduleDivision order numVars f divs)

-- | Computes the S-polynomial S(g_a, g_b) between two module elements.
-- If the leading terms of g_a and g_b are in DIFFERENT module components,
-- their S-polynomial is identically zero: S(g_a, g_b) = 0.
-- If they share the same component e_i:
--   S(g_a, g_b) = (LCM(m_a, m_b) / LT(g_a)) * g_a - (LCM(m_a, m_b) / LT(g_b)) * g_b.
moduleSPolynomial :: forall c. (Eq c, Fractional c)
                  => ModuleOrder
                  -> Int
                  -> ModPoly c
                  -> ModPoly c
                  -> ModPoly c
moduleSPolynomial order numVars gA gB =
  case (modLeadingTerm gA, modLeadingTerm gB) of
    (Just (ModTerm cA mA compA), Just (ModTerm cB mB compB))
      | compA == compB ->
          let lcmM = sparseLCM mA mB
              pOrder = extractPolyOrder order
              multA = Poly [Term (1 / cA) (sparseDiv lcmM mA)] pOrder
              multB = Poly [Term (1 / cB) (sparseDiv lcmM mB)] pOrder
              termA = mulPolyModPoly multA gA
              termB = mulPolyModPoly multB gB
          in canonicalizeModPoly order numVars (subModPoly termA termB)
      | otherwise -> zeroModPoly
    _ -> zeroModPoly
  where
    extractPolyOrder (TOP po _) = po
    extractPolyOrder (POT po) = po
    extractPolyOrder (SchreyerOrder bo _) = extractPolyOrder bo

-- | Submodule Buchberger algorithm for computing a Groebner basis of a submodule M subset R^m.
moduleBuchberger :: forall c. (Eq c, Fractional c)
                 => ModuleOrder
                 -> Int
                 -> [ModPoly c]
                 -> [ModPoly c]
moduleBuchberger order numVars initialGens =
  let canonGens = filter (not . isZeroModPoly) (map (canonicalizeModPoly order numVars) initialGens)
      initialPairs = [ (i, j)
                     | i <- [0 .. length canonGens - 1]
                     , j <- [i + 1 .. length canonGens - 1]
                     , shareComponent (canonGens !! i) (canonGens !! j)
                     ]
  in go canonGens initialPairs
  where
    shareComponent g1 g2 =
      case (modLeadingMonomial g1, modLeadingMonomial g2) of
        (Just (ModMonomial _ c1), Just (ModMonomial _ c2)) -> c1 == c2
        _ -> False

    go !basis [] = basis
    go !basis ((i, j):remPairs) =
      let gI = basis !! i
          gJ = basis !! j
          sPoly = moduleSPolynomial order numVars gI gJ
          remPoly = moduleNormalForm order numVars sPoly basis
      in if isZeroModPoly remPoly
         then go basis remPairs
         else
           let newIdx = length basis
               newBasis = basis ++ [remPoly]
               newPairs = [ (k, newIdx)
                          | k <- [0 .. newIdx - 1]
                          , shareComponent (newBasis !! k) remPoly
                          ]
           in go newBasis (remPairs ++ newPairs)

-- | Inter-reduces a submodule Groebner basis into a minimal and reduced Groebner basis.
moduleReduce :: forall c. (Eq c, Fractional c)
             => ModuleOrder
             -> Int
             -> [ModPoly c]
             -> [ModPoly c]
moduleReduce order numVars gb =
  let nonZeroGB = filter (not . isZeroModPoly) gb
      monicGB = map makeMonic nonZeroGB
      minimalGB = pruneRedundant monicGB
  in fullyReduce minimalGB
  where
    makeMonic p =
      case modLeadingCoeff p of
        Nothing -> p
        Just c  -> scaleModPoly (recip c) p

    pruneRedundant [] = []
    pruneRedundant (p:ps) =
      if any (`divides` p) ps || any (`divides` p) (pruneRedundant ps)
      then pruneRedundant ps
      else p : pruneRedundant [ q | q <- ps, not (divides p q) ]

    divides g f =
      case (modLeadingMonomial g, modLeadingMonomial f) of
        (Just (ModMonomial mG cG), Just (ModMonomial mF cF)) ->
          cG == cF && sparseDivides mG mF
        _ -> False

    fullyReduce [] = []
    fullyReduce basis =
      [ reduceOne b (List.delete b basis) | b <- basis ]

    reduceOne p others =
      let (_, remPoly) = moduleDivision order numVars p others
      in makeMonic remPoly
