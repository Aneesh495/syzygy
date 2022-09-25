{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Pretty
  ( prettyPoly
  , prettyUnicode
  , prettyLaTeX
  , prettyANSI
  ) where

import Syzygy.Core.Rational (Rat(..))
import Syzygy.Monomial.Sparse (sparseOne, sparseSupport, sparseGetExp)
import Syzygy.Polynomial.Base (Term(..), Poly(..), isZeroPoly, polyTerms)

-- | Converts an exponent integer to Unicode superscript characters: 0 -> ⁰, 2 -> ², etc.
toUnicodeSuperscript :: Int -> String
toUnicodeSuperscript n
  | n < 0     = '⁻' : toUnicodeSuperscript (-n)
  | otherwise = map digitToSuper (show n)
  where
    digitToSuper '0' = '⁰'
    digitToSuper '1' = '¹'
    digitToSuper '2' = '²'
    digitToSuper '3' = '³'
    digitToSuper '4' = '⁴'
    digitToSuper '5' = '⁵'
    digitToSuper '6' = '⁶'
    digitToSuper '7' = '⁷'
    digitToSuper '8' = '⁸'
    digitToSuper '9' = '⁹'
    digitToSuper c   = c

-- | Resolves a 0-indexed variable ID to a variable name string.
varName :: [String] -> Int -> String
varName names idx
  | idx < length names = names !! idx
  | otherwise          = "x_" ++ show idx

-- | Pretty-prints a polynomial in standard ASCII using named variables.
prettyPoly :: [String] -> Poly Rat -> String
prettyPoly names p
  | isZeroPoly p = "0"
  | otherwise =
      let terms = polyTerms p
      in case terms of
           [] -> "0"
           (t0:ts) -> formatTerm True t0 ++ concatMap (formatTerm False) ts
  where
    formatTerm isFirst (Term (Rat num den) m) =
      let signStr = if isFirst
                    then (if num < 0 then "-" else "")
                    else (if num < 0 then " - " else " + ")
          absNum = abs num
          coeffStr
            | den == 1 && absNum == 1 && m /= sparseOne = ""
            | den == 1                                  = show absNum
            | otherwise                                 = show absNum ++ "/" ++ show den
          monoStr = formatMonomial m
          sep = if not (null coeffStr) && not (null monoStr) then "*" else ""
      in signStr ++ coeffStr ++ sep ++ monoStr

    formatMonomial m
      | m == sparseOne = ""
      | otherwise =
          let vars = [ varName names i ++ (if e == 1 then "" else "^" ++ show e)
                     | i <- sparseSupport m
                     , let e = sparseGetExp i m
                     , e > 0
                     ]
          in unwords vars

-- | Pretty-prints a polynomial using Unicode mathematical superscripts (e.g. x² + 2xy - y³).
prettyUnicode :: [String] -> Poly Rat -> String
prettyUnicode names p
  | isZeroPoly p = "0"
  | otherwise =
      let terms = polyTerms p
      in case terms of
           [] -> "0"
           (t0:ts) -> formatTerm True t0 ++ concatMap (formatTerm False) ts
  where
    formatTerm isFirst (Term (Rat num den) m) =
      let signStr = if isFirst
                    then (if num < 0 then "-" else "")
                    else (if num < 0 then " - " else " + ")
          absNum = abs num
          coeffStr
            | den == 1 && absNum == 1 && m /= sparseOne = ""
            | den == 1                                  = show absNum
            | otherwise                                 = show absNum ++ "/" ++ show den
          monoStr = formatUnicodeMonomial m
      in signStr ++ coeffStr ++ monoStr

    formatUnicodeMonomial m
      | m == sparseOne = ""
      | otherwise =
          let vars = [ varName names i ++ (if e == 1 then "" else toUnicodeSuperscript e)
                     | i <- sparseSupport m
                     , let e = sparseGetExp i m
                     , e > 0
                     ]
          in concat vars

-- | Pretty-prints a polynomial in LaTeX math syntax (e.g. x^2 + 2xy - y^3).
prettyLaTeX :: [String] -> Poly Rat -> String
prettyLaTeX names p
  | isZeroPoly p = "0"
  | otherwise =
      let terms = polyTerms p
      in case terms of
           [] -> "0"
           (t0:ts) -> formatTerm True t0 ++ concatMap (formatTerm False) ts
  where
    formatTerm isFirst (Term (Rat num den) m) =
      let signStr = if isFirst
                    then (if num < 0 then "-" else "")
                    else (if num < 0 then " - " else " + ")
          absNum = abs num
          coeffStr
            | den == 1 && absNum == 1 && m /= sparseOne = ""
            | den == 1                                  = show absNum
            | otherwise                                 = "\\frac{" ++ show absNum ++ "}{" ++ show den ++ "}"
          monoStr = formatLaTeXMonomial m
      in signStr ++ coeffStr ++ monoStr

    formatLaTeXMonomial m
      | m == sparseOne = ""
      | otherwise =
          let vars = [ varName names i ++ (if e == 1 then "" else "^{" ++ show e ++ "}")
                     | i <- sparseSupport m
                     , let e = sparseGetExp i m
                     , e > 0
                     ]
          in concat vars

-- | Pretty-prints a polynomial with ANSI terminal colors for high-contrast mathematical REPLs.
prettyANSI :: [String] -> Poly Rat -> String
prettyANSI names p
  | isZeroPoly p = "\ESC[33m0\ESC[0m"
  | otherwise =
      let terms = polyTerms p
      in case terms of
           [] -> "\ESC[33m0\ESC[0m"
           (t0:ts) -> formatTerm True t0 ++ concatMap (formatTerm False) ts
  where
    cyan = "\ESC[36m"
    green = "\ESC[32m"
    yellow = "\ESC[33m"
    magenta = "\ESC[35m"
    reset = "\ESC[0m"

    formatTerm isFirst (Term (Rat num den) m) =
      let signStr = if isFirst
                    then (if num < 0 then magenta ++ "-" ++ reset else "")
                    else (if num < 0 then magenta ++ " - " ++ reset else magenta ++ " + " ++ reset)
          absNum = abs num
          coeffStr
            | den == 1 && absNum == 1 && m /= sparseOne = ""
            | den == 1                                  = green ++ show absNum ++ reset
            | otherwise                                 = green ++ show absNum ++ "/" ++ show den ++ reset
          monoStr = formatANSIMonomial m
      in signStr ++ coeffStr ++ monoStr

    formatANSIMonomial m
      | m == sparseOne = ""
      | otherwise =
          let vars = [ cyan ++ varName names i ++ reset ++ (if e == 1 then "" else yellow ++ "^" ++ show e ++ reset)
                     | i <- sparseSupport m
                     , let e = sparseGetExp i m
                     , e > 0
                     ]
          in unwords vars
