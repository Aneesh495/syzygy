{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.UI.BettiTable
Description : Pretty-printing and multi-format export for graded Betti tables
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Provides renderers for graded Betti tables in multiple formats:
- ANSI colored terminal output for interactive REPL use
- Classic Macaulay2 and Singular text format
- LaTeX matrix environment for academic publications
- GitHub-flavored Markdown tables
- JSON encoding for WebGL variety visualizer integration
-}
module Syzygy.UI.BettiTable
  ( renderBettiTableClassic
  , renderBettiTableAnsi
  , renderBettiTableLatex
  , renderBettiTableMarkdown
  , renderBettiTableJson
  , renderBettiDiagnostics
  ) where

import Data.List (intercalate)
import Syzygy.Resolution.Betti
  ( BettiTable(..)
  )

-- | Classic Macaulay2 style text output with dots for zeroes.
renderBettiTableClassic :: BettiTable -> String
renderBettiTableClassic bt =
  let totals = bettiTotal bt
      pDim = bettiProjDim bt
      colWidth = 5
      padLeft n s = replicate (max 0 (n - length s)) ' ' ++ s

      colIndices = concatMap (padLeft colWidth . show) [0 .. pDim]
      headerTop  = "             " ++ colIndices
      headerTot  = "       total:" ++ concatMap (padLeft colWidth . show) totals
      separator  = "       " ++ replicate (7 + (pDim + 1) * colWidth) '-'

      formatRow (rowIdx, entries) =
        let rowLabel = padLeft 5 (show rowIdx) ++ ": "
            rowVals  = concatMap (\v -> padLeft colWidth (if v == 0 then "." else show v)) entries
        in rowLabel ++ rowVals

      rows = zip [bettiMinRow bt .. bettiMaxRow bt] (bettiMatrix bt)
      rowStrs = map formatRow rows
  in unlines ([headerTop, headerTot, separator] ++ rowStrs)

-- | ANSI 256-color enhanced terminal output.
renderBettiTableAnsi :: BettiTable -> String
renderBettiTableAnsi bt =
  let totals = bettiTotal bt
      pDim = bettiProjDim bt
      colWidth = 6
      padLeft n s = replicate (max 0 (n - length s)) ' ' ++ s

      bold = "\ESC[1m"
      cyan = "\ESC[36m"
      green = "\ESC[32m"
      yellow = "\ESC[33m"
      gray = "\ESC[90m"
      reset = "\ESC[0m"

      colIndices = concatMap (\i -> cyan ++ padLeft colWidth (show i) ++ reset) [0 .. pDim]
      headerTop  = "             " ++ colIndices
      headerTot  = bold ++ "       total:" ++ reset ++ concatMap (\t -> bold ++ green ++ padLeft colWidth (show t) ++ reset) totals
      separator  = gray ++ "       " ++ replicate (7 + (pDim + 1) * colWidth) '-' ++ reset

      formatRow (rowIdx, entries) =
        let rowLabel = cyan ++ padLeft 5 (show rowIdx) ++ ": " ++ reset
            rowVals  = concatMap (\v ->
                         if v == 0
                         then gray ++ padLeft colWidth "." ++ reset
                         else bold ++ yellow ++ padLeft colWidth (show v) ++ reset
                       ) entries
        in rowLabel ++ rowVals

      rows = zip [bettiMinRow bt .. bettiMaxRow bt] (bettiMatrix bt)
      rowStrs = map formatRow rows
  in unlines ([headerTop, headerTot, separator] ++ rowStrs)

-- | Formats the Betti table as a LaTeX matrix suitable for insertion into papers.
renderBettiTableLatex :: BettiTable -> String
renderBettiTableLatex bt =
  let totals = bettiTotal bt
      pDim = bettiProjDim bt
      colsHeader = " & " ++ intercalate " & " (map show [0 .. pDim]) ++ " \\\\"
      totRow     = "\\text{total} & " ++ intercalate " & " (map show totals) ++ " \\\\"
      formatRow (rowIdx, entries) =
        show rowIdx ++ " & " ++ intercalate " & " (map (\v -> if v == 0 then "\\cdot" else show v) entries) ++ " \\\\"
      rows = zip [bettiMinRow bt .. bettiMaxRow bt] (bettiMatrix bt)
      bodyRows = map formatRow rows
      colSpec = "r|" ++ replicate (pDim + 1) 'r'
  in unlines $
       [ "\\begin{array}{" ++ colSpec ++ "}"
       , colsHeader
       , "\\hline"
       , totRow
       , "\\hline"
       ] ++ bodyRows ++
       [ "\\end{array}" ]

-- | Formats the Betti table as a GitHub-flavored Markdown table.
renderBettiTableMarkdown :: BettiTable -> String
renderBettiTableMarkdown bt =
  let totals = bettiTotal bt
      pDim = bettiProjDim bt
      header = "| Degree | " ++ intercalate " | " (map (\i -> "Tor_" ++ show i) [0 .. pDim]) ++ " |"
      divider = "|:---|" ++ concat (replicate (pDim + 1) ":---:|")
      totRow = "| **Total** | " ++ intercalate " | " (map (\t -> "**" ++ show t ++ "**") totals) ++ " |"
      formatRow (rowIdx, entries) =
        "| " ++ show rowIdx ++ " | " ++ intercalate " | " (map (\v -> if v == 0 then "-" else show v) entries) ++ " |"
      rows = zip [bettiMinRow bt .. bettiMaxRow bt] (bettiMatrix bt)
      bodyRows = map formatRow rows
  in unlines ([header, divider, totRow] ++ bodyRows)

-- | Serializes the Betti table into JSON format for front-end visualizers.
renderBettiTableJson :: BettiTable -> String
renderBettiTableJson bt =
  let totals = bettiTotal bt
      matrix = bettiMatrix bt
      showList' xs = "[" ++ intercalate ", " (map show xs) ++ "]"
      showMatrix m = "[" ++ intercalate ", " (map showList' m) ++ "]"
  in "{\n" ++
     "  \"projectiveDimension\": " ++ show (bettiProjDim bt) ++ ",\n" ++
     "  \"regularity\": " ++ show (bettiRegularity bt) ++ ",\n" ++
     "  \"minRow\": " ++ show (bettiMinRow bt) ++ ",\n" ++
     "  \"maxRow\": " ++ show (bettiMaxRow bt) ++ ",\n" ++
     "  \"totals\": " ++ showList' totals ++ ",\n" ++
     "  \"matrix\": " ++ showMatrix matrix ++ "\n" ++
     "}"

-- | Produces a comprehensive homological diagnostic summary for the module.
renderBettiDiagnostics :: BettiTable -> String
renderBettiDiagnostics bt =
  unlines
    [ "=== Homological Invariants ==="
    , "  Projective Dimension : " ++ show (bettiProjDim bt)
    , "  Castelnuovo-Mumford  : " ++ show (bettiRegularity bt)
    , "  Total Syzygy Ranks   : " ++ intercalate " -> " (map show (reverse (bettiTotal bt)))
    , "  Degree Range         : [" ++ show (bettiMinRow bt) ++ " .. " ++ show (bettiMaxRow bt) ++ "]"
    ]
