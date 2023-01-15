{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

{- |
Module      : Syzygy.UI.CLI
Description : Command line interface, argument parser, and subcommand dispatcher
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Command line interface for Syzygy.
Supports non-interactive batch commands for automated pipelines and an
interactive REPL:
  syzygy repl                          Launch interactive session
  syzygy grobner [OPTIONS] <POLYS...>  Compute Groebner basis
  syzygy betti [OPTIONS] <POLYS...>    Compute minimal resolution and Betti table
  syzygy raymarch <SURFACE>            Render 3D algebraic variety in terminal
  syzygy sheaf <SPACE>                 Compute cellular sheaf cohomology
-}
module Syzygy.UI.CLI
  ( mainCLI
  ) where

import System.Environment (getArgs)
import System.Exit (exitFailure, exitSuccess)
import Data.Char (toLower)
import Data.List (isPrefixOf)
import Syzygy.Core.Rational (Rat)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Polynomial.Base (Poly(..))
import Syzygy.Polynomial.Parser (parsePoly, defaultVarNames)
import Syzygy.Polynomial.Pretty (prettyUnicode)
import Syzygy.Grobner.Buchberger (buchberger)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Resolution.FreeResolution (freeResolution)
import Syzygy.Resolution.Minimization (minimizeComplex)
import Syzygy.Resolution.Betti (computeBettiTable)
import Syzygy.UI.BettiTable (renderBettiTableAnsi, renderBettiDiagnostics)
import Syzygy.UI.VarietyVisualizer
  ( raymarchVariety3D
  , whitneyUmbrella
  , clebschCubic
  , cayleyCubic
  , kummerSurface
  , torusVariety
  , sphereVariety
  )
import Syzygy.Sheaf.Complex
  ( circleComplex
  , sphere2Complex
  , torusComplex
  )
import Syzygy.Sheaf.CellularSheaf (constantSheaf1D)
import Syzygy.Sheaf.Cohomology (allSheafBettiNumbers, sheafEulerCharacteristic)
import Syzygy.UI.REPL (runREPL)

-- | Main entry point for command-line execution.
mainCLI :: IO ()
mainCLI = do
  args <- getArgs
  case args of
    [] -> runREPL
    ("repl" : _) -> runREPL
    ("help" : _) -> printUsage >> exitSuccess
    ("--help" : _) -> printUsage >> exitSuccess
    ("-h" : _) -> printUsage >> exitSuccess
    ("grobner" : rest) -> runGrobner rest
    ("gb" : rest) -> runGrobner rest
    ("betti" : rest) -> runBetti rest
    ("res" : rest) -> runBetti rest
    ("raymarch" : rest) -> runRaymarch rest
    ("sheaf" : rest) -> runSheaf rest
    unknown -> do
      putStrLn $ "Error: Unknown subcommand '" ++ unwords unknown ++ "'"
      printUsage
      exitFailure

-- | Prints CLI usage documentation.
printUsage :: IO ()
printUsage = putStrLn $ unlines
  [ "Syzygy - Computational Algebraic Geometry Engine"
  , "Usage: syzygy <command> [options] [arguments]"
  , ""
  , "Commands:"
  , "  repl                     Start interactive algebraic geometry session (default)"
  , "  grobner [polys...]       Compute reduced Groebner basis"
  , "  betti [polys...]         Compute minimal graded free resolution & Betti table"
  , "  raymarch <surface>       Raymarch 3D variety (whitney|clebsch|cayley|kummer|torus|sphere)"
  , "  sheaf <space>            Compute sheaf cohomology & Betti numbers (circle|sphere|torus)"
  , ""
  , "Options:"
  , "  --order <lex|deglex|grevlex>  Set monomial ordering (default: grevlex)"
  , "  --vars x,y,z                  Set variable names (default: x,y,z,w)"
  ]

-- | Parses command line options and separates polynomial strings.
parseOptions :: [String] -> (MonomialOrder, [String], [String])
parseOptions args = go DegRevLex (take 4 defaultVarNames) [] args
  where
    go ord vars polys [] = (ord, vars, reverse polys)
    go _ vars polys ("--order" : o : rest) =
      let newOrd = case map toLower o of
                     "lex"     -> Lex
                     "deglex"  -> DegLex
                     _         -> DegRevLex
      in go newOrd vars polys rest
    go ord _ polys ("--vars" : v : rest) =
      let newVars = splitCommas v
      in go ord newVars polys rest
    go ord vars polys (p : rest)
      | "--" `isPrefixOf` p = go ord vars polys rest
      | otherwise = go ord vars (p : polys) rest

    splitCommas [] = []
    splitCommas s =
      let (item, rest) = break (== ',') s
      in item : case rest of
                  [] -> []
                  (_:xs) -> splitCommas xs

-- | Batch execution: computes and prints reduced Groebner basis.
runGrobner :: [String] -> IO ()
runGrobner args = do
  let (ord, vars, polyStrs) = parseOptions args
  if null polyStrs
    then putStrLn "Error: No polynomials provided for Groebner basis computation." >> exitFailure
    else case mapM (parsePoly ord vars) polyStrs of
           Left err -> putStrLn ("Parse error: " ++ err) >> exitFailure
           Right (polys :: [Poly Rat]) -> do
             let gb = reducedGrobnerBasis (buchberger polys)
             putStrLn "\ESC[1mReduced Groebner Basis:\ESC[0m"
             mapM_ (\(i, g) -> putStrLn $ "  g_" ++ show (i :: Int) ++ " = " ++ prettyUnicode vars g) (zip [1..] gb)

-- | Batch execution: computes minimal resolution, Betti table, and homological invariants.
runBetti :: [String] -> IO ()
runBetti args = do
  let (ord, vars, polyStrs) = parseOptions args
      n = length vars
  if null polyStrs
    then putStrLn "Error: No polynomials provided for Betti table computation." >> exitFailure
    else case mapM (parsePoly ord vars) polyStrs of
           Left err -> putStrLn ("Parse error: " ++ err) >> exitFailure
           Right (polys :: [Poly Rat]) -> do
             let fullRes = freeResolution ord n polys
                 minRes = minimizeComplex ord fullRes
                 bt = computeBettiTable minRes
             putStrLn "\ESC[1mMinimal Graded Free Resolution Betti Table:\ESC[0m"
             putStr (renderBettiTableAnsi bt)
             putStr (renderBettiDiagnostics bt)

-- | Batch execution: raymarches a 3D algebraic variety.
runRaymarch :: [String] -> IO ()
runRaymarch args = do
  let name = case args of (x:_) -> map toLower x; [] -> "whitney"
      surface = case name of
        "clebsch" -> clebschCubic
        "cayley"  -> cayleyCubic
        "kummer"  -> kummerSurface
        "torus"   -> torusVariety 1.2 0.5
        "sphere"  -> sphereVariety 1.2
        _         -> whitneyUmbrella
  putStrLn $ "\ESC[1mRaymarching Algebraic Variety: " ++ name ++ "\ESC[0m"
  putStr (raymarchVariety3D surface 70 30)

-- | Batch execution: computes sheaf cohomology.
runSheaf :: [String] -> IO ()
runSheaf args = do
  let name = case args of (x:_) -> map toLower x; [] -> "torus"
      (label, complex) = case name of
        "sphere" -> ("2-Sphere S^2", sphere2Complex)
        "circle" -> ("Circle S^1", circleComplex 4)
        _        -> ("2-Torus T^2", torusComplex)
      sheaf = constantSheaf1D complex
      bettis = allSheafBettiNumbers sheaf
      chi = sheafEulerCharacteristic sheaf
  putStrLn $ "\ESC[1mCellular Sheaf Cohomology on " ++ label ++ ":\ESC[0m"
  putStrLn $ "  Betti numbers [b_0, b_1, ...]: " ++ show bettis
  putStrLn $ "  Euler characteristic chi      : " ++ show chi
