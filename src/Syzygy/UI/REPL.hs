{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

{- |
Module      : Syzygy.UI.REPL
Description : Interactive algebraic geometry Read-Eval-Print Loop
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Interactive REPL providing real-time Gröbner basis computation, graded free
resolutions, Betti table visualization, Hilbert series analysis, algebraic surface
raymarching, and cellular sheaf cohomology.
-}
module Syzygy.UI.REPL
  ( runREPL
  , ReplState(..)
  , initialReplState
  ) where

import Data.Char (toLower)
import Data.List (isPrefixOf)
import System.IO (hFlush, stdout)
import Syzygy.Core.Rational (Rat)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Polynomial.Base (Poly(..))
import Syzygy.Polynomial.Parser (parsePoly, defaultVarNames)
import Syzygy.Polynomial.Pretty (prettyUnicode)
import Syzygy.Grobner.Reduce (reducedGrobnerBasis)
import Syzygy.Resolution.FreeResolution (freeResolution)
import Syzygy.Resolution.Minimization (minimizeComplex)
import Syzygy.Resolution.Betti (computeBettiTable)
import Syzygy.Resolution.HilbertSeries (hilbertSeriesFromResolution, formatHilbertSeries)
import Syzygy.Geometry.Dimension (idealDimension)
import Syzygy.UI.BettiTable (renderBettiTableAnsi)
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
import Syzygy.Sheaf.Cohomology (allSheafBettiNumbers)

-- | State maintained across the interactive session.
data ReplState = ReplState
  { replOrder    :: !MonomialOrder
  , replVarNames :: ![String]
  , replIdeal    :: ![Poly Rat]
  }

-- | Default initial state: Grevlex order, vars x, y, z, w, empty ideal.
initialReplState :: ReplState
initialReplState = ReplState
  { replOrder    = DegRevLex
  , replVarNames = take 4 defaultVarNames
  , replIdeal    = []
  }

-- | Welcome ASCII banner printed on startup.
welcomeBanner :: String
welcomeBanner = unlines
  [ "\ESC[1;36m"
  , "   _____                           "
  , "  / ____|                          "
  , " | (___  _   _ _____   _  __ _ _   _ "
  , "  \\___ \\| | | |_  / | | |/ _` | | | |"
  , "  ____) | |_| |/ /| |_| | (_| | |_| |"
  , " |_____/ \\__, /___|\\__, |\\__, |\\__, |"
  , "          __/ |     __/ | __/ | __/ |"
  , "         |___/     |___/ |___/ |___/ "
  , "\ESC[0m"
  , "\ESC[1mComputational Algebraic Geometry & Graded Free Resolutions\ESC[0m"
  , "Type \ESC[33m:help\ESC[0m for command reference, or enter polynomials directly."
  , ""
  ]

-- | Help documentation string.
helpText :: String
helpText = unlines
  [ "=== Syzygy REPL Commands ==="
  , "  :order <lex|deglex|grevlex>  Set monomial ordering (current: lex/grevlex)"
  , "  :vars <v1> <v2> ...          Set coordinate variable names"
  , "  :ideal <poly1>, <poly2> ...  Set generators of the current ideal"
  , "  :add <poly>                  Add a polynomial generator to current ideal"
  , "  :clear                       Clear current ideal"
  , "  :show                        Print current generators"
  , "  :gb                          Compute reduced Groebner basis"
  , "  :res                         Compute minimal graded free resolution"
  , "  :betti                       Display ANSI graded Betti table"
  , "  :hilbert                     Compute Hilbert-Poincare series"
  , "  :dim                         Compute Krull dimension of V(I)"
  , "  :raymarch <surface>          Raymarch 3D variety: whitney|clebsch|cayley|kummer|torus|sphere"
  , "  :sheaf <space>               Compute sheaf cohomology: circle|sphere|torus"
  , "  :quit, :q                    Exit the interactive session"
  ]

-- | Main REPL entry point.
runREPL :: IO ()
runREPL = do
  putStr welcomeBanner
  loop initialReplState
  where
    loop state = do
      putStr "\ESC[1;32msyzygy>\ESC[0m "
      hFlush stdout
      line <- getLine
      let trimmed = dropWhile (== ' ') line
      if null trimmed
        then loop state
        else case trimmed of
          ":quit" -> putStrLn "Exiting Syzygy."
          ":q"    -> putStrLn "Exiting Syzygy."
          ":help" -> putStr helpText >> loop state
          ":?"    -> putStr helpText >> loop state
          cmd | ":" `isPrefixOf` cmd -> handleCommand state (drop 1 cmd) >>= loop
          polyStr -> handleAddPoly state polyStr >>= loop

-- | Command handler for colon directives.
handleCommand :: ReplState -> String -> IO ReplState
handleCommand state cmdLine =
  let (cmd, rest) = break (== ' ') cmdLine
      arg = dropWhile (== ' ') rest
  in case map toLower cmd of
       "order" ->
         let newOrd = case map toLower arg of
                        "lex"     -> Lex
                        "deglex"  -> DegLex
                        _         -> DegRevLex
         in do
           putStrLn $ "Monomial order set to " ++ show newOrd
           return state { replOrder = newOrd }

       "vars" ->
         let vars = words arg
         in if null vars
            then do
              putStrLn $ "Variables: " ++ unwords (replVarNames state)
              return state
            else do
              putStrLn $ "Variables set to: " ++ unwords vars
              return state { replVarNames = vars }

       "clear" -> do
         putStrLn "Current ideal cleared."
         return state { replIdeal = [] }

       "show" -> do
         let vars = replVarNames state
         putStrLn "Current generators:"
         mapM_ (\(i, p) -> putStrLn $ "  f_" ++ show (i :: Int) ++ " = " ++ prettyUnicode vars p) (zip [1..] (replIdeal state))
         return state

       "gb" -> do
         let gens = replIdeal state
         if null gens
           then putStrLn "Current ideal is empty."
           else do
             let gb = reducedGrobnerBasis gens
                 vars = replVarNames state
             putStrLn "\ESC[1mReduced Groebner Basis:\ESC[0m"
             mapM_ (\(i, g) -> putStrLn $ "  g_" ++ show (i :: Int) ++ " = " ++ prettyUnicode vars g) (zip [1..] (gb :: [Poly Rat]))
         return state

       "res" -> do
         let gens = replIdeal state
             ord = replOrder state
             n = length (replVarNames state)
         if null gens
           then putStrLn "Current ideal is empty."
           else do
             let fullRes = freeResolution ord n gens
                 minRes = minimizeComplex ord fullRes
                 bt = computeBettiTable minRes
             putStrLn "\ESC[1mMinimal Graded Free Resolution:\ESC[0m"
             putStrLn (renderBettiTableAnsi bt)
         return state

       "betti" -> do
         let gens = replIdeal state
             ord = replOrder state
             n = length (replVarNames state)
         if null gens
           then putStrLn "Current ideal is empty."
           else do
             let fullRes = freeResolution ord n gens
                 minRes = minimizeComplex ord fullRes
                 bt = computeBettiTable minRes
             putStrLn (renderBettiTableAnsi bt)
         return state

       "hilbert" -> do
         let gens = replIdeal state
             ord = replOrder state
             n = length (replVarNames state)
         if null gens
           then putStrLn "Current ideal is empty."
           else do
             let fullRes = freeResolution ord n gens
                 minRes = minimizeComplex ord fullRes
                 hs = hilbertSeriesFromResolution minRes n
             putStrLn "\ESC[1mHilbert-Poincare Series:\ESC[0m"
             putStrLn (formatHilbertSeries hs)
         return state

       "dim" -> do
         let gens = replIdeal state
             n = length (replVarNames state)
         if null gens
           then putStrLn "Dimension of empty ideal: ambient space"
           else do
             let gb = reducedGrobnerBasis gens
                 d = idealDimension n gb
             putStrLn $ "Krull dimension dim(V(I)) = " ++ show d
         return state

       "raymarch" -> do
         let surface = case map toLower arg of
               "clebsch" -> clebschCubic
               "cayley"  -> cayleyCubic
               "kummer"  -> kummerSurface
               "torus"   -> torusVariety 1.2 0.5
               "sphere"  -> sphereVariety 1.2
               _         -> whitneyUmbrella
         putStrLn "\ESC[1mRaymarching Algebraic Variety (Taubin sphere tracing):\ESC[0m"
         putStrLn (raymarchVariety3D surface 64 28)
         return state

       "sheaf" -> do
         let (name, complex) = case map toLower arg of
               "sphere" -> ("2-Sphere S^2 (Octahedron)", sphere2Complex)
               "torus"  -> ("2-Torus T^2 (Moebius Triangulation)", torusComplex)
               _        -> ("Circle S^1", circleComplex 4)
             sheaf = constantSheaf1D complex
             bettis = allSheafBettiNumbers sheaf
         putStrLn $ "\ESC[1mCellular Sheaf Cohomology on " ++ name ++ ":\ESC[0m"
         putStrLn $ "  Betti numbers [b_0, b_1, ...]: " ++ show bettis
         return state

       unknown -> do
         putStrLn $ "Unknown command: :" ++ unknown ++ ". Type :help for list of commands."
         return state

-- | Adds a parsed polynomial to the current ideal.
handleAddPoly :: ReplState -> String -> IO ReplState
handleAddPoly state input =
  let vars = replVarNames state
      ord = replOrder state
  in case parsePoly ord vars input of
       Left err -> do
         putStrLn $ "\ESC[31mParse error:\ESC[0m " ++ err
         return state
       Right p -> do
         let newIdeal = replIdeal state ++ [p]
             idx = length newIdeal
         putStrLn $ "Added generator f_" ++ show idx ++ " = " ++ prettyUnicode vars p
         return state { replIdeal = newIdeal }
