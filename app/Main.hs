{- |
Module      : Main
Description : Main executable entry point for the Syzygy engine
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT
-}
module Main (main) where

import Syzygy.UI.CLI (mainCLI)

main :: IO ()
main = mainCLI
