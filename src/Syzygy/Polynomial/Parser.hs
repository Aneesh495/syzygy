{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Syzygy.Polynomial.Parser
  ( parsePoly
  , parsePolyDefault
  , parseSystem
  , defaultVarNames
  ) where

import Data.Char (isAlpha, isAlphaNum, isDigit, isSpace)
import qualified Data.Map.Strict as Map
import Syzygy.Core.Rational (Rat, rat)
import Syzygy.Monomial.Order (MonomialOrder(..))
import Syzygy.Monomial.Sparse (sparseVar)
import Syzygy.Polynomial.Base
  ( Term(..)
  , Poly(..)
  , constantPoly
  )
import Syzygy.Polynomial.Arithmetic (addPoly, subPoly, mulPoly, scalePoly, powPoly)

-- | Default variable names used when none are explicitly provided.
defaultVarNames :: [String]
defaultVarNames = ["x", "y", "z", "w", "u", "v", "a", "b", "c", "d"]

-- | Token representation for polynomial expressions.
data Token
  = TokNum !Integer !(Maybe Integer) -- ^ Number n or fraction n/d
  | TokVar !String                  -- ^ Variable name identifier
  | TokPlus                         -- ^ +
  | TokMinus                        -- ^ -
  | TokStar                         -- ^ *
  | TokCaret                        -- ^ ^
  | TokLParen                       -- ^ (
  | TokRParen                       -- ^ )
  | TokEOF
  deriving (Eq, Show)

-- | Lexer: converts an input string into a list of tokens.
tokenize :: String -> Either String [Token]
tokenize [] = Right [TokEOF]
tokenize (c:cs)
  | isSpace c = tokenize cs
  | c == '+'  = (TokPlus :) <$> tokenize cs
  | c == '-'  = (TokMinus :) <$> tokenize cs
  | c == '*'  = (TokStar :) <$> tokenize cs
  | c == '^'  = (TokCaret :) <$> tokenize cs
  | c == '('  = (TokLParen :) <$> tokenize cs
  | c == ')'  = (TokRParen :) <$> tokenize cs
  | isDigit c =
      let (digits, rest) = span isDigit (c:cs)
          num = read digits :: Integer
      in case rest of
           ('/' : d : ds) | isDigit d ->
             let (denDigits, rest2) = span isDigit (d:ds)
                 den = read denDigits :: Integer
             in if den == 0
                then Left "Division by zero in rational constant"
                else (TokNum num (Just den) :) <$> tokenize rest2
           _ -> (TokNum num Nothing :) <$> tokenize rest
  | isAlpha c =
      let (name, rest) = span isAlphaNum (c:cs)
      in (TokVar name :) <$> tokenize rest
  | otherwise = Left $ "Unexpected character: " ++ [c]

-- | Recursive descent parser for polynomial expressions:
-- Expr   ::= Term (('+' | '-') Term)*
-- Term   ::= Factor ('*' Factor | Factor)*
-- Factor ::= Power
-- Power  ::= Atom ('^' Integer)?
-- Atom   ::= Number | Variable | '(' Expr ')'
parseExpression :: MonomialOrder -> Map.Map String Int -> [Token] -> Either String (Poly Rat, [Token])
parseExpression ord varMap tokens = parseExpr tokens
  where
    parseExpr toks = do
      (left, toks1) <- parseTerm toks
      go left toks1
      where
        go !acc (TokPlus : ts) = do
          (right, ts') <- parseTerm ts
          go (addPoly acc right) ts'
        go !acc (TokMinus : ts) = do
          (right, ts') <- parseTerm ts
          go (subPoly acc right) ts'
        go !acc ts = Right (acc, ts)

    parseTerm toks = do
      (firstFactor, toks1) <- parseFactor toks
      go firstFactor toks1
      where
        go !acc (TokStar : ts) = do
          (factor, ts') <- parseFactor ts
          go (mulPoly acc factor) ts'
        go !acc (t@(TokVar _) : _) = do
          -- Implicit multiplication: 2*x or x*y without explicit star
          (factor, ts') <- parseFactor (t : drop 1 toks)
          go (mulPoly acc factor) ts'
        go !acc (TokLParen : ts) = do
          (factor, ts') <- parseFactor (TokLParen : ts)
          go (mulPoly acc factor) ts'
        go !acc ts = Right (acc, ts)

    parseFactor toks = do
      (base, toks1) <- parseAtom toks
      case toks1 of
        (TokCaret : TokNum expVal Nothing : ts2) ->
          Right (powPoly base expVal, ts2)
        (TokCaret : _) ->
          Left "Exponent must be a non-negative integer"
        _ -> Right (base, toks1)

    parseAtom (TokNum n Nothing : ts) =
      Right (constantPoly (rat n 1) ord, ts)
    parseAtom (TokNum n (Just d) : ts) =
      Right (constantPoly (rat n d) ord, ts)
    parseAtom (TokVar name : ts) =
      case Map.lookup name varMap of
        Just varId ->
          let p = Poly [Term 1 (sparseVar varId 1)] ord
          in Right (p, ts)
        Nothing -> Left $ "Unknown variable: " ++ name
    parseAtom (TokLParen : ts) = do
      (inner, toks1) <- parseExpr ts
      case toks1 of
        (TokRParen : ts2) -> Right (inner, ts2)
        _                 -> Left "Expected closing parenthesis ')'"
    parseAtom (TokMinus : ts) = do
      (atom, ts1) <- parseAtom ts
      Right (scalePoly (-1) atom, ts1)
    parseAtom (tok : _) = Left $ "Unexpected token in atom: " ++ show tok
    parseAtom [] = Left "Unexpected end of expression"

-- | Parses a multivariate polynomial expression string with explicit variable mapping.
parsePoly :: MonomialOrder -> [String] -> String -> Either String (Poly Rat)
parsePoly ord varNames input = do
  tokens <- tokenize input
  let varMap = Map.fromList (zip varNames [0 :: Int ..])
  (polyVal, rest) <- parseExpression ord varMap tokens
  case rest of
    [TokEOF] -> Right polyVal
    []       -> Right polyVal
    _        -> Left "Trailing unparsed tokens at end of expression"

-- | Parses a polynomial using default variables ["x", "y", "z", "w", "u", "v"] and DegRevLex order.
parsePolyDefault :: String -> Either String (Poly Rat)
parsePolyDefault = parsePoly DegRevLex defaultVarNames

-- | Parses a system of polynomial equations/expressions.
parseSystem :: [String] -> [String] -> Either String [Poly Rat]
parseSystem varNames inputs =
  mapM (parsePoly DegRevLex varNames) inputs
