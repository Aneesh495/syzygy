{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- |
Module      : Syzygy.Sheaf.Complex
Description : Regular cell complexes, simplicial complexes, and boundary operators
Copyright   : (c) 2026 Aneesh Krishna
License     : MIT

Regular cellular and simplicial complexes form the topological substrate
for cellular sheaves and cosheaves.

A regular cell complex X is a graded poset of cells stratified by dimension:
  X = union_{k=0}^d X_k

where each k-cell sigma in X_k is homeomorphic to an open k-ball.
For codimension-1 incidence sigma in X_{k-1} and tau in X_k (written sigma <| tau),
the incidence number [sigma : tau] in (-1, +1) satisfies the fundamental boundary relation:

  sum_{tau in X_k} [sigma : tau] * [tau : gamma] = 0

for any (k-1)-cell sigma and (k+1)-cell gamma. This guarantees that the boundary
operator satisfies partial_k . partial_{k+1} = 0.
-}
module Syzygy.Sheaf.Complex
  ( Cell(..)
  , CellComplex(..)
  , emptyComplex
  , addCell
  , complexDimension
  , cellsOfDim
  , boundaryMatrix
  , verifyBoundarySquareZero
  , fromSimplices
  , circleComplex
  , sphere2Complex
  , torusComplex
  , projectivePlaneComplex
  ) where

import GHC.Generics (Generic)
import Data.List (sort, nub)
import qualified Data.Map.Strict as Map
import Syzygy.Linear.Matrix (Matrix(..), matMul, matZero)

-- | A cell in a regular cell complex or simplicial complex.
data Cell = Cell
  { cellId       :: !Int       -- ^ Unique integer identifier
  , cellDim      :: !Int       -- ^ Dimension k of the cell (0 for vertex, 1 for edge, etc.)
  , cellVertices :: ![Int]     -- ^ Sorted list of vertex IDs spanning the cell
  } deriving (Eq, Ord, Show, Read, Generic)

-- | A regular cell complex containing cells of various dimensions and signed facet incidences.
data CellComplex = CellComplex
  { ccCells      :: ![Cell]                         -- ^ All cells in the complex
  , ccIncidences :: !(Map.Map (Int, Int) Int)      -- ^ Map from (facetId, cellId) to sign in (-1 or +1)
  } deriving (Eq, Show, Read, Generic)

-- | The empty cell complex.
emptyComplex :: CellComplex
emptyComplex = CellComplex [] Map.empty

-- | Adds a cell to the complex.
addCell :: Cell -> CellComplex -> CellComplex
addCell c cc = cc { ccCells = c : ccCells cc }

-- | Maximum dimension of any cell in the complex.
complexDimension :: CellComplex -> Int
complexDimension cc =
  if null (ccCells cc)
  then -1
  else maximum (map cellDim (ccCells cc))

-- | Retrieves all cells of a specific dimension k, sorted by cellId.
cellsOfDim :: CellComplex -> Int -> [Cell]
cellsOfDim cc k =
  let matching = filter (\c -> cellDim c == k) (ccCells cc)
  in sort matching

-- | Constructs the boundary matrix partial_k: C_k(X; Q) -> C_{k-1}(X; Q).
-- Rows correspond to (k-1)-cells, columns correspond to k-cells.
boundaryMatrix :: CellComplex -> Int -> Matrix Rational
boundaryMatrix cc k =
  let kCells = cellsOfDim cc k
      km1Cells = cellsOfDim cc (k - 1)
      numCols = length kCells
      numRows = length km1Cells
      colIndices = Map.fromList (zip (map cellId kCells) [0 ..])
      rowIndices = Map.fromList (zip (map cellId km1Cells) [0 ..])
      -- Build matrix entries
      rowList rId =
        let rIdx = rowIndices Map.! rId
        in [ case Map.lookup (rId, cId) (ccIncidences cc) of
               Just s  -> (rIdx, colIndices Map.! cId, fromIntegral s)
               Nothing -> (rIdx, colIndices Map.! cId, 0)
           | cId <- map cellId kCells
           ]
      entryTriples = concatMap (rowList . cellId) km1Cells
      entryMap = Map.fromList [ ((r, c), v) | (r, c, v) <- entryTriples ]
      entries = [ [ Map.findWithDefault 0 (r, c) entryMap | c <- [0 .. numCols - 1] ]
                | r <- [0 .. numRows - 1]
                ]
  in Matrix numRows numCols entries

-- | Verifies that partial_k . partial_{k+1} = 0 as matrices over Q.
verifyBoundarySquareZero :: CellComplex -> Int -> Bool
verifyBoundarySquareZero cc k =
  let bk = boundaryMatrix cc k
      bkp1 = boundaryMatrix cc (k + 1)
      numR = case bk of Matrix r _ _ -> r
      numC = case bkp1 of Matrix _ c _ -> c
  in if numR == 0 || numC == 0
     then True
     else
       let prod = matMul bk bkp1
           zeroM = matZero numR numC
       in prod == zeroM

-- | Constructs a complete simplicial complex from a list of maximal simplices.
-- Automatically derives all lower-dimensional faces with canonical alternating signs:
-- [d_i(tau) : tau] = (-1)^i.
fromSimplices :: [[Int]] -> CellComplex
fromSimplices maxSimplices =
  let sortedSimps = map sort maxSimplices
      -- Subsets of size k
      subsets 0 _ = [[]]
      subsets _ [] = []
      subsets k (x:xs) = map (x:) (subsets (k - 1) xs) ++ subsets k xs

      allFaces = nub $ concat
        [ [ s | k <- [1 .. length simp], s <- subsets k simp ]
        | simp <- sortedSimps
        ]
      -- Group faces by dimension (length - 1)
      facesByDim = sort allFaces
      faceToId = Map.fromList (zip facesByDim [0 ..])
      cells = [ Cell (faceToId Map.! f) (length f - 1) f | f <- facesByDim ]

      -- Compute codimension-1 facet incidences with alternating signs
      buildIncidences face =
        let k = length face - 1
            fId = faceToId Map.! face
        in if k <= 0
           then []
           else
             [ let facet = take i face ++ drop (i + 1) face
                   facetId = faceToId Map.! facet
                   signVal = if even i then 1 else -1
               in ((facetId, fId), signVal)
             | i <- [0 .. k]
             ]
      incMap = Map.fromList (concatMap buildIncidences facesByDim)
  in CellComplex cells incMap

-- | Constructs the circle S^1 as an n-gon simplicial complex (default: triangle n = 3).
circleComplex :: Int -> CellComplex
circleComplex n =
  let count = max 3 n
      edges = [ [i, (i + 1) `mod` count] | i <- [0 .. count - 1] ]
  in fromSimplices edges

-- | Constructs the 2-sphere S^2 as an octahedron simplicial complex (6 vertices, 8 triangular faces).
sphere2Complex :: CellComplex
sphere2Complex =
  let faces =
        [ [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]
        , [5, 2, 1], [5, 3, 2], [5, 4, 3], [5, 1, 4]
        ]
  in fromSimplices faces

-- | Constructs the 2-torus T^2 as Moebius' minimal 7-vertex triangulation (14 triangles).
torusComplex :: CellComplex
torusComplex =
  let faces =
        [ [0, 1, 3], [1, 3, 4], [1, 2, 4], [2, 4, 5], [2, 0, 5], [0, 5, 3]
        , [3, 4, 6], [4, 5, 6], [5, 0, 6], [0, 1, 6], [1, 2, 6], [2, 3, 6]
        , [0, 2, 4], [1, 5, 3]
        ]
  in fromSimplices faces

-- | Constructs the real projective plane RP^2 as a minimal 6-vertex triangulation (10 triangles).
projectivePlaneComplex :: CellComplex
projectivePlaneComplex =
  let faces =
        [ [0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 0], [4, 0, 1]
        , [0, 2, 5], [2, 4, 5], [4, 1, 5], [1, 3, 5], [3, 0, 5]
        ]
  in fromSimplices faces
