{- Math.Graph.Types
Gregory W. Schwartz

Types used for graph sections of the program.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module Math.Graph.Types where

-- Remote
import Data.List (sort)
-- import qualified Data.Eigen.SparseMatrix as E
import qualified Data.Map.Strict as Map
import qualified Data.Graph.Inductive as G
import qualified Data.Sparse.Common as S
import qualified Data.Vector.Unboxed as V
import qualified Math.Clustering.Spectral.Dense as D
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Sparse as S

-- Local
import qualified Math.Clustering.Spectral.Dense as D
-- import qualified Math.Clustering.Spectral.Eigen.AdjacencyMatrix as E
import qualified Math.Clustering.Spectral.Sparse as S

-- | Get a re-mapped edge list with nodes ordered from 0 to the number of nodes
-- in the graph.
orderedEdges :: G.Gr Int a -> [(Int, Int, a)]
orderedEdges gr =
    fmap (\(!i,  !j,  !v) -> (getUpdate i, getUpdate j, v)) . G.labEdges $ gr
  where
    nodeMap = Map.fromList . flip zip [0..] . sort . G.nodes $ gr
    getUpdate = flip ( Map.findWithDefault
                        (error "Unexpected missing node in orderedEdges.")
                     )
                     nodeMap

-- | Graphable class for converting matrices to and from a graph.
class Graphable a where
  toGraph :: a -> G.Gr Int Double
  fromGraph :: G.Gr Int Double -> a

instance Graphable D.AdjacencyMatrix where
  toGraph mat = G.mkGraph (zip [0 .. H.rows mat - 1] [0 .. H.rows mat - 1])
              . filter (\(_, _, x) -> x /= 0)
              . concatMap (\(!i, !xs) -> fmap (\(!j, !v) -> (i, j, v)) xs)
              . zip [0..]
              . fmap (zip [0..] . H.toList)
              . H.toRows
              $ mat
  fromGraph gr = H.assoc (G.noNodes gr, G.noNodes gr) 0
               . fmap (\(!i, !j, !v) -> ((i, j), v))
               . orderedEdges
               $ gr

instance Graphable S.AdjacencyMatrix where
  toGraph mat = G.mkGraph (zip [0 .. S.nrows mat - 1] [0 .. S.nrows mat - 1])
              . filter (\(_, _, x) -> x /= 0)
              . S.toListSM
              $ mat
  fromGraph gr = S.fromListSM (G.noNodes gr, G.noNodes gr) . orderedEdges $ gr

-- instance Graphable E.AdjacencyMatrix where
--   toGraph mat = G.mkGraph (zip [0 .. E.rows mat - 1] [0 .. E.rows mat - 1])
--               . filter (\(_, _, x) -> x /= 0)
--               . E.toList
--               $ mat
--   fromGraph gr = E.fromList (G.noNodes gr) (G.noNodes gr) . orderedEdges $ gr
