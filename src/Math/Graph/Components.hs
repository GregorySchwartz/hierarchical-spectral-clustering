{- Math.Graph.Components
Gregory W. Schwartz

Find connected components of matrices.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Graph.Components
  ( getComponentMats
  , getComponentMatsItems
  ) where

-- Remote
import Data.List (sort)
import qualified Data.Graph.Inductive as G
import qualified Data.Vector as V

-- Local
import Math.Graph.Types
import Math.Clustering.Hierarchical.Spectral.Utility

-- | Get the components of a graphable object, an adjacency matrix.
getComponentMats :: (Graphable a) => a -> [a]
getComponentMats mat = fmap (fromGraph . flip G.subgraph gr) . G.components $ gr
  where
    gr = toGraph mat

-- | Get the components of a graphable object, an adjacency matrix with the
-- associated items.
getComponentMatsItems :: (Graphable b) => V.Vector a -> b -> [(V.Vector a, b)]
getComponentMatsItems items mat =
    fmap (\ !xs -> (subsetVector items . sort $ xs, fromGraph . G.subgraph xs $ gr))
      . G.components
      $ gr
  where
    gr = toGraph mat
