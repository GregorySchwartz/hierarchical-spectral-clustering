{- Math.Clustering.Hierarchical.Spectral
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Dense
    ( hierarchicalSpectralCluster
    , getClusterItems
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Dense (spectralClusterNorm)
import Math.Modularity.Dense (getModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Foldable as F
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

type AdjacencyMatrix = H.Matrix Double
type Items a         = V.Vector a

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria.
hierarchicalSpectralCluster :: Items a
                            -> AdjacencyMatrix
                            -> ClusteringTree a AdjacencyMatrix
hierarchicalSpectralCluster !items !adjMat =
    if ngMod > Q 0 && H.rows adjMat > 1
        then
            Node { rootLabel = vertex
                 , subForest =
                    [ hierarchicalSpectralCluster (getItems leftIdxs) left
                    , hierarchicalSpectralCluster (getItems rightIdxs) right
                    ]
                 }
        else
            Node {rootLabel = vertex, subForest = []}
  where
    vertex      = ClusteringVertex { _clusteringItems = items
                                   , _clusteringMatrix = adjMat
                                   , _ngMod = ngMod
                                   }
    clusters    = spectralClusterNorm adjMat
    ngMod       = getModularity clusters $ adjMat
    getIdxs val = VS.ifoldl' (\ !acc !i -> bool acc (i:acc) . (== val)) []
    leftIdxs    = getIdxs 0 $ clusters
    rightIdxs   = getIdxs 1 $ clusters
    left        = adjMat H.?? (H.Pos (H.idxs leftIdxs), H.Pos (H.idxs leftIdxs))
    right       =
        adjMat H.?? (H.Pos (H.idxs rightIdxs), H.Pos (H.idxs rightIdxs))
    getItems    =
        V.fromList . F.foldl' (\ !acc !i -> (items V.! i) : acc) [] . V.fromList

-- | Gather clusters (leaves) from the tree.
getClusterItems :: Foldable t => t (Items a, b) -> [V.Vector a]
getClusterItems = fmap fst . F.toList
