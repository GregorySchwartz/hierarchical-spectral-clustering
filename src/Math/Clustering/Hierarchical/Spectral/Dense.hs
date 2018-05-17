{- Math.Clustering.Hierarchical.Spectral
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Dense
    ( hierarchicalSpectralCluster
    , AdjacencyMatrix (..)
    , Items (..)
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Maybe (fromMaybe)
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Dense (spectralClusterNorm, spectralClusterKNorm)
import Math.Modularity.Dense (getModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Foldable as F
import qualified Data.Set as Set
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

type AdjacencyMatrix = H.Matrix Double
type Items a         = V.Vector a

-- | Check if there is more than one cluster.
hasMultipleClusters :: H.Vector Double -> Bool
hasMultipleClusters = (> 1) . Set.size . Set.fromList . H.toList

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can also use minimum number
-- of observations in a cluster as the stopping criteria.
hierarchicalSpectralCluster :: EigenGroup
                            -> Maybe Int
                            -> Items a
                            -> AdjacencyMatrix
                            -> ClusteringTree a
hierarchicalSpectralCluster eigenGroup !minSizeMay !items !adjMat =
    if ngMod > Q 0
        && H.rows adjMat > 1
        && H.rows left >= minSize
        && H.rows right >= minSize
        && hasMultipleClusters clusters
        then
            Node { rootLabel = vertex
                 , subForest =
                    [ hierarchicalSpectralCluster eigenGroup minSizeMay (getItems leftIdxs) left
                    , hierarchicalSpectralCluster eigenGroup minSizeMay (getItems rightIdxs) right
                    ]
                 }
        else
            Node {rootLabel = vertex, subForest = []}
  where
    minSize     = fromMaybe 1 minSizeMay
    vertex      = ClusteringVertex { _clusteringItems = items
                                   , _ngMod = ngMod
                                   }
    clusters    = spectralClustering eigenGroup adjMat
    spectralClustering :: EigenGroup -> AdjacencyMatrix -> H.Vector Double
    spectralClustering SignGroup   = spectralClusterNorm
    spectralClustering KMeansGroup = spectralClusterKNorm 2
    ngMod       = getModularity clusters $ adjMat
    getIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
    leftIdxs    = getIdxs 0 $ clusters
    rightIdxs   = getIdxs 1 $ clusters
    left        = adjMat H.?? (H.Pos (H.idxs leftIdxs), H.Pos (H.idxs leftIdxs))
    right       =
        adjMat H.?? (H.Pos (H.idxs rightIdxs), H.Pos (H.idxs rightIdxs))
    getItems    =
        V.fromList . F.foldr' (\ !i !acc -> (items V.! i) : acc) [] . V.fromList
