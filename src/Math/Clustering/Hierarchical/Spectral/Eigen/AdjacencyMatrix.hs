{- Math.Clustering.Hierarchical.Spectral.Eigen.AdjacencyMatrix
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Eigen.AdjacencyMatrix
    ( hierarchicalSpectralCluster
    , AdjacencyMatrix (..)
    , Items (..)
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Maybe (fromMaybe)
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Eigen.AdjacencyMatrix (spectralClusterNorm, spectralClusterKNorm)
import Math.Modularity.Eigen.Sparse (getModularity)
import Math.Modularity.Types (Q (..))
import Safe (headMay)
import qualified Data.Foldable as F
import qualified Data.Set as Set
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Eigen.SparseMatrix.Utility as S
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

type AdjacencyMatrix = S.SparseMatrixXd
type Items a         = V.Vector a

-- | Check if there is more than one cluster.
hasMultipleClusters :: S.SparseMatrixXd -> Bool
hasMultipleClusters = (> 1)
                    . Set.size
                    . Set.fromList
                    . fromMaybe (error "No rows in \"vector\".")
                    . headMay
                    . S.toDenseList

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
        && S.rows adjMat > 1
        && S.rows left >= minSize
        && S.rows right >= minSize
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
    spectralClustering :: EigenGroup -> AdjacencyMatrix -> S.SparseMatrixXd
    spectralClustering SignGroup   = spectralClusterNorm
    spectralClustering KMeansGroup = spectralClusterKNorm 2
    ngMod       = getModularity clusters $ adjMat
    getIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
                . VS.fromList
                . fromMaybe (error "No rows in \"vector\".")
                . headMay
                . S.toDenseList
    leftIdxs    = getIdxs 0 clusters
    rightIdxs   = getIdxs 1 clusters
    left        = S.squareSubset leftIdxs adjMat
    right       = S.squareSubset rightIdxs adjMat
    getItems    = V.fromList . F.foldr' (\ !i !acc -> (items V.! i) : acc) []
