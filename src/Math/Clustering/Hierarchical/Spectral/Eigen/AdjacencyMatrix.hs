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
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS

-- Local
import Math.Clustering.Hierarchical.Spectral.Types
import Math.Clustering.Hierarchical.Spectral.Utility

type AdjacencyMatrix = S.SparseMatrixXd
type Items a         = V.Vector a

-- | Check if there is more than one cluster.
hasMultipleClusters :: S.SparseMatrixXd -> Bool
hasMultipleClusters = (> 1)
                    . Set.size
                    . Set.fromList
                    . concat
                    . S.toDenseList

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can also use minimum number
-- of observations in a cluster as the stopping criteria.
hierarchicalSpectralCluster :: EigenGroup
                            -> Maybe NumEigen
                            -> Maybe Int
                            -> Maybe Q
                            -> Items a
                            -> AdjacencyMatrix
                            -> ClusteringTree a
hierarchicalSpectralCluster !eigenGroup !numEigenMay !minSizeMay !minModMay !items !adjMat =

    if S.rows adjMat > 1
        && hasMultipleClusters clusters
        && ngMod > minMod
        && S.rows left >= minSize
        && S.rows right >= minSize
        then do
            Node { rootLabel = vertex
                 , subForest = [ hierarchicalSpectralCluster
                                  eigenGroup
                                  numEigenMay
                                  minSizeMay
                                  minModMay
                                  (subsetVector items leftIdxs)
                                  left
                               , hierarchicalSpectralCluster
                                  eigenGroup
                                  numEigenMay
                                  minSizeMay
                                  minModMay (subsetVector items rightIdxs)
                                  right
                               ]
                 }
        else
            Node {rootLabel = vertex, subForest = []}
  where
    clusters = spectralClustering eigenGroup adjMat
    spectralClustering :: EigenGroup -> AdjacencyMatrix -> S.SparseMatrixXd
    spectralClustering SignGroup   = spectralClusterNorm
    spectralClustering KMeansGroup = spectralClusterKNorm numEigen 2
    minMod      = fromMaybe (Q 0) minModMay
    minSize     = fromMaybe 1 minSizeMay
    numEigen    = fromMaybe 1 numEigenMay
    vertex      = ClusteringVertex { _clusteringItems = items
                                   , _ngMod = ngMod
                                   }
    ngMod       = getModularity clusters adjMat
    getIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
                . VS.fromList
                . concat
                . S.toDenseList
    leftIdxs    = getIdxs 0 clusters
    rightIdxs   = getIdxs 1 clusters
    left        = S.squareSubset leftIdxs adjMat
    right       = S.squareSubset rightIdxs adjMat
