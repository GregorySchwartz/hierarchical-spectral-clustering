{- Math.Clustering.Hierarchical.Spectral.Sparse
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering for sparse
data.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Sparse
    ( hierarchicalSpectralCluster
    , hierarchicalSpectralClusterAdj
    , FeatureMatrix (..)
    , B (..)
    , Items (..)
    , ShowB (..)
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Maybe (fromMaybe)
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Sparse (B (..), AdjacencyMatrix (..), getB, spectralCluster, spectralClusterK, spectralClusterNorm, spectralClusterKNorm)
import Math.Modularity.Sparse (getBModularity, getModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Foldable as F
import qualified Data.Set as Set
import qualified Data.Sparse.Common as S
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra.Sparse as S

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

type FeatureMatrix   = S.SpMatrix Double
type Items a         = V.Vector a
type ShowB           = ((Int, Int), [(Int, Int, Double)])
type NormalizeFlag   = Bool

-- | Check if there is more than one cluster.
hasMultipleClusters :: S.SpVector Double -> Bool
hasMultipleClusters = (> 1) . Set.size . Set.fromList . S.toDenseListSV

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can use minimum number of
-- observations in a cluster as a stopping criteria. Assumes the feature matrix
-- has column features and row observations. Items correspond to rows. Can
-- use FeatureMatrix or a pre-generated B matrix. See Shu et al., "Efficient
-- Spectral Neighborhood Blocking for Entity Resolution", 2011.
hierarchicalSpectralCluster :: EigenGroup
                            -> NormalizeFlag
                            -> Maybe Int
                            -> Items a
                            -> Either FeatureMatrix B
                            -> ClusteringTree a
hierarchicalSpectralCluster eigenGroup normFlag initMinSizeMay initItems initMat =
    go initMinSizeMay initItems initB
  where
    initB = either (getB normFlag) id $ initMat
    go :: Maybe Int -> Items a -> B -> ClusteringTree a
    go !minSizeMay !items !b =
        if ngMod > Q 0
            && (S.nrows $ unB b) > 1
            && S.nrows (unB left) >= minSize
            && S.nrows (unB right) >= minSize
            && hasMultipleClusters clusters
            then
                Node { rootLabel = vertex
                     , subForest = [ go minSizeMay (getItems leftIdxs) left
                                   , go minSizeMay (getItems rightIdxs) right
                                   ]
                     }

            else
                Node {rootLabel = vertex, subForest = []}
      where
        minSize     = fromMaybe 1 minSizeMay
        vertex      = ClusteringVertex
                        { _clusteringItems = items
                        , _ngMod = ngMod
                        }
        clusters :: S.SpVector Double
        clusters = spectralClustering eigenGroup b
        spectralClustering :: EigenGroup -> B -> S.SpVector Double
        spectralClustering SignGroup   = spectralCluster
        spectralClustering KMeansGroup = spectralClusterK 2
        ngMod :: Q
        ngMod       = getBModularity clusters $ b
        getSortedIdxs :: Double -> S.SpVector Double -> [Int]
        getSortedIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
                          . VS.fromList
                          . S.toDenseListSV
        leftIdxs :: [Int]
        leftIdxs    = getSortedIdxs 0 $ clusters
        rightIdxs :: [Int]
        rightIdxs   = getSortedIdxs 1 $ clusters
        left :: B
        left        = B $ extractRows (unB b) leftIdxs
        right :: B
        right       = B $ extractRows (unB b) rightIdxs
        extractRows :: S.SpMatrix Double -> [Int] -> S.SpMatrix Double
        extractRows mat = S.transposeSM . S.fromColsL . fmap (S.extractRow mat)
        getItems    =
            V.fromList
                . F.foldr' (\ !i !acc
                           -> ( fromMaybe (error "Matrix has more rows than items.")
                              $ items V.!? i
                              ) : acc
                           ) []

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can use minimum number of
-- observations in a cluster as a stopping criteria. Uses an adjacency matrix.
-- Items correspond to rows.
hierarchicalSpectralClusterAdj :: EigenGroup
                               -> Maybe Int
                               -> Items a
                               -> AdjacencyMatrix
                               -> ClusteringTree a
hierarchicalSpectralClusterAdj eigenGroup initMinSizeMay initItems initMat =
    go initMinSizeMay initItems initMat
  where
    go :: Maybe Int -> Items a -> AdjacencyMatrix -> ClusteringTree a
    go !minSizeMay !items !mat =
        if ngMod > Q 0
            && (S.nrows $ mat) > 1
            && S.nrows mat >= minSize
            && S.nrows mat >= minSize
            && hasMultipleClusters clusters
            then
                Node { rootLabel = vertex
                     , subForest = [ go minSizeMay (getItems leftIdxs) left
                                   , go minSizeMay (getItems rightIdxs) right
                                   ]
                     }

            else
                Node {rootLabel = vertex, subForest = []}
      where
        minSize     = fromMaybe 1 minSizeMay
        vertex      = ClusteringVertex
                        { _clusteringItems = items
                        , _ngMod = ngMod
                        }
        clusters :: S.SpVector Double
        clusters = spectralClustering eigenGroup mat
        spectralClustering :: EigenGroup -> AdjacencyMatrix -> S.SpVector Double
        spectralClustering SignGroup   = spectralClusterNorm
        spectralClustering KMeansGroup = spectralClusterKNorm 2
        ngMod :: Q
        ngMod = getModularity clusters mat
        getSortedIdxs :: Double -> S.SpVector Double -> [Int]
        getSortedIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
                          . VS.fromList
                          . S.toDenseListSV
        leftIdxs :: [Int]
        leftIdxs    = getSortedIdxs 0 $ clusters
        rightIdxs :: [Int]
        rightIdxs   = getSortedIdxs 1 $ clusters
        left :: AdjacencyMatrix
        left        = getSubMat mat leftIdxs
        right :: AdjacencyMatrix
        right       = getSubMat mat rightIdxs
        getSubMat :: S.SpMatrix Double -> [Int] -> S.SpMatrix Double
        getSubMat mat is = S.fromColsL
                         . (\x -> fmap (S.extractCol x) is)
                         . S.transposeSM
                         . S.fromColsL
                         . fmap (S.extractRow mat)
                         $ is
        getItems    =
            V.fromList
                . F.foldr' (\ !i !acc
                           -> ( fromMaybe (error "Matrix has more rows than items.")
                              $ items V.!? i
                              ) : acc
                           ) []
