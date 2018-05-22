{- Math.Clustering.Hierarchical.Spectral.Eigen.FeatureMatrix
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering for
feature matrices.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Eigen.FeatureMatrix
    ( hierarchicalSpectralCluster
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
import Math.Clustering.Spectral.Eigen.FeatureMatrix (B (..), getB, spectralCluster, spectralClusterK)
import Math.Modularity.Eigen.Sparse (getBModularity)
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

type FeatureMatrix   = S.SparseMatrixXd
type Items a         = V.Vector a
type ShowB           = ((Int, Int), [(Int, Int, Double)])
type NormalizeFlag   = Bool

-- | Check if there is more than one cluster.
hasMultipleClusters :: S.SparseMatrixXd -> Bool
hasMultipleClusters = (> 1)
                    . Set.size
                    . Set.fromList
                    . fromMaybe (error "No rows in \"vector\".")
                    . headMay
                    . S.toDenseList

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
            && (S.rows $ unB b) > 1
            && S.rows (unB left) >= minSize
            && S.rows (unB right) >= minSize
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
        clusters :: S.SparseMatrixXd
        clusters = spectralClustering eigenGroup b
        spectralClustering :: EigenGroup -> B -> S.SparseMatrixXd
        spectralClustering SignGroup   = spectralCluster
        spectralClustering KMeansGroup = spectralClusterK 2
        ngMod :: Q
        ngMod       = getBModularity clusters $ b
        getSortedIdxs :: Double -> S.SparseMatrixXd -> [Int]
        getSortedIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
                          . VS.fromList
                          . fromMaybe (error "No rows in \"vector\".")
                          . headMay
                          . S.toDenseList
        leftIdxs :: [Int]
        leftIdxs    = getSortedIdxs 0 $ clusters
        rightIdxs :: [Int]
        rightIdxs   = getSortedIdxs 1 $ clusters
        left :: B
        left        = B $ extractRows (unB b) leftIdxs
        right :: B
        right       = B $ extractRows (unB b) rightIdxs
        extractRows :: S.SparseMatrixXd -> [Int] -> S.SparseMatrixXd
        extractRows mat = S.fromRows . fmap (flip S.getRow mat)
        getItems    =
            V.fromList . F.foldr' (\ !i !acc -> (items V.! i) : acc) []
