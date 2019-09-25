{- Math.Clustering.Hierarchical.Spectral
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Dense
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
import Math.Clustering.Spectral.Dense (B (..), AdjacencyMatrix (..), getB, spectralCluster, spectralClusterK, spectralClusterNorm, spectralClusterKNorm)
import Math.Modularity.Dense (getModularity, getBModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Foldable as F
import qualified Data.Set as Set
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Clustering.Hierarchical.Spectral.Types
import Math.Clustering.Hierarchical.Spectral.Utility

type FeatureMatrix   = H.Matrix Double
type Items a         = V.Vector a
type ShowB           = ((Int, Int), [(Int, Int, Double)])
type NormalizeFlag   = Bool

-- | Check if there is more than one cluster.
hasMultipleClusters :: H.Vector Double -> Bool
hasMultipleClusters = (> 1) . Set.size . Set.fromList . H.toList

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can use minimum number of
-- observations in a cluster as a stopping criteria. Assumes the feature matrix
-- has column features and row observations. Items correspond to rows. Can
-- use FeatureMatrix or a pre-generated B matrix. See Shu et al., "Efficient
-- Spectral Neighborhood Blocking for Entity Resolution", 2011.
hierarchicalSpectralCluster :: EigenGroup
                            -> NormalizeFlag
                            -> Maybe NumEigen
                            -> Maybe Int
                            -> Maybe Q
                            -> Maybe Runs
                            -> Items a
                            -> Either FeatureMatrix B
                            -> IO (ClusteringTree a)
hierarchicalSpectralCluster eigenGroup normFlag numEigenMay minSizeMay minModMay runsMay initItems initMat =
    go initItems initB
  where
    initB = either (getB normFlag) id $ initMat
    minMod   = fromMaybe (Q 0) minModMay
    minSize  = fromMaybe 1 minSizeMay
    numEigen = fromMaybe 1 numEigenMay
    go :: Items a -> B -> IO (ClusteringTree a)
    go !items !b =
        if (H.rows $ unB b) > 1
            && hasMultipleClusters clusters
            && ngMod > minMod
            && H.rows (unB left) >= minSize
            && H.rows (unB right) >= minSize
            then do
                goLeft <- go (subsetVector items leftIdxs) left
                goRight <- go (subsetVector items rightIdxs) right
                vertex <- getVertex
                return $ Node { rootLabel = vertex
                              , subForest = [goLeft, goRight]
                              }
            else do
                vertex <- getVertex
                return $ Node {rootLabel = vertex, subForest = []}
      where
        getVertex = do
          pVal <- mapM (\ x -> permutationTest
                                (\l m -> unQ $ getBModularity l m)
                                x
                                clusters
                                b
                       )
                  runsMay
          return $ ClusteringVertex
                 { _clusteringItems = items
                 , _ngMod = ngMod
                 , _pVal = pVal
                 }
        clusters :: H.Vector Double
        clusters = spectralClustering eigenGroup b
        spectralClustering :: EigenGroup -> B -> H.Vector Double
        spectralClustering SignGroup   = spectralCluster
        spectralClustering KMeansGroup = spectralClusterK numEigen 2
        ngMod :: Q
        ngMod       = getBModularity clusters $ b
        getIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
        leftIdxs    = getIdxs 0 $ clusters
        rightIdxs   = getIdxs 1 $ clusters
        left        = B $ (unB b) H.? leftIdxs
        right       = B $ (unB b) H.? rightIdxs

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Can also use minimum number
-- of observations in a cluster as the stopping criteria.
hierarchicalSpectralClusterAdj :: (Show a) => EigenGroup
                            -> Maybe NumEigen
                            -> Maybe Int
                            -> Maybe Q
                            -> Maybe Runs
                            -> Items a
                            -> AdjacencyMatrix
                            -> IO (ClusteringTree a)
hierarchicalSpectralClusterAdj !eigenGroup !numEigenMay !minSizeMay !minModMay !runsMay !items !adjMat =
    if H.rows adjMat > 1
        && hasMultipleClusters clusters
        && ngMod > minMod
        && H.rows left >= minSize
        && H.rows right >= minSize
        then do
            goLeft <- hierarchicalSpectralClusterAdj
                        eigenGroup
                        numEigenMay
                        minSizeMay
                        minModMay
                        runsMay
                        (subsetVector items leftIdxs)
                        left
            goRight <- hierarchicalSpectralClusterAdj
                        eigenGroup
                        numEigenMay
                        minSizeMay
                        minModMay
                        runsMay
                        (subsetVector items rightIdxs)
                        right
            vertex <- getVertex
            return $ Node { rootLabel = vertex , subForest = [goLeft, goRight] }
        else do
            vertex <- getVertex
            return $ Node {rootLabel = vertex, subForest = []}
  where
    minMod      = fromMaybe (Q 0) minModMay
    minSize     = fromMaybe 1 minSizeMay
    numEigen    = fromMaybe 1 numEigenMay
    getVertex = do
      pVal <- mapM
                (\ x -> permutationTest
                        (\l m -> unQ $ getModularity l m)
                        x
                        clusters
                        adjMat
                )
                runsMay
      return $ ClusteringVertex
             { _clusteringItems = items
             , _ngMod = ngMod
             , _pVal = pVal
             }
    clusters    = spectralClustering eigenGroup adjMat
    spectralClustering :: EigenGroup -> AdjacencyMatrix -> H.Vector Double
    spectralClustering SignGroup   = spectralClusterNorm
    spectralClustering KMeansGroup = spectralClusterKNorm numEigen 2
    ngMod       = getModularity clusters $ adjMat
    getIdxs val = VS.ifoldr' (\ !i !v !acc -> bool acc (i:acc) $ v == val) []
    leftIdxs    = getIdxs 0 $ clusters
    rightIdxs   = getIdxs 1 $ clusters
    left        = adjMat H.?? (H.Pos (H.idxs leftIdxs), H.Pos (H.idxs leftIdxs))
    right       =
        adjMat H.?? (H.Pos (H.idxs rightIdxs), H.Pos (H.idxs rightIdxs))
