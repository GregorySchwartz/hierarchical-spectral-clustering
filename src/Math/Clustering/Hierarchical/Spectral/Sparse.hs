{- Math.Clustering.Hierarchical.Spectral.Sparse
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering for sparse
data.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Sparse
    ( hierarchicalSpectralCluster
    , FeatureMatrix (..)
    , Items (..)
    , ShowB (..)
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Sparse (B (..), getB, spectralCluster)
import Math.Modularity.Sparse (getBModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Foldable as F
import qualified Data.Sparse.Common as S
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra.Sparse as S

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

type FeatureMatrix   = S.SpMatrix Double
type Items a         = V.Vector a
type ShowB           = ((Int, Int), [(Int, Int, Double)])

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria. Assumes the feature matrix
-- has column features and row observations. Can use FeatureMatrix or a
-- pre-generated B matrix. See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
hierarchicalSpectralCluster :: Items a
                            -> Either FeatureMatrix B
                            -> ClusteringTree a ShowB
hierarchicalSpectralCluster !initItems = go initItems . either getB id
  where
    go :: Items a -> B -> ClusteringTree a ShowB
    go items b =
        if ngMod > Q 0 && (S.nrows $ unB b) > 1
            then
                Node { rootLabel = vertex
                     , subForest = [ go (getItems leftIdxs) left
                                   , go (getItems rightIdxs) right
                                   ]
                     }

            else
                Node {rootLabel = vertex, subForest = []}
      where
        vertex      = ClusteringVertex
                        { _clusteringItems = items
                        , _clusteringMatrix =
                            (S.dimSM . unB $ b, S.toListSM . unB $ b)
                        , _ngMod = ngMod
                        }
        clusters :: S.SpVector Double
        clusters    = spectralCluster b
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
            V.fromList . F.foldr' (\ !i !acc -> (items V.! i) : acc) []
