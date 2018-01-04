{- Math.Clustering.Hierarchical.Sparse.Spectral
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering for sparse
data.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Sparse
    ( hierarchicalSpectralCluster
    , getClusterItems
    , test
    , test2
    , items2
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Tree (Tree (..))
import Math.Clustering.Spectral.Sparse (B (..), getB, spectralCluster)
import Math.Modularity.Sparse (getBModularity)
import Math.Modularity.Types (Q (..))
import qualified Data.Discrimination as D
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
        getSortedIdxs val =
            D.sort
                . VS.ifoldl' (\ !acc !i -> bool acc (i:acc) . (== val)) []
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
        extractRows mat = S.fromRowsL . fmap (S.extractRow mat)
        getItems    =
            V.fromList . F.foldl' (\ !acc !i -> (items V.! i) : acc) []

-- | Gather clusters (leaves) from the tree.
getClusterItems :: Foldable t => t (Items a, b) -> [V.Vector a]
getClusterItems = fmap fst . F.toList

test :: IO B
test = fmap B . S.chol . S.fromListDenseSM 7 $ ([1,0.423,0.56,0.004,0.004,0.24,0.013,0.423,1,0.466,0.00,50.004,0.013,0.043,0.56,0.466,1,0.004,0.003,0.01,0.035,0.004,0.005,0.004,1,0.71,0.008,0.009,0.004,0.004,0.003,0.71,1,0.008,0.008,0.24,0.013,0.01,0.008,0.008,1,0.513,0.013,0.043,0.035,0.009,0.008,0.513,1] :: [Double])

test2 = B $ S.fromListDenseSM 4 ([1,0, 1,0, 0,1, 0,1] :: [Double])
items2 = V.fromList ["a", "a", "b", "b"]
