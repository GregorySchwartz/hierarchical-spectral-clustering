{- Math.Clustering.Hierarchical.Types
Gregory W. Schwartz

Collects the types used in hierarchical clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Types
    ( ClusteringTree (..)
    , ClusteringVertex (..)
    , clusteringTreeToDendrogram
    ) where

-- Remote
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Monoid ((<>))
import Data.Tree (Tree (..))
import Math.Modularity.Types (Q (..))
import qualified Data.Vector as V

-- Local

type Items a         = V.Vector a
type ClusteringTree a b = Tree (ClusteringVertex a b)

data ClusteringVertex a b = ClusteringVertex
    { _clusteringItems  :: Items a
    , _clusteringMatrix :: b
    , _ngMod            :: Q
    } deriving (Eq, Ord, Read, Show)

-- | Convert a ClusteringTree to a Dendrogram.
clusteringTreeToDendrogram :: ClusteringTree a b -> Dendrogram (Items a, b)
clusteringTreeToDendrogram (Node { rootLabel = !n, subForest = []}) =
    Leaf (_clusteringItems n, _clusteringMatrix n)
clusteringTreeToDendrogram (Node { rootLabel = !n, subForest = [x, y]}) =
    Branch
        (unQ . _ngMod $ n)
        (clusteringTreeToDendrogram x)
        (clusteringTreeToDendrogram y)
clusteringTreeToDendrogram (Node { subForest = xs}) =
    error $ "Clustering tree has "
         <> (show $ length xs)
         <> " children. Requires two or none."
