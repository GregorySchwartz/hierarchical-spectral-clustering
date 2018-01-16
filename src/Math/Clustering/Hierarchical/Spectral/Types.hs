{- Math.Clustering.Hierarchical.Types
Gregory W. Schwartz

Collects the types used in hierarchical clustering.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module Math.Clustering.Hierarchical.Spectral.Types
    ( ClusteringTree (..)
    , ClusteringVertex (..)
    , clusteringTreeToDendrogram
    , getClusterItemsDend
    , getClusterItemsTree
    , Q (..)
    ) where

-- Remote
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Monoid ((<>))
import Data.Tree (Tree (..))
import GHC.Generics (Generic)
import Math.Modularity.Types (Q (..))
import Math.TreeFun.Tree (leaves)
import qualified Data.Foldable as F
import qualified Data.Vector as V

-- Local

type Items a         = V.Vector a
type ClusteringTree a = Tree (ClusteringVertex a)

data ClusteringVertex a = ClusteringVertex
    { _clusteringItems  :: !(Items a)
    , _ngMod            :: !Q
    } deriving (Eq, Ord, Read, Show, Generic)

-- | Convert a ClusteringTree to a Dendrogram. Modularity is the distance.
clusteringTreeToDendrogram :: ClusteringTree a -> Dendrogram (Items a)
clusteringTreeToDendrogram = fst . go
  where
    go (Node { rootLabel = !n, subForest = []}) =
        (Leaf (_clusteringItems n), 0)
    go (Node { rootLabel = !n, subForest = [x, y]}) =
        (Branch newD l r, newD)
      where
        newD = (unQ . _ngMod $ n) + max lDist rDist
        (!l, !lDist) = go x
        (!r, !rDist) = go y
    go (Node { subForest = xs}) =
        error $ "Clustering tree has "
            <> (show $ length xs)
            <> " children. Requires two or none."

-- | Gather clusters (leaves) from the dendrogram.
getClusterItemsDend :: Foldable t => t (Items a) -> [Items a]
getClusterItemsDend = F.toList

-- | Gather clusters (leaves) from the tree.
getClusterItemsTree :: ClusteringTree a -> [Items a]
getClusterItemsTree = fmap _clusteringItems . leaves
