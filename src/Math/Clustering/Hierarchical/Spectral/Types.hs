{- Math.Clustering.Hierarchical.Types
Gregory W. Schwartz

Collects the types used in hierarchical clustering.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE StandaloneDeriving #-}

module Math.Clustering.Hierarchical.Spectral.Types
    ( ClusteringTree (..)
    , ClusteringVertex (..)
    , GenericClusteringTree (..)
    , GenericClusteringVertex (..)
    , EigenGroup (..)
    , clusteringTreeToDendrogram
    , clusteringTreeToDendrogramCumulative
    , clusteringTreeToGenericClusteringTree
    , getClusterItemsDend
    , getClusterItemsTree
    , getClusterItemsGenericTree
    , Q (..)
    , NumEigen (..)
    ) where

-- Remote
import Data.Clustering.Hierarchical (Dendrogram (..))
import Data.Maybe (catMaybes)
import Data.Monoid ((<>))
import Data.Tree (Tree (..))
import GHC.Generics (Generic)
import Math.Modularity.Types (Q (..))
import Math.TreeFun.Tree (leaves)
import qualified Data.Aeson as A
import qualified Data.Foldable as F
import qualified Data.Vector as V

-- Local

type Items a         = V.Vector a
type ClusteringTree a = Tree (ClusteringVertex a)
type GenericClusteringTree a = Tree (GenericClusteringVertex a)
type NumEigen = Int

data EigenGroup = SignGroup | KMeansGroup deriving (Read, Show, Generic)

data ClusteringVertex a = ClusteringVertex
    { _clusteringItems  :: !(Items a)
    , _ngMod            :: !Q
    } deriving (Eq, Ord, Read, Show, Generic)

data GenericClusteringVertex a = GenericClusteringVertex
    { _item :: !(Maybe (Items a))
    , _distance :: !(Maybe Double)
    } deriving (Eq, Ord, Read, Show, Generic)

-- | Convert a ClusteringTree to a Dendrogram. Modularity is the distance.
clusteringTreeToDendrogram :: ClusteringTree a -> Dendrogram (Items a)
clusteringTreeToDendrogram = go
  where
    go (Node { rootLabel = !n, subForest = []}) = Leaf $ _clusteringItems n
    go (Node { rootLabel = !n, subForest = [x, y]}) =
        Branch (unQ . _ngMod $ n) (go x) (go y)
    go (Node { subForest = xs}) =
        error $ "Clustering tree has "
            <> (show $ length xs)
            <> " children. Requires two or none."

-- | Convert a ClusteringTree to a Dendrogram. Modularity is the distance, such
-- that the distance is the modularity plus the maximum distance of each branch.
clusteringTreeToDendrogramCumulative :: ClusteringTree a -> Dendrogram (Items a)
clusteringTreeToDendrogramCumulative = fst . go
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

-- | Convert a ClusteringTree to a GenericClusteringVertex tree (more
-- standardized for our purposes).
clusteringTreeToGenericClusteringTree :: ClusteringTree a
                                      -> GenericClusteringTree a
clusteringTreeToGenericClusteringTree = go
  where
    go (Node { rootLabel = !n, subForest = []}) =
      Node { rootLabel = ( GenericClusteringVertex
                             { _item = Just $ _clusteringItems n
                             , _distance = Nothing
                             }
                         )
           , subForest = []
           }
    go (Node { rootLabel = !n, subForest = xs }) =
      Node { rootLabel = ( GenericClusteringVertex
                            { _item = Nothing
                            , _distance = Just . unQ . _ngMod $ n
                            }
                         )
           , subForest = fmap go xs
           }

-- | Gather clusters (leaves) from the dendrogram.
getClusterItemsDend :: Foldable t => t (Items a) -> [Items a]
getClusterItemsDend = F.toList

-- | Gather clusters (leaves) from the tree.
getClusterItemsTree :: ClusteringTree a -> [Items a]
getClusterItemsTree = fmap _clusteringItems . leaves
                    
-- | Gather clusters (leaves) from the generic tree.
getClusterItemsGenericTree :: GenericClusteringTree a -> [Items a]
getClusterItemsGenericTree = catMaybes . fmap _item . leaves

deriving instance (Read a) => Read (Dendrogram a)
deriving instance Generic (Dendrogram a)

instance (A.ToJSON a) => A.ToJSON (Dendrogram a) where
    toEncoding = A.genericToEncoding A.defaultOptions
instance (A.FromJSON a) => A.FromJSON (Dendrogram a)

instance (A.ToJSON a) => A.ToJSON (GenericClusteringVertex a) where
    toEncoding = A.genericToEncoding A.defaultOptions
instance (A.FromJSON a) => A.FromJSON (GenericClusteringVertex a)
