{- Math.Clustering.Hierarchical.Spectral
Gregory W. Schwartz

Collects the functions pertaining to hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral
    ( hierarchicalSpectralCluster
    , getClusterItems
    ) where

-- Remote
import Data.Bool (bool)
import Data.Clustering.Hierarchical (Dendrogram (..))
import Math.Clustering.Spectral (spectralClusterNorm)
import Math.Modularity (getModularity, Q (..))
import qualified Data.Foldable as F
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H

-- Local

type AdjacencyMatrix = H.Matrix Double
type Items a         = V.Vector a

-- | Generates a tree through divisive hierarchical clustering using
-- Newman-Girvan modularity as a stopping criteria.
hierarchicalSpectralCluster :: Items a
                            -> AdjacencyMatrix
                            -> Dendrogram (V.Vector a, AdjacencyMatrix)
hierarchicalSpectralCluster !items !adjMat =
    if ngMod > 0 && H.rows adjMat > 1
        then Branch
                ngMod
                (hierarchicalSpectralCluster (getItems leftIdxs) left)
                (hierarchicalSpectralCluster (getItems rightIdxs) right)
        else 
            Leaf (items, adjMat)
  where
    clusters    = spectralClusterNorm adjMat
    ngMod       = unQ . getModularity clusters $ adjMat
    getIdxs val = VS.ifoldl' (\ !acc !i -> bool acc (i:acc) . (== val)) []
    leftIdxs    = getIdxs 0 $ clusters
    rightIdxs   = getIdxs 1 $ clusters
    left        = adjMat H.?? (H.Pos (H.idxs leftIdxs), H.Pos (H.idxs leftIdxs))
    right       =
        adjMat H.?? (H.Pos (H.idxs rightIdxs), H.Pos (H.idxs rightIdxs))
    getItems    =
        V.fromList . F.foldl' (\ !acc !i -> (items V.! i) : acc) [] . V.fromList
    
-- | Gather clusters (leaves) from the dendrogram.
getClusterItems :: Dendrogram (V.Vector a, AdjacencyMatrix) -> [V.Vector a]
getClusterItems = fmap fst . F.toList
