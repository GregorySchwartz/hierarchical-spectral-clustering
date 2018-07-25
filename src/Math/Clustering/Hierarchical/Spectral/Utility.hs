{- Math.Clustering.Hierarchical.Spectral.Utility
Gregory W. Schwartz

Collects utility functions for the clustering section of the program.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Utility
  ( subsetVector
  ) where

-- Remote
import Data.Maybe (fromMaybe)
import qualified Data.Foldable as F
import qualified Data.Vector as V

-- Local


subsetVector :: V.Vector a -> [Int] -> V.Vector a
subsetVector xs =
    V.fromList
        . F.foldr' (\ !i !acc
                    -> ( fromMaybe (error "Out of bounds in subsetVector.")
                      $ xs V.!? i
                      ) : acc
                    ) []
