{- Math.Clustering.Hierarchical.Spectral.Utility
Gregory W. Schwartz

Collects utility functions for the clustering section of the program.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Utility
  ( subsetVector
  , permutationTest
  , permutationTestSparse
  ) where

-- Remote
import Data.List (genericLength)
import Data.Maybe (fromMaybe)
import qualified Data.Foldable as F
import qualified Data.Sparse.Common as S
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G
import qualified System.Random.MWC as R
import qualified System.Random.MWC.Distributions as R

-- Local


subsetVector :: V.Vector a -> [Int] -> V.Vector a
subsetVector xs =
    V.fromList
        . F.foldr' (\ !i !acc
                    -> ( fromMaybe (error "Out of bounds in subsetVector (do the number of observations and features match the size of the matrix? Is the input format correct?).")
                      $ xs V.!? i
                      ) : acc
                    ) []


-- | Permutation test returning p-value.
permutationTest :: (G.Vector v a)
                => (v a -> b -> Double)
                -> Int
                -> v a
                -> b
                -> IO Double
permutationTest f runs labs samples = do
  let obs = f labs samples
      run r = do
        g <- R.initialize . V.singleton $ fromIntegral r
        ls <- R.uniformShuffle labs g
        return $ f ls samples

  n <- genericLength . filter (>= obs) <$> mapM run [1..runs]

  return $ n / fromIntegral runs

-- | Permutation test returning p-value from sparse vector.
permutationTestSparse :: (S.SpVector Double -> b -> Double)
                -> Int
                -> S.SpVector Double
                -> b
                -> IO Double
permutationTestSparse f runs labsInit samples = do
  let obs = f labsInit samples
      labs = U.fromList $ S.toDenseListSV labsInit
      run r = do
        g <- R.initialize . V.singleton $ fromIntegral r
        ls <- S.vr . U.toList <$> R.uniformShuffle labs g
        return $ f ls samples

  n <- genericLength . filter (>= obs) <$> mapM run [1..runs]

  return $ n / fromIntegral runs
