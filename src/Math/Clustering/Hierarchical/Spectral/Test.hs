{- Math.Clustering.Hierarchical.Spectral.Test
Gregory W. Schwartz

Collects the functions pertaining to testing hierarchical spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Hierarchical.Spectral.Test where

-- Remote
import Data.List (tails)
import Data.Maybe (fromMaybe)
import Data.Monoid ((<>))
import Math.Clustering.Spectral.Sparse (getB, B (..))
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.Sparse.Common as S
import qualified Data.Eigen.SparseMatrix as E
import qualified Data.Vector as V
import qualified Math.Clustering.Spectral.Eigen.FeatureMatrix as EF
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Clustering.Hierarchical.Spectral.Types
import Math.Clustering.Hierarchical.Spectral.Sparse
import qualified Math.Clustering.Hierarchical.Spectral.Dense as Dense
import qualified Math.Clustering.Hierarchical.Spectral.Eigen.FeatureMatrix as EF
import qualified Math.Clustering.Hierarchical.Spectral.Eigen.AdjacencyMatrix as EA

newtype QGram = QGram { unQGram :: String } deriving (Eq, Ord, Read, Show)
newtype QGramMap = QGramMap
    { unQGramMap :: Map.Map QGram Int
    } deriving (Eq,Ord,Read,Show)

exampleData :: [String]
exampleData = [ "600 MOUNTAIN AVENUE"
              , "700 MOUNTAIN AVE"
              , "600-700 MOUNTAIN AVE"
              , "100 DIAMOND HILL RD"
              , "100 DIAMOND HILL ROAD"
              , "123 SPRINGFIELD AVENUE"
              , "123 SPRINFGIELD AVE"
              ]

exampleItems :: V.Vector String
exampleItems = V.fromList exampleData

-- | Add beginning and ending symbols.
addBorders :: String -> String
addBorders x = "##" <> x <> "$$"

-- | Generate qgrams for a string.
getQGrams :: Int -> String -> [QGram]
getQGrams n =
    fmap QGram . filter ((== n) . length) . fmap (take n) . tails . addBorders

-- | Get mapping of qgrams to indices.
getQGramMap :: [QGram] -> QGramMap
getQGramMap = QGramMap . Map.fromList . flip zip [0,1..]

-- | Convert a record to a vector.
recordToRow :: Int -> QGramMap -> String -> S.SpVector Double
recordToRow n (QGramMap qgramMap) = S.fromListSV (Map.size qgramMap)
                                  . Map.toAscList
                                  . Map.fromListWith (+)
                                  . flip zip [1,1..]
                                  . fromMaybe (error "Invalid qgram.")
                                  . mapM (flip Map.lookup qgramMap)
                                  . getQGrams n

-- | Generate the matrix of qgrams from a list of records and qgram length.
exampleMatrix :: Int -> [String] -> S.SpMatrix Double
exampleMatrix n records =
    S.transposeSM . S.fromColsL . fmap (recordToRow n qgramMap) $ records
  where
    qgramMap = getQGramMap
             . Set.toList
             . Set.fromList
             . concatMap (getQGrams n)
             $ records

clusterExample = hierarchicalSpectralCluster
                    SignGroup
                    True
                    Nothing
                    Nothing
                    Nothing
                    exampleItems
                    (Left $ exampleMatrix 3 exampleData)

clusterKExample = hierarchicalSpectralCluster
                    KMeansGroup
                    True
                    (Just 2)
                    Nothing
                    Nothing
                    exampleItems
                    (Left $ exampleMatrix 3 exampleData)

clusterAdjExample = hierarchicalSpectralClusterAdj
                        SignGroup
                        Nothing
                        Nothing
                        Nothing
                        exampleItems
                        adjacencyExample

clusterKAdjExample = hierarchicalSpectralClusterAdj
                        KMeansGroup
                        (Just 2)
                        Nothing
                        Nothing
                        exampleItems
                        adjacencyExample

adjacencyExample :: S.SpMatrix Double
adjacencyExample = S.filterSM (\i j _ -> i /= j) $ (unB b) S.##^ (unB b)
  where
    b = getB True $ exampleMatrix 3 exampleData

denseAdjacencyExample :: H.Matrix Double
denseAdjacencyExample = H.assoc (S.dimSM adjacencyExample) 0
                      . fmap (\(!x, !y, !z) -> if x == y then ((x, y), 0) else ((x, y), z))
                      . S.toListSM
                      $ adjacencyExample

denseClusterExample = Dense.hierarchicalSpectralCluster
                        SignGroup
                        Nothing
                        Nothing
                        Nothing
                        exampleItems
                        denseAdjacencyExample

denseClusterKExample = Dense.hierarchicalSpectralCluster
                        KMeansGroup
                        (Just 2)
                        Nothing
                        Nothing
                        exampleItems
                        denseAdjacencyExample

-- | Generate the matrix of qgrams from a list of records and qgram length.
exampleEigenMatrix :: Int -> [String] -> E.SparseMatrixXd
exampleEigenMatrix n records = E.fromList (S.nrows mat) (S.ncols mat) . S.toListSM $ mat
  where
    mat = exampleMatrix n records

adjacencyEigenExample :: E.SparseMatrixXd
adjacencyEigenExample = E._imap (\i j v -> if i == j then 0 else v)
                      $ (EF.unB b) * E.transpose (EF.unB b)
  where
    b = EF.getB True $ exampleEigenMatrix 3 exampleData

clusterEigenExample = EF.hierarchicalSpectralCluster
                        SignGroup
                        True
                        Nothing
                        Nothing
                        Nothing
                        exampleItems
                        (Left $ exampleEigenMatrix 3 exampleData)

clusterKEigenExample = EF.hierarchicalSpectralCluster
                        KMeansGroup
                        True
                        (Just 2)
                        Nothing
                        Nothing
                        exampleItems
                        (Left $ exampleEigenMatrix 3 exampleData)

clusterAdjEigenExample = EA.hierarchicalSpectralCluster
                            SignGroup
                            Nothing
                            Nothing
                            Nothing
                            exampleItems
                            adjacencyEigenExample

clusterKAdjEigenExample = EA.hierarchicalSpectralCluster
                        KMeansGroup
                        (Just 2)
                        Nothing
                        Nothing
                        exampleItems
                        adjacencyEigenExample
