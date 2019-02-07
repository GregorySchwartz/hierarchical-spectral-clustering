{- cluster-tree
Gregory W. Schwartz

Hierarchical spectral clustering of data.
-}

{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeOperators     #-}
{-# LANGUAGE StandaloneDeriving #-}

module Main where

-- Standard
import Data.Maybe (fromMaybe, catMaybes)
import GHC.Generics

-- Cabal
import Data.Char (ord)
import Data.List (intercalate)
import Data.Monoid ((<>))
import Options.Generic
import Safe (atMay)
import System.IO (stdin)
import Text.Read (readMaybe)
import TextShow (showt)
import qualified Control.Lens as L
import qualified Data.Aeson as A
import qualified Data.Aeson.Encode.Pretty as A
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Csv as CSV
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Numeric.LinearAlgebra as H
import qualified System.FilePath as File
import Math.Graph.Types

-- Local
import Math.Clustering.Hierarchical.Spectral.Load
import Math.Clustering.Hierarchical.Spectral.Types
import Math.Graph.Components
import qualified Math.Clustering.Hierarchical.Spectral.Dense as HD
import qualified Math.Clustering.Hierarchical.Spectral.Eigen.AdjacencyMatrix as HS

newtype Delimiter  = Delimiter { unDelimiter :: Char } deriving (Read, Show)
newtype Row        = Row { unRow :: Int } deriving (Eq, Ord, Read, Show)
newtype Column     = Column { unColumn :: Int } deriving (Eq, Ord, Read, Show)
newtype OutputTree = OutputTree { unOutputTree :: String } deriving (Read, Show)
newtype MinSize    = MinSize { unMinSize :: Int } deriving (Read, Show)
newtype NumEigen   = NumEigen { unNumEigen :: Int } deriving (Read, Show)

data ClusteringType = Sparse | Dense deriving (Read, Show)
data Components a = Single (ClusteringTree a) | Multiple [ClusteringTree a]

instance A.ToJSON Q where
      toEncoding = A.genericToEncoding A.defaultOptions
instance A.FromJSON Q

instance (A.ToJSON a) => A.ToJSON (ClusteringVertex a) where
      toEncoding = A.genericToEncoding A.defaultOptions
instance (A.FromJSON a) => A.FromJSON (ClusteringVertex a)

-- | Command line arguments
data Options = Options { clusteringType :: Maybe String
                                       <?> "([Sparse] | Dense) Method for clustering data."
                       , delimiter      :: Maybe Char
                                       <?> "([,] | CHAR) The delimiter of the CSV file. Format is row,column,value with no header."
                       , minSize        :: Maybe Int
                                       <?> "([Nothing] | INT) Minimum size of a cluster."
                       , numEigen       :: Maybe Int
                                       <?> "([1] | INT) Number of eigenvectors to use while clustering with kmeans. Takes from the first eigenvector. Recommended to start at 2 and work up from there if needed."
                       , minModularity  :: Maybe Double
                                       <?> "([0] | DOUBLE) Minimum modularity to be over to continue recursion."
                       , eigenGroup     :: Maybe String
                                       <?> "([SignGroup] | KMeansGroup) Whether to group the eigenvector using the sign or kmeans while clustering. While the default is sign, kmeans may be more accurate (but starting points are arbitrary)."
                       , separateComponents :: Bool
                                           <?> "Whether to first separate connected components of the graph first. Will output a dendrogram for each component with the name of the tree and the number of nodes within the tree, along with the base set by --output-tree."
                       , outputTree     :: Maybe String
                                       <?> "([Nothing] | FILE) The name of the file to output the tree in JSON format."
                       }
               deriving (Generic)

modifiers :: Modifiers
modifiers = lispCaseModifiers { shortNameModifier = short }
  where
    short "minSize" = Just 'S'
    short x         = firstLetter x

instance ParseRecord Options where
    parseRecord = parseRecordWithModifiers modifiers

main :: IO ()
main = do
    opts <- getRecord "cluster-tree, Gregory W. Schwartz.\
                      \ Hierarchical spectral clustering of data Computes real\
                      \ symmetric part of matrix, so ensure the input is real\
                      \ and symmetric. Diagonal should be 0s for\
                      \ adjacency matrix.\
                      \ Format is an edge list of node1,node2,value with no header, fed in through STDIN.\
                      \ Must end with a newline."

    let readOrErr err       = fromMaybe (error err) . readMaybe
        clusteringType'     =
          maybe Sparse (readOrErr "Cannot read --clustering-type")
            . unHelpful
            . clusteringType
            $ opts
        delim'              =
            Delimiter . fromMaybe ',' . unHelpful . delimiter $ opts
        minSize'            = fmap MinSize . unHelpful . minSize $ opts
        numEigen'           = fmap NumEigen . unHelpful . numEigen $ opts
        minModularity'      = fmap Q . unHelpful . minModularity $ opts
        eigenGroup'         =
          maybe SignGroup (readOrErr "Cannot read --eigen-group")
            . unHelpful
            . eigenGroup
            $ opts
        separateComponents' = unHelpful . separateComponents $ opts
        outputTree'         = fmap OutputTree . unHelpful . outputTree $ opts
        decodeOpt       = CSV.defaultDecodeOptions
                            { CSV.decDelimiter =
                                fromIntegral (ord . unDelimiter $ delim')
                            }
        encodeOpt       = CSV.defaultEncodeOptions
                            { CSV.encDelimiter =
                                fromIntegral (ord . unDelimiter $ delim')
                            }

    clusteringTree <-
        case clusteringType' of
            Dense -> do
                (items, mat) <- readDenseAdjMatrix decodeOpt stdin

                let cluster items = HD.hierarchicalSpectralClusterAdj
                                      eigenGroup'
                                      (fmap unNumEigen numEigen')
                                      (fmap unMinSize minSize')
                                      minModularity'
                                      items

                return $
                    if separateComponents'
                        then Multiple
                           . fmap (uncurry cluster)
                           . getComponentMatsItems items
                           $ mat
                        else Single $ cluster items mat
            Sparse -> do
                (items, mat) <- readEigenSparseAdjMatrix decodeOpt stdin

                let cluster items = HS.hierarchicalSpectralCluster
                                      eigenGroup'
                                      (fmap unNumEigen numEigen')
                                      (fmap unMinSize minSize')
                                      minModularity'
                                      items

                return $
                    if separateComponents'
                        then Multiple
                           . fmap (uncurry cluster)
                           . getComponentMatsItems items
                           $ mat
                        else Single $ cluster items mat

    body <- case clusteringTree of
        (Single ct) -> do
            let clustering = zip ([1..] :: [Int]) . getClusterItemsTree $ ct
                body :: [(T.Text, T.Text)]
                body = concatMap
                        (\(!c, xs) -> fmap (\ !x -> (x, showt c)) . V.toList $ xs)
                        clustering

            case outputTree' of
                Nothing                -> return ()
                Just (OutputTree file) -> B.writeFile file
                                        . A.encodePretty
                                        . clusteringTreeToDendrogram
                                        $ ct

            return body
        (Multiple cts) -> do
            let clustering =
                    fmap (L.over L._2 (zip ([1..] :: [Int]) . getClusterItemsTree))
                        . zip ([1..] :: [Int])
                        $ cts
                body :: [(T.Text, T.Text)]
                body = concatMap
                        (\ (!t, xs)
                        -> concatMap (\ (!c, ys)
                                     -> fmap (\ x -> (x, showt t <> "/" <> showt c))
                                      . V.toList
                                      $ ys
                                     )
                           xs
                        )
                        clustering

            case outputTree' of
                Nothing                -> return ()
                Just (OutputTree file) -> do
                    let getSize :: ClusteringTree a -> Int
                        getSize = sum . fmap V.length . getClusterItemsTree
                        getFileName :: Int -> Int -> String
                        getFileName t n =
                            uncurry (File.</>)
                                . L.over L._2 (\x -> intercalate "_" ["tree", show t, "size", show n, x])
                                . File.splitFileName
                                $ file
                        write :: (A.ToJSON a) => Int -> ClusteringTree a -> IO ()
                        write t tree = do
                            B.writeFile (getFileName t (getSize tree))
                                . A.encodePretty
                                . clusteringTreeToDendrogram
                                $ tree

                    mapM_ (uncurry write) . zip [1..] $ cts

            return body

    -- | Print final result.
    B.putStr
        . (<>) "item,cluster\n"
        . CSV.encodeWith encodeOpt
        $ body

    return ()
