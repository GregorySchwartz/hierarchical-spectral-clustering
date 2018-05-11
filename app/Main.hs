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
import Data.Monoid ((<>))
import Options.Generic
import Safe (atMay)
import System.IO (stdin)
import qualified Data.Aeson as A
import qualified Data.Aeson.Encode.Pretty as A
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Csv as CSV
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Clustering.Hierarchical.Spectral.Dense
import Math.Clustering.Hierarchical.Spectral.Load
import Math.Clustering.Hierarchical.Spectral.Types


newtype Delimiter  = Delimiter { unDelimiter :: Char } deriving (Read, Show)
newtype Row        = Row { unRow :: Int } deriving (Eq, Ord, Read, Show)
newtype Column     = Column { unColumn :: Int } deriving (Eq, Ord, Read, Show)
newtype OutputTree = OutputTree { unOutputTree :: String } deriving (Read, Show)
newtype MinSize    = MinSize { unMinSize :: Int } deriving (Read, Show)

data ClusteringType = Dense deriving (Read, Show)

instance A.ToJSON Q where
      toEncoding = A.genericToEncoding A.defaultOptions
instance A.FromJSON Q

instance (A.ToJSON a) => A.ToJSON (ClusteringVertex a) where
      toEncoding = A.genericToEncoding A.defaultOptions
instance (A.FromJSON a) => A.FromJSON (ClusteringVertex a)

-- | Command line arguments
data Options = Options { clusteringType :: Maybe String
                                       <?> "([Dense]) Method for clustering data. Dense only so far."
                       , delimiter      :: Maybe Char
                                       <?> "([,] | CHAR) The delimiter of the CSV file. Format is row,column,value with no header."
                       , minSize        :: Maybe Int
                                       <?> "([Nothing] | INT) Minimum size of a cluster."
                       , outputTree     :: Maybe String
                                       <?> "([Nothing] | FILE) The name of the file to output the tree in JSON format."
                       }
               deriving (Generic)

instance ParseRecord Options

main :: IO ()
main = do
    opts <- getRecord "cluster-tree, Gregory W. Schwartz.\
                      \ Hierarchical spectral clustering of data Computes real\
                      \ symmetric part of matrix, so ensure the input is real\
                      \ and symmetric. Diagonal should be 0s for\
                      \ adjacency matrix.\
                      \ Format is row,column,value with no header."

    let clusteringType' = maybe Dense read . unHelpful . clusteringType $ opts
        delim'          =
            Delimiter . fromMaybe ',' . unHelpful . delimiter $ opts
        minSize'        = fmap MinSize . unHelpful . minSize $ opts
        outputTree'     = fmap OutputTree . unHelpful . outputTree $ opts
        decodeOpt       = CSV.defaultDecodeOptions
                            { CSV.decDelimiter =
                                fromIntegral (ord . unDelimiter $ delim')
                            }
        encodeOpt       = CSV.defaultEncodeOptions
                            { CSV.encDelimiter =
                                fromIntegral (ord . unDelimiter $ delim')
                            }

    clusteringTree <- do
        (items, mat) <- readDenseAdjMatrix decodeOpt stdin

        return $ hierarchicalSpectralCluster (fmap unMinSize minSize') items mat

    let clustering = zip [1..] . getClusterItemsTree $ clusteringTree
        body :: [(T.Text, Double)]
        body = concatMap
                (\(!c, xs) -> fmap (\ !x -> (x, c)) . V.toList $ xs)
                clustering

    case outputTree' of
        Nothing                -> return ()
        Just (OutputTree file) -> B.writeFile file
                                . A.encodePretty
                                . clusteringTreeToDendrogram
                                $ clusteringTree

    -- | Print final result.
    B.putStr
        . (<>) "item,cluster\n"
        . CSV.encodeWith encodeOpt
        $ body

    return ()
