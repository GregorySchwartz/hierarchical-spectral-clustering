{- Math.Clustering.Hierarchical.Spectral.Load
Gregory W. Schwartz

Collects the functions pertaining to loading a matrix.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Math.Clustering.Hierarchical.Spectral.Load
    ( readDenseAdjMatrix
    , readSparseAdjMatrix
    , readEigenSparseAdjMatrix
    ) where

-- Remote
import Control.Monad.Except (runExceptT, ExceptT (..))
import Control.Monad.Managed (with, liftIO, Managed (..))
import Data.Maybe (fromMaybe, catMaybes)
import System.IO (Handle (..))
import qualified Data.ByteString.Streaming.Char8 as BS
import qualified Data.Csv as CSV
import qualified Data.Eigen.SparseMatrix as E
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.Sparse.Common as SH
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Numeric.LinearAlgebra as H
import qualified Streaming as S
import qualified Streaming.Cassava as S
import qualified Streaming.Prelude as S
import qualified Streaming.With.Lifted as SW

-- Local
import Math.Clustering.Hierarchical.Spectral.Types

-- | Generic error message.
errorMsg = error "Not correct format (requires row,column,value)"

-- | Parse a row of a label index file.
parseRow :: (T.Text, T.Text, Double) -> ((T.Text, T.Text), Double)
parseRow (i, j, v) = ((i, j), v)

-- | Ignore the disconnected vertices, not used (rather use very small weight).
ignoreDisconnected :: V.Vector T.Text
                   -> H.Matrix Double
                   -> (V.Vector T.Text, H.Matrix Double)
ignoreDisconnected items mat = (newItems, newMat)
  where
    newItems = V.fromList $ fmap ((V.!) items) valid
    newMat = mat H.?? (H.Pos $ H.idxs valid, H.Pos $ H.idxs valid)
    valid = catMaybes
          . zipWith (\x xs -> if sum xs > 0 then Just x else Nothing) [0..]
          . H.toLists
          $ mat

-- | Ensure symmetry.
symmetric :: [((Int, Int), Double)] -> [((Int, Int), Double)]
symmetric = concatMap (\((!i, !j), v) -> [((i, j), v), ((j, i), v)])

-- | Ensure zeros on diagonal.
zeroDiag :: [((Int, Int), Double)] -> [((Int, Int), Double)]
zeroDiag = filter (\((!i, !j), _) -> i /= j)

-- | Get the translated matrix indices.
getNewIndices
    -- :: (Eq a, Ord a)
    -- => [((a, a), Double)] -> [((Int, Int), Double)]
    :: [((T.Text, T.Text), Double)] -> [((Int, Int), Double)]
getNewIndices xs =
    fmap
        (\((!i,!j),!v) ->
              ( ( Map.findWithDefault eMsg i idxMap
                , Map.findWithDefault eMsg j idxMap
                )
              , v
              )
        )
        xs
  where
    eMsg     = error "Index not found during index conversion."
    indices  = getAllIndices xs
    idxMap   = Map.fromList $ zip indices [0 ..]

-- | Get the list of all indices.
getAllIndices :: (Eq a, Ord a) => [((a, a), Double)] -> [a]
getAllIndices xs = Set.toAscList . Set.union (getSet fst) $ getSet snd
  where
    getSet f = Set.fromList . fmap (f . fst) $ xs

-- | Get a dense adjacency matrix from a handle.
readDenseAdjMatrix :: CSV.DecodeOptions
                   -> Handle
                   -> IO (V.Vector T.Text, H.Matrix Double)
readDenseAdjMatrix decodeOpt handle = flip with return $ do
    let getAssocList = S.toList_ . S.map parseRow

    assocList <-
        fmap (either (error . show) id)
            . runExceptT
            . getAssocList
            . S.decodeWith decodeOpt S.NoHeader
            $ (BS.hGetContents handle :: BS.ByteString (ExceptT S.CsvParseException Managed) ())

    let items = V.fromList $ getAllIndices assocList
        mat   = H.assoc (V.length items, V.length items) 0
              . Set.toList
              . Set.fromList -- Ensure no duplicates.
              . symmetric -- Ensure symmetry.
              . zeroDiag -- Ensure zeros on diagonal.
              . getNewIndices -- Only look at present rows by converting indices.
              $ assocList

    return (items, mat)

-- | Get a sparse adjacency matrix from a handle.
readSparseAdjMatrix :: CSV.DecodeOptions
                    -> Handle
                    -> IO (V.Vector T.Text, SH.SpMatrix Double)
readSparseAdjMatrix decodeOpt handle = flip with return $ do
    let getAssocList = S.toList_ . S.map parseRow

    assocList <-
        fmap (either (error . show) id)
            . runExceptT
            . getAssocList
            . S.decodeWith decodeOpt S.NoHeader
            $ (BS.hGetContents handle :: BS.ByteString (ExceptT S.CsvParseException Managed) ())

    let items = V.fromList $ getAllIndices assocList
        mat   = SH.fromListSM (V.length items, V.length items)
              . Set.toList
              . Set.fromList -- Ensure no duplicates.
              . fmap (\((i, j), v) -> (i, j, v))
              . symmetric -- Ensure symmetry.
              . zeroDiag -- Ensure zeros on diagonal.
              . getNewIndices -- Only look at present rows by converting indices.
              $ assocList

    return (items, mat)

-- | Get a sparse adjacency matrix from a handle.
readEigenSparseAdjMatrix :: CSV.DecodeOptions
                    -> Handle
                    -> IO (V.Vector T.Text, E.SparseMatrixXd)
readEigenSparseAdjMatrix decodeOpt handle = flip with return $ do
    let getAssocList = S.toList_ . S.map parseRow

    assocList <-
        fmap (either (error . show) id)
            . runExceptT
            . getAssocList
            . S.decodeWith decodeOpt S.NoHeader
            $ (BS.hGetContents handle :: BS.ByteString (ExceptT S.CsvParseException Managed) ())

    let items = V.fromList $ getAllIndices assocList
        mat   = E.fromList (V.length items) (V.length items)
              . Set.toList
              . Set.fromList -- Ensure no duplicates.
              . fmap (\((i, j), v) -> (i, j, v))
              . symmetric -- Ensure symmetry.
              . zeroDiag -- Ensure zeros on diagonal.
              . getNewIndices -- Only look at present rows by converting indices.
              $ assocList

    return (items, mat)
