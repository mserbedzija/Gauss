module Gauss.Prime
    (
      -- * General number-theoretic functions
      coprime
    , units
    , extractExponent
    , maxExponent
    
    -- * Arithmetic functions
    , totient
    
    -- * Sieve algorithms
    , sieve
    , sieveUsing
    , eratosthenes
    , filterMultiples
    
    -- * Coprime embeddings
    , CoprimeEmbedding(..)
    , identityEmbedding
    , embeddingMod2
    , embeddingMod6
    ) where

import Control.Monad       ( forM_, when )
import Control.Monad.ST    ( ST )
import Data.Array.ST       ( Ix, STUArray, runSTUArray, newArray, readArray, writeArray )
import Data.Array.Unboxed  ( assocs )
import Data.List           ( nub, sort )

import Gauss.Utils     ( integralSqrt, cdiv )

-- | Evaluates Euler's totient function for n.
totient :: (Integral a) => a -> a
totient n = fromIntegral . length $ filter (coprime n) [1 .. n - 1]

-- | Evaluates the exponentiation n^e in (Z/Zm)* using binary exponentiation.
modpow :: (Integral a, Num a) => a -> a -> a -> a
modpow n e m
    | m < 0     = error "modular base is undefined."
    | e < 0     = case inverseMod n m of
        Just n' -> modpow n' (-e) m
        Nothing -> error "inverse power does not exist."
    | e == 0    = 1
    | e == 1    = n `mod` m
    | odd e     = (n*r) `mod` m
    | otherwise = r
  where
    r = modpow ((n*n) `mod` m) (e `div` 2) m

-- | Evaluates the multiplicative inverse of n in Z/Zp for p a prime.
inverseModP :: (Integral a) => a -> a -> a
inverseModP n p = modpow n (p - 2) p

-- | Evaluates the multiplicative inverse of a in Z/Zn. If a and n
--   are not coprime, it returns Nothing.
inverseMod :: (Integral a) => a -> a -> Maybe a
inverseMod a n | coprime a n = Just (modpow a ((totient n) - 1) n)
               | otherwise   = Nothing

-- | A list containing the divisors of n; note that n is required to
--   be non-negative here.
factors :: (Integral a) => a -> [a]
factors 0 = []
factors 1 = [1]
factors n
    | n < 0     = error "`factors` expects a non-negative input."
    | otherwise = sort . ([1, n] ++) $ factors' 2
      where limit = floor . sqrt . fromIntegral $ n
            factors' d
                | d > limit      = []
                | n `mod` d == 0 = nub [d, div n d] ++ factors' (d + 1)
                | otherwise      = factors' (d + 1)

-- | A list of each distinct prime that divides n
primeFactors :: (Integral a) => a -> [a]
primeFactors n | n == 0 || n == 1 = []
               | otherwise        = filterMultiples (drop 1 $ factors n)

-- | The prime factorization of n; a list of pairs (p_i, e_i) satisfying
--   p_1^e_1 * ... * p_k^e_k = n.
primeFactorization :: (Integral a) => a -> [(a, a)]
primeFactorization n = map (\p -> (p, maxExponent p n)) (primeFactors n)

-- | Finds integers v and m satisfying (d^v)*m = n.
extractExponent :: (Integral a) => a -> a -> (a, a)
extractExponent d n
    | n `rem` d == 0 =
        let (v, m) = extractExponent d (n `quot` d)
        in  (v + 1, m)
    | otherwise      = (0, n)

-- | The maximum integer e >= 0 such that d^e divides n.
maxExponent :: (Integral a) => a -> a -> a
maxExponent d = fst . extractExponent d

-- | Checks for primality of n.
isPrime :: (Integral a) => a -> Bool
isPrime n 
    | n <= 1    = False
    | otherwise = isPrime' 2
      where
        limit = floor . sqrt . fromIntegral $ n
        isPrime' d | d > limit      = True
                   | n `mod` d == 0 = False
                   | otherwise      = isPrime' (d + 1)

-- | Determines coprimality between a and b
coprime :: (Integral a) => a -> a -> Bool
coprime a b = gcd a b == 1

-- | The unit (i.e., multiplicatively invertible) elements of Z/Zn.
units :: (Integral a) => a -> [a]
units 2 = [1]
units n = let units' = 1 : filter (coprime n) [2 .. (n + 1) `div` 2]
          in  (units' ++) . reverse . map ((-) n) $ units'

-- | Filters a list of integers such that a mod b != 0 for all
--   a b in the new list. The list is assumed to be given pre-sorted
--   in ascending order.
filterMultiples :: (Integral a) => [a] -> [a]
filterMultiples []     = []
filterMultiples (n:xs) = n : filterMultiples xs'
  where
    xs' = filter ((0 /= ) . (flip mod) n) xs

-- SIEVE ALGORITHMS

sieve :: (Integral a) => a -> [a]
sieve = sieveUsing embeddingMod6

sieveUsing :: (Integral a) => CoprimeEmbedding -> a -> [a]
sieveUsing f n
    | hasEnough = filter (<= n) base_
    | otherwise =
        let array = assocs $ runSTUArray $ eratosthenes f (fromIntegral n)
        in  base_ ++ [fromIntegral $ f <== p | (p, True) <- array]
  where
    base_     = map fromIntegral $ base f
    hasEnough = length base_ > 0 && n < last base_

-- | The sieve of Eratosthenes. It requires a coprime embedding to determine
-- the optimisation level.
eratosthenes :: CoprimeEmbedding -> Int -> ST s (STUArray s Int Bool)
eratosthenes f n = do
    let baseval_ = product $ base f
    -- Generate an array storing the encoded integers.
    sieve <- newArray (1, f ==> n) True
    forM_ [1 .. f ==> (floor . integralSqrt $ n)] $ \i -> do
        -- If the current integer is prime
        isPrime <- readArray sieve i
        when isPrime $ do
            -- Fetch the decoded prime
            let p = f <== i
            forM_ (branches f) $ \branch -> do
                -- Fetch the arithemtic branch's starting point
                let s = branch p
                -- Mark each multiple of p as non-prime.
                forM_ [f ==> s, f ==> (s + p*baseval_) .. f ==> n] $ \cmp -> do
                    writeArray sieve cmp False
    return sieve

-- | Given a prime p, let M_p denote the set of all numbers that are coprime
-- to p!. A bijection M_p -> N that comes equipped with maps {f_i : N -> M_p}
-- such that
--      1) {f_i(n) | n in N} is an arithmetic progression in M_p, and
--      2) The union of {f_i(N)} covers M_p,
-- will be called a `coprime embedding`.
data CoprimeEmbedding = CoprimeEmbedding
    { base     :: [Int]          -- ^ Base of the embedding.
    , (==>)    ::  Int -> Int    -- ^ Image function.
    , (<==)    ::  Int -> Int    -- ^ Inverse function
    , branches :: [Int -> Int]   -- ^ Arithmetic branches of the inverse.
    }

-- | The identity embedding N -> N; this runs a standard sieve of Eratosthenes
-- without any substantial optimisations. Note that one @embeddingMod2@ and
-- @embeddingMod6@ are faster, with @embeddingMod2@ being about twice as fast,
-- and @embeddingMod6@ being about 4 times as fast.
identityEmbedding :: CoprimeEmbedding
identityEmbedding = CoprimeEmbedding
    { base     = []
    , (==>)    = (subtract 1)
    , (<==)    = (+ 1)
    , branches = [\p -> p*p]
    }

-- | The embedding M_2 -> N; this is faster than @identityEmbedding@, although
-- slower than @embeddingMod6@.
embeddingMod2 :: CoprimeEmbedding
embeddingMod2 = CoprimeEmbedding
    { base     = [2]
    , (==>)    = (`div` 2) . (subtract 1)
    , (<==)    = (+ 1) . (* 2)
    , branches = [\p -> p*p]
    }

-- | The embedding M_3 -> N; this is faster than both @identityEmbedding@
-- and @embeddingMod2@.
embeddingMod6 :: CoprimeEmbedding
embeddingMod6 = CoprimeEmbedding
    { base     = [2, 3]
    , (==>)    = \n -> ((n - 1) `div` 6) + ((n + 1) `div` 6)
    , (<==)    = \n -> 3*n + 1 + (n `rem` 2)
    , branches =
        [ \p -> p*p
        , \p -> case p `rem` 6 of
                    1 -> 6*j*p + 5*p where j = (((p - 5) `div` 2) + 2) `div` 3
                    _ -> 6*j*p + p   where j = (((p - 1) `div` 2) + 2) `div` 3
        ]
    }