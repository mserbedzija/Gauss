module Gauss.Utils (
    integralSqrt,
    cdiv
) where

integralSqrt :: (Integral a) => a -> Double
integralSqrt = sqrt . fromIntegral

cdiv :: (Integral a) => a -> a -> a
cdiv a b = let (q, r) = quotRem a b in  q + (if r /= 0 then 1 else 0)