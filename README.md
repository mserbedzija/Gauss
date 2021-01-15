# Gauss
A collection of number-theoretic algorithms in Haskell. It features an optimised version of the sieve of Eratosthenes that relies on _coprime embeddings_.

### Examples
```haskell
-- How many times does 2 divide 4024050954978504420349706240?
extractExponent 2 4024050954978504420349706240

-- Give me the prime factorization of 1051050
primeFactorization 1051050

-- What are the unit elements of Z/Z(210)?
units 210

-- Fetch all primes <= 10^7!
sieve (10^7)

-- Use the M_2 -> N embedding to find all primes <= 100!
sieveUsing embeddingMod2 100
```
