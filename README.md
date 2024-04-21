HashChain
=========

HashChain is a family of very fast factor-based sublinear exact-matching search algorithms.  

They work by building a bloom-filter based on hashes of q-grams within the
pattern to be searched, and hashes of their adjacent q-gram.  This permits
the algorithm to efficiently reject non-adjacent q-grams in the text with
very high probability, allowing it to skip ahead of the mis-matching factor.

### Algorithms ###
Algorithms in the HashChain family include:

* HashChain - the original HashChain algorithm.
* SentinelHashChain - a faster HashChain using a search text modification hack.
* WeakerHashChain - a faster HashChain algorithm which does not re-scan data during filtering.
* LinearHashChain - HashChain with a guaranteed linear worst-case, based on Linear WFR.

### Similar algorithms ###

Similar algorithms include:

* Weak Factor Recognition (WFR) by Simone Faro, 
* Linear WFR (LWFR) by Domenico Cantone, Simone Faro and Arianna Pavone,
* Q-gram filtering (QF) by Branislav Ďurian, Hannu Peltola, Leena Salmela,
and Jorma Tarhio.

WFR uses a rolling hash function that gives it a longer pre-processing stage 
and which sets more bits in its filter, giving it a higher false positive 
rate.  It is especially effective on low entropy data as the rolling hash
expands the accessible space in the filter table.

LWFR uses two techniques to achieve linearity in the worst case, while
remaining sublinear on average.  First that it does not re-scan data it
has previously scanned during the filtering phrase, and
second, that it uses a forward linear matching algorithm (KMP) for
verification that likewise guarantees linearity.

QF uses an alignment concept; that successive q-grams of size Q 
read back in the text should be aligned in the pattern in the same way.  Either
they belong to the chain of q-grams that start from the end of the
pattern, or the second from the end, and so on up to Q chains. It appears
to be an efficient filter, but the space it has to set bits is
limited by Q rather than the size of a word in memory.


