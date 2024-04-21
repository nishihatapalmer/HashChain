LinearHashChain
===============

LinearHashChain is a variant of the HashChain search algorithm 
that remains linear in the worst case, instead of quadratic as for HashChain.

It achieves linearity by ensuring we don't re-scan bytes in both the filtering 
and verification phases.

In the filtering phase, we remember the rightmost position previously 
identified as a possible factor, and don't re-scan past that in a 
further filtering scan.  We assume that the previous matches will
still hold up, which is a weaker form of recognition than scanning an
entire chain from end to start.  This is the same concept used in the
WeakerHashChain algorithm.

In the verification phase, we use a linear forward matching algorithm (KMP)
to identify matches which will not rescan previously verified positions.

The combination of both techniques ensures that performance remains linear
in the worst case, but is very fast and sublinear on average.
