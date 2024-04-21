WeakerHashChain
===============

WeakerHashChain is a variation on the HashChain search algorithm.
It remembers the rightmost position it has previously identified a
possible factor, and will not rescan past that in a new filtering phase.
This is a weaker form of recognition than scanning an entire chain, 
but means we don't have to scan data we've already looked at.

Experimentally it seems to be reasonably faster than HashChain,
except on long, low entropy data (e.g. genome).
