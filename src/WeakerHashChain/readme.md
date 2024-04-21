WeakerHashChain
===============

WeakerHashChain is a variation on the HashChain search algorithm.
It remembers the rightmost position it has previously identified a
possible factor, and will not rescan past that in a new filtering phase.
This is a weaker form of recognition than scanning an entire chain, 
but means we don't have to scan data we've already looked at.

The technique was taken from 
*"Linear and Efficient String Matching Algorithms Based on Weak Factor Recognition"*
by DOMENICO CANTONE and SIMONE FARO, University of Catania
ARIANNA PAVONE, University of Messina.  

The idea to try using just one of the two techniques in that paper is mine.  
This algorithm is still quadratic in the worst case and is not linear.
However, experimentally it seems to be reasonably faster than HashChain,
except on long, low entropy data (e.g. genome).
