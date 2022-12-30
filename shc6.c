/*
 * Copyright 2022 Matt Palmer.  All rights reserved.
 *
 * This is an implementation of the simple HashChain algorithm, (currently unpublished) by Matt Palmer.
 * It is a factor search similar to WFR or the QF family of algorithms.
 *
 * It builds a hash table containing entries for chains of hashes.  Hashes are chained together by
 * placing a fingerprint of the *next* hash into the entry for the *current* hash.  This enables
 * a check for the second hash value to be performed without requiring a second lookup in the hash table.
 *
 * It only creates Q chains of hashes from the end of the pattern back to the start.
 * This makes pre-processing much faster than the general hash chain algorithm.
 * It performs badly on low alphabets as it does not use a rolling hash to expand the effective alphabet.
 * It also performs worse on longer patterns than the general algorithm.
 * It can perform very well on higher alphabet searches, for e.g. up to 2048 in length.
 *
 * This implementation is written to integrate with the SMART string search benchmarking tool,
 * by Simone Faro, Matt Palmer, Stefano Stefano Scafiti and Thierry Lecroq.
*/

#include "include/define.h"
#include "include/main.h"

/*
 * Hash table size.  Must be a power of two, minimum size is 32.
 */
#define ASIZE 4096

/*
 * Bit shift for each of the chain hash byte components.
 */
#define S3    2

/*
 * Number of bytes in a q-gram.
 * Chain hash functions defined below must be written to process this number of bytes.
 */
#define	Q     6

/*
 * Functions and calculated parameters.
 * Hash functions must be written to use the number of bytes defined in Q. They scan backwards from the initial position.
 */
#define HASH(x, p, s)     ((((((((((x[p] << (s)) + x[p - 1]) << (s)) + x[p - 2]) << (s)) + x[p - 3]) << (s)) + x[p - 4]) << (s)) + x[p - 5])
#define CHAIN_HASH(x, p)  HASH((x), (p), (S3))                     // Hash function for chain hashes, using the S3 bitshift.
#define FINGERPRINT(H)    (1U << ((H) & 0x1F))                     // Hash fingerprint, taking low 5 bits of the hash to set one of 32 bits.
#define TABLE_MASK        ((ASIZE) - 1)                            // Mask for table is one less than the power of two size.
#define Q2                (Q + Q)                                  // 2 Qs.
#define END_FIRST_QGRAM   (Q - 1)                                  // Position of the end of the first q-gram.
#define END_SECOND_QGRAM  (Q2 - 1)                                 // Position of the end of the second q-gram.

/*
 * Builds the hash table B of size ASIZE for a string x of length m.
 * Returns the 32-bit hash value of matching the entire pattern.
 */
unsigned int preprocessing(const unsigned char *x, int m, unsigned int *B) {

    // 0. Zero out the hash table.
    for (int i = 0; i < ASIZE; i++) B[i] = 0;

    // 1. Calculate all the chain hashes, ending with processing the entire pattern so H has the cumulative value.
    unsigned int H;
    for (int chain_no = Q; chain_no >= 1; chain_no--)
    {
        H = CHAIN_HASH(x, m - chain_no);
        for (int chain_pos = m - chain_no - Q; chain_pos >= END_FIRST_QGRAM; chain_pos -=Q)
        {
            unsigned int H_last = H;
            H = CHAIN_HASH(x, chain_pos);
            B[H_last & TABLE_MASK] |= FINGERPRINT(H);
        }
    }

    // 2. Add in hashes for the first qgrams that have no preceding value.  Only set a value if there is nothing there already.
    unsigned int F;
    int stop = MIN(m, END_SECOND_QGRAM);
    for (int chain_pos = END_FIRST_QGRAM; chain_pos < stop; chain_pos++)
    {
        F = CHAIN_HASH(x, chain_pos);
        if (!B[F & TABLE_MASK]) B[F & TABLE_MASK] = FINGERPRINT(~F);
    }

    return H; // Return 32-bit hash value for processing the entire pattern.
}

/*
 * Searches the chain backwards and verifies a match if it scans back to the start of the pattern.
 */
inline void search_chain(int *pos, int *count, unsigned int H, unsigned int V, const unsigned int B[ASIZE],
                         const unsigned char *x, const int m, const unsigned char *y, const unsigned int Hm)
{
    const int end_second_qgram_pos = *pos - m + Q2;
    while (*pos >= end_second_qgram_pos)
    {
        *pos -= Q;
        H = CHAIN_HASH(y, *pos);
        if (!(V & FINGERPRINT(H))) return;
        V = B[H & TABLE_MASK];
    }

    *pos = end_second_qgram_pos - Q;
    if (H == Hm && memcmp(y + *pos - END_FIRST_QGRAM, x, m) == 0) {
        (*count)++;
    }
}

/*
 * Searches for a pattern x of length m in a text y of length n and reports the number of occurrences found.
 */
int search(unsigned char *x, int m, unsigned char *y, int n) {
    if (m < Q) return -1;  // have to be at least Q in length to search.
    unsigned int H, V, B[ASIZE];

    /* Preprocessing */
    BEGIN_PREPROCESSING
    const int MQ1 = m - Q + 1;
    const unsigned int Hm = preprocessing(x, m, B);
    END_PREPROCESSING

    /* Searching */
    BEGIN_SEARCHING
    int count = 0;
    int pos = m - 1;
    while (pos < n) {
        H = CHAIN_HASH(y, pos);
        V = B[H & TABLE_MASK];
        if (V) search_chain(&pos,  &count, H, V, B, x, m, y, Hm);
        pos += MQ1;
    }
    END_SEARCHING

    return count;
}

