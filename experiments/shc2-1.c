/*
 * Copyright 2022 Matt Palmer.  All rights reserved.
 *
 * This is an implementation of the HashChain algorithm, (currently unpublished) by Matt Palmer.
 *
 * It is written to integrate with the SMART string search benchmarking tool, by Simone Faro and Thierry Lecroq.
 *
 * //TODO; unclear explanation.
 * It is a factor search similar to WFR or QF family of algorithms.
 * It builds a hash table containing an entry for each of the qgrams and for any prior qgrams in the pattern
 * leading up to it.  Subsequent q-grams are chained to each preceding hash entry by a rolling hash.
 * We put a fingerprint of the next qgram hash value in the preceding (from right to left) qgram hash value.
 * This creates a chain of hash values, where each value has to have a fingerprint of the *next* hash value
 * in the chain.  We can then check whether the next hash probably exists in the chain without incurring another index lookup.
 */

#include "../include/define.h"
#include "../include/main.h"
#include "math.h"

//TODO: quantify limits on pattern sizes etc. by which these values were derived.
//TODO: validate or caveat any statements of "fact".

/*  Tuning algorithm notes.
 *
 *  - Shorter patterns often benefit most from larger table sizes, as they rely more on blank hash table entries.
 *  - Longer patterns often work better with smaller hash table sizes.  They are not primarily relying on blank
 *    entries for their speed, and can benefit more from having better cache-hits on the hash table.
 */

/*
 * Hash table size.  Must be a power of two, minimum size is 32.
 */
#define ASIZE 2048

/*
 * Bit shift for each of the anchor hash byte components.
 * We want to ensure a reasonably good spread and mixing of initial values over the hash table given Q bytes.
 */
#define S1    3

/*
 * Rolling hash bit-shift.  This shifts the previous hash value over by this number of bits.
 * Lower values give longer hash chains (more entries in the hash table).
 * We find that very long chains do not improve performance, but neither do the shortest chains.
 * Setting to 4 seems to work well.
 */
//#define S2    4

/*
 * Bit shift for each of the chain hash byte components.
 * This is added to the anchor hash, which should already have a fairly good spread of initial values.
 * We find that very low values of this bit shift work best in general, and higher values usually don't.
 */
#define S3    3

/*
 * Number of bytes in a q-gram.
 * Anchor and chain hash functions defined below must be written to process this number of bytes.
 */
#define	Q     2

/*
 * Functions and calculated parameters.
 * Hash functions must be written to use the number of bytes defined in Q. They scan backwards from the initial position.
 */
#define HASH(x, p, s)     ((((x)[(p)]) << (s)) + ((x)[(p) - 1]))   // General hash function using a bitshift for each byte added.
//#define ANCHOR_HASH(x, p) HASH((x), (p), (S1))                     // Hash function for anchor hashes, using the S1 bitshift.
#define CHAIN_HASH(x, p)  HASH((x), (p), (S3))                     // Hash function for chain hashes, using the S3 bitshift.
#define FINGERPRINT(H)    (1U << ((H) & 0x1F))                     // Hash fingerprint, taking low 5 bits of the hash to set one of 32 bits.
#define TABLE_MASK        ((ASIZE) - 1)                            // Mask for table is one less than the power of two size.
#define END_FIRST_QGRAM   ((Q) - 1)                                // Position of the end of the first q-gram.
#define END_SECOND_QGRAM  (2 * (Q) - 1)                            // Position of the end of the second q-gram.
//#define CEIL_DIV(n, d)    ((int) ceil((double) (n) / (d)))         // Returns the integer ceiling division of numerator and denominator.
//#define CHAIN_LENGTH      ((CEIL_DIV(log2(ASIZE), (S2)) + 1) * (Q) // Length required to synchronise with the rolling hash chain.
//#define HM_LENGTH         ((CEIL_DIV(32, (S2)) * Q)                // Length required to obtain Hm value for a 32 bit rolling hash.

//TODO: this is the simplest form of hash chain without a rolling hash and just one hash function.
//      Each hash gets linked to the next qgram via the fingerprint, that's it.

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
 * Searches for a pattern x of length m in a text y of length n and reports the number of occurences found.
 */
int search(unsigned char *x, int m, unsigned char *y, int n) {
    if (m < Q) return -1;  // have to be at least Q in length to search.

    const int MQ  = m - Q;
    const int MQQ = MQ - Q;
    const int MQ1 = MQ + 1;
    unsigned int H, Hm, V, B[ASIZE];

    /* Preprocessing */
    BEGIN_PREPROCESSING
    Hm = preprocessing(x, m, B); // Hm is the rolling hash value we see after processing the entire pattern.
    END_PREPROCESSING

    /* Searching */
    BEGIN_SEARCHING
    int count = 0;
    int pos = m - 1;
    while (pos < n) {

        H = CHAIN_HASH(y, pos);
        V = B[H & TABLE_MASK];

        if (V) { // If the hash entry is not empty, we have a potential match for the anchor hash

            //TODO: second qgram construction any faster?
            const int end_second_qgram_pos = pos - MQQ;
            while (1)
            {
                if (pos >= end_second_qgram_pos) {
                    pos -= Q;
                    H = CHAIN_HASH(y, pos);
                    if (!(V & FINGERPRINT(H))) break;  // no fingerprint - end chain and continue main loop.
                    V = B[H & TABLE_MASK]; // get the next value.
                }
                else // We read back as far as we can.  Check that the rolling hash equals Hm and if so, verify a match.
                {
                    pos = end_second_qgram_pos - Q;
                    if (H == Hm && memcmp(y + pos - Q + 1, x, m) == 0) {
                        count++;
                    }
                    break;
                }
            }
        }

        pos += MQ1;
    }
    END_SEARCHING

    return count;
}
