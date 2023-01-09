/*
 * Copyright 2022 Matt Palmer.  All rights reserved.
 *
 * This is an implementation of the HashChain algorithm, (currently unpublished) by Matt Palmer.
 * It is a factor search similar to WFR or the QF family of algorithms.
 *
 * It builds a hash table containing entries for chains of hashes.  Hashes are chained together by
 * placing a fingerprint of the *next* hash into the entry for the *current* hash.  This enables
 * a check for the second hash value to be performed without requiring a second lookup in the hash table.
 *
 * A rolling hash function is used for each iteration of the chain, which has three main effects:
 *
 * 1) It expands the effective alphabet of the pattern, which makes it work better on lower alphabet data.
 * 2) It creates multiple chains of hashes, which increases the load on the hash table.
 *    Each position in the pattern is an "anchor" hash, which is linked to a chain of q-grams back to the start of the pattern.
 * 3) Because it is a rolling hash, eventually the hash values converge on the same sequences, so we only have to process
 *    part of each chain to ensure a complete set of entries all the way back to the start of the pattern in the hash table.
 *
 * This implementation is written to integrate with the SMART string search benchmarking tool,
 * by Simone Faro, Matt Palmer, Stefano Stefano Scafiti and Thierry Lecroq.
 */

#include "include/define.h"
#include "include/main.h"
#include "math.h"

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
#define S2    4

/*
 * Bit shift for each of the chain hash byte components.
 * This is added to the anchor hash, which should already have a fairly good spread of initial values.
 * We find that very low values of this bit shift work best in general, and higher values usually don't.
 */
#define S3    1

/*
 * Number of bytes in a q-gram.
 * Anchor and chain hash functions defined below must be written to process this number of bytes.
 */
#define	Q     3

/*
 * Functions and calculated parameters.
 * Hash functions must be written to use the number of bytes defined in Q. They scan backwards from the initial position.
 */
#define HASH(x, p, s)     ((((x[p] << (s)) + x[p - 1]) << (s)) + x[p - 2])  // General hash function using a bitshift for each byte added.
#define ANCHOR_HASH(x, p) HASH((x), (p), (S1))                      // Hash function for anchor hashes, using the S1 bitshift.
#define CHAIN_HASH(x, p)  HASH((x), (p), (S3))                      // Hash function for chain hashes, using the S3 bitshift.
#define LINK_HASH(H)    (1U << ((H) & 0x1F))                        // Hash fingerprint, taking low 5 bits of the hash to set one of 32 bits.
#define TABLE_MASK        (ASIZE - 1)                               // Mask for table is one less than the power of two size.
#define Q2                (Q + Q)                                   // 2 Qs.
#define END_FIRST_QGRAM   (Q - 1)                                   // Position of the end of the first q-gram.
#define END_SECOND_QGRAM  (Q2 - 1)                                  // Position of the end of the second q-gram.
#define CEIL_DIV(n, d)    ((int) ceil((double) (n) / (d)))          // Returns the integer ceiling division of numerator and denominator.
#define CHAIN_LENGTH      ((CEIL_DIV(log2(ASIZE), (S2)) + 1) * (Q)) // Length required to synchronise with the rolling hash chain.

/*
 * Builds the hash table B of size ASIZE for a string x of length m.
 * Returns the 32-bit hash value of matching the entire pattern.
 */
unsigned int preprocessing(const unsigned char *x, int m, unsigned int *B) {

    // 0. Zero out the hash table.
    for (int i = 0; i < ASIZE; i++) B[i] = 0;

    // 1. Process all the anchor q-grams with q-grams before them.
    unsigned int H;
    for (int anchor_pos = END_SECOND_QGRAM; anchor_pos < m; anchor_pos++) {
        H = ANCHOR_HASH(x, anchor_pos);
        int start_chain = anchor_pos - Q;
        int stop_chain = MAX(END_FIRST_QGRAM, start_chain - CHAIN_LENGTH);
        for (int chain_pos = start_chain; chain_pos >= stop_chain; chain_pos -= Q) {
            unsigned int H_last = H;
            H = (H << S2) + CHAIN_HASH(x, chain_pos);
            B[H_last & TABLE_MASK] |= LINK_HASH(H);
        }
    }

    // 2. Process the first q-grams at the start of the pattern that have no preceding q-grams.
    //    There is no q-gram before them that we can calculate a fingerprint of, to put in their hash table entry.
    //    However, there is equally no check on its content other than it not being zero.
    //    If it is currently empty, set it to the fingerprint of the inverse of the current hash value, to avoid pointing back to ourselves.
    int stop = MIN(m, END_SECOND_QGRAM);
    for (int anchor = END_FIRST_QGRAM; anchor < stop; anchor++) {
        H = ANCHOR_HASH(x, anchor);
        if (!(B[H & TABLE_MASK])) B[H & TABLE_MASK] = LINK_HASH(~H);
    }

    // 3. Calculate the 32-bit hash value we check when we need to verify a match.
    //    This is the total 32-bit rolling hash value we would see if processing the entire pattern back to the start.
    int final_pos = m - 1;
    H = ANCHOR_HASH(x, final_pos);
    for (int chain_pos = final_pos - Q; chain_pos >= END_FIRST_QGRAM; chain_pos -=Q)
        H = (H << S2) + CHAIN_HASH(x, chain_pos);

    return H; // Return 32-bit hash value for processing the entire pattern.
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
    // While within the search text:
    while (pos < n) {

        // If there is a bit set for the anchor hash:
        H = ANCHOR_HASH(y, pos);
        V = B[H & TABLE_MASK];
        if (V) {

            // Look at the chain of q-grams that precede it:
            const int end_second_qgram_pos = pos - m + Q2;
            while (pos >= end_second_qgram_pos)
            {
                pos -= Q;
                H = (H << S2) + CHAIN_HASH(y, pos);
                // If we have no match for this chain q-gram, shift and go around the main loop again.
                // C does not have an explicit "while...else" construct, so we implement it here with a goto
                // to break out of the loop, avoiding the subsequent verification stage, to proceed straight to shifting.
                if (!(V & LINK_HASH(H))) goto shift;
                V = B[H & TABLE_MASK];
            }

            // Matched the chain all the way back to the start - verify the pattern if the total hash Hm matches as well:
            pos = end_second_qgram_pos - Q;
            if (H == Hm && memcmp(y + pos - END_FIRST_QGRAM, x, m) == 0) {
                (count)++;
            }
        }

        // Shift by MQ1 and go around the main loop looking for another anchor hash.
        shift:
        pos += MQ1;
    }
    END_SEARCHING

    return count;
}

