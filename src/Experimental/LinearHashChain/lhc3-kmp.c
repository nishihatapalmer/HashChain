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
 * It creates Q chains of hashes from the end of the pattern back to the start.
 *
 * This implementation is written to integrate with the SMART string search benchmarking tool,
 * by Simone Faro, Matt Palmer, Stefano Stefano Scafiti and Thierry Lecroq.
*/

#include "include/define.h"
#include "include/main.h"

/*
 * Alpha - the number of bits in the hash table.
 */
#define ALPHA 11

/*
 * Number of bytes in a q-gram.
 * Chain hash functions defined below must be written to process this number of bytes.
 */
#define	Q     3

/*
 * Functions and calculated parameters.
 * Hash functions must be written to use the number of bytes defined in Q. They scan backwards from the initial position.
 */
#define S                 ((ALPHA) / (Q))                          // Bit shift for each of the chain hash byte components.
#define HASH(x, p, s)     ((((x[p] << (s)) + x[p - 1]) << (s)) + x[p - 2])  // General hash function using a bitshift for each byte added.
#define CHAIN_HASH(x, p)  HASH((x), (p), (S))                      // Hash function for chain hashes, using the S3 bitshift.
#define LINK_HASH(H)      (1U << ((H) & 0x1F))                     // Hash fingerprint, taking low 5 bits of the hash to set one of 32 bits.
#define ASIZE             (1 << (ALPHA))                           // Hash table size.
#define TABLE_MASK        ((ASIZE) - 1)                            // Mask for table is one less than the power of two size.
#define Q2                (Q + Q)                                  // 2 Qs.
#define END_FIRST_QGRAM   (Q - 1)                                  // Position of the end of the first q-gram.
#define END_SECOND_QGRAM  (Q2 - 1)                                 // Position of the end of the second q-gram.

/*
 * Calculates the KMP next table of size m + 1 for a pattern x of length m.
 * It is implemented from the description in "Fast Pattern Matching in Strings" by Knuth, Morris and Pratt, June 1977,
 * with the following differences:
 *  1. Indices are adjusted to be zero-indexed rather than 1-indexed, as in the original paper.
 *  2. The KMP table has an extra entry at position m, so we can continue searching after a match.
 */
void pre_kmp(unsigned char *x, int m, int KMP[])
{
    int j = 0;
    int t = -1;
    KMP[0] = -1;
    while (j < m - 1) {
        while (t > -1 && x[j] != x[t]) {
            t = KMP[t];
        }
        j++;
        t++;
        if (x[j] != x[t]) {
            KMP[j] = KMP[t];
        }
        else {
            KMP[j] = t;
        }
    }
    KMP[m] = t;
}

void preKmp(unsigned char *x, int m, int kmpNext[]) {
    int i, j;
    i = 0;
    j = kmpNext[0] = -1;
    while (i < m) {
        while (j > -1 && x[i] != x[j])
            j = kmpNext[j];
        i++;
        j++;
        if (i<m && x[i] == x[j])
            kmpNext[i] = kmpNext[j];
        else
            kmpNext[i] = j; // i == m, requires m + 1 elements in kmpNext.
    }
}

/*
 * Builds the hash table B of size ASIZE for a string x of length m.
 * Returns the 32-bit hash value of matching the entire pattern.
 */
unsigned int preprocessing(const unsigned char *x, int m, unsigned int *B) {

    // 0. Zero out the hash table.
    for (int i = 0; i < ASIZE; i++) B[i] = 0;

    // 1. Calculate all the chain hashes, ending with processing the entire pattern so H has the cumulative value.
    unsigned int H;
    int start = m < Q2 ? m - END_FIRST_QGRAM : Q;
    for (int chain_no = start; chain_no >= 1; chain_no--)
    {
        H = CHAIN_HASH(x, m - chain_no);
        for (int chain_pos = m - chain_no - Q; chain_pos >= END_FIRST_QGRAM; chain_pos -=Q)
        {
            unsigned int H_last = H;
            H = CHAIN_HASH(x, chain_pos);
            B[H_last & TABLE_MASK] |= LINK_HASH(H);
        }
    }

    // 2. Add in hashes for the first qgrams that have no preceding value.  Only set a value if there is nothing there already.
    unsigned int F;
    int stop = MIN(m, END_SECOND_QGRAM);
    for (int chain_pos = END_FIRST_QGRAM; chain_pos < stop; chain_pos++)
    {
        F = CHAIN_HASH(x, chain_pos);
        if (!B[F & TABLE_MASK]) B[F & TABLE_MASK] = LINK_HASH(~F);
    }

    return H; // Return 32-bit hash value for processing the entire pattern.
}


/*
 * Searches for a pattern x of length m in a text y of length n and reports the number of occurrences found.
 */
int search(unsigned char *x, int m, unsigned char *y, int n) {
    if (m < Q) return -1;  // have to be at least Q in length to search.
    unsigned int H, V, B[ASIZE];
    int KMP[m + 1];

    /* Preprocessing */
    BEGIN_PREPROCESSING
    const int MQ1 = m - Q + 1;
    preprocessing(x, m, B);
    preKmp(x, m, KMP);
    END_PREPROCESSING

    /* Searching */
    BEGIN_SEARCHING
    int count = 0;
    int pos = m - 1;
    int next_verify_pos = 0;
    int pattern_pos = 0;

    // While within the search text:
    while (pos < n) {

        // If there is a bit set for the hash:
        H = CHAIN_HASH(y, pos);
        V = B[H & TABLE_MASK];
        if (V) {

            // Look at the chain of q-grams that precede it:
            const int end_second_qgram_pos = pos - m + Q2;
            const int scan_back_limit = MAX(end_second_qgram_pos, next_verify_pos - 1 + Q);
            while (pos >= scan_back_limit)
            {
                pos -= Q;
                H = CHAIN_HASH(y, pos);
                // If we have no match for this chain q-gram, break out and go around the main loop again:
                if (!(V & LINK_HASH(H))) goto shift;
                V = B[H & TABLE_MASK];
            }

            // Matched the chain all the way back to the start - verify the pattern:
            const int window_start_pos = end_second_qgram_pos - Q2 + 1;

            // Check if we need to re-start KMP if our window start is after last results.
            if (window_start_pos > next_verify_pos) {
                next_verify_pos = window_start_pos;
                pattern_pos = 0;
            }

            //TODO: what does this condition signify?  It is unclear even if theoretically sound...
            //      can it overflow next_verify_pos beyond n?
            while (pattern_pos >= next_verify_pos - window_start_pos) {

                // Naive string matching - how many characters do we match...
                while (pattern_pos < m && x[pattern_pos] == y[next_verify_pos]) {
                    pattern_pos++;
                    next_verify_pos++;
                }

                // If we matched the whole length of the pattern (and we're still inside the text), increase match count.
                if (pattern_pos == m) count++;

                // Get the next matching pattern position.
                pattern_pos = KMP[pattern_pos];
                if (pattern_pos < 0) {
                    pattern_pos++;
                    next_verify_pos++;
                }
            }

            pos = next_verify_pos + m - 1 - pattern_pos;
            continue;
        }

        // Go around the main loop looking for another hash, incrementing the pos by MQ1.
        shift:
        pos += MQ1;
    }
    END_SEARCHING

    return count;
}
