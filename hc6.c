/*
 * Copyright 2022 Matt Palmer.  All rights reserved.
 *
 * This is an implementation of the HashChain algorithm, (currently unpublished) by Matt Palmer.
 *
 * It is written to integrate with the SMART string search benchmarking tool, by Simone Faro and Thierry Lecroq.
 *
 * It is a factor search similar to WFR or QF family of algorithms.
 * It builds a hash table containing an entry for each of the qgrams and for any prior qgrams in the pattern
 * leading up to it.  Subsequent q-grams are chained to each preceding hash entry by a rolling hash.
 * We put a fingerprint of the next qgram hash value in the preceding (from right to left) qgram hash value.
 * This creates a chain of hash values, where each value has to have a fingerprint of the *next* hash value
 * in the chain.  We can then check whether the next hash probably exists in the chain without incurring another index lookup.
 */

#include "include/define.h"
#include "include/log2.h"
#include "include/main.h"

#define	Q     6     // qgram length.
#define S1    2     // bit shift for each of the anchor hash bytes.
#define S2    4     // rolling hash shift.  smaller values create longer hash chains.
#define S3    1     // bit shift for each of the chain hash bytes.
#define ASIZE 4096 // Must be a power of two size, minimum size is 32.

#define TABLE_MASK (ASIZE - 1) // mask for table is one less than the power of two size.

#define HASH(x, p, s)     ((((((((((x[p] << (s)) + x[p - 1]) << (s)) + x[p - 2]) << (s)) + x[p - 3]) << (s)) + x[p - 4]) << (s)) + x[p - 5])
#define ANCHOR_HASH(x, p) HASH((x), (p), (S1))        // Hash function for anchor hashes, using the S1 bitshift.
#define CHAIN_HASH(x, p)  HASH((x), (p), (S3))        // Hash function for chain hashes, using the S3 bitshift.
#define FINGERPRINT(x)    (1U << ((x) & 0x1F))        // Hash fingerprint, taking low 5 bits of the hash to set one of 32 bits.

/*
 * Builds the hash table B of size ASIZE for a string x of length m.
 */
unsigned int preprocessing(unsigned char *x, int m, unsigned int *B) {
    unsigned int H, H2;

    const int rollingHashLimit = (LOG2(ASIZE) * Q / S2) + (Q * 2); // TODO: FIX - only need one more qgram - but the loop below takes one away so we add two here...
    const int fact = m < rollingHashLimit ? m : rollingHashLimit;

    // 0. Zero out the entries in the table.
    for (int i = 0; i < ASIZE; i++) B[i] = 0;

    // 1. Process all the qgrams that will have a qgram preceding it
    //    and limiting how long the hash chain is by the rolling hash limit.
    for (int anchor = (2*Q) - 1; anchor < m - 1; anchor++) {
        H = ANCHOR_HASH(x, anchor);
        int stop = (anchor - fact + 1) > (Q - 1) ? (anchor - fact + 1) :  (Q - 1);
        for (int chainPos = anchor - Q; chainPos >= stop; chainPos -= Q) {
            H2 = H;
            H = (H << S2) + CHAIN_HASH(x, chainPos);
            B[H2 & TABLE_MASK] |= FINGERPRINT(H);
        }
    }

    // 2. Now continue the hash chain for the last qgram down to q-1 to get the cumulative matching hash value.
    //TODO: we only need to process final chain for a rolling hash limit of 32 bits - this will be a lot smaller than a big string.
    H = ANCHOR_HASH(x, m - 1);
    for (int chainPos = m - 1 - Q; chainPos >= Q - 1; chainPos -= Q) {
        H2 = H;
        H = (H << S2) + CHAIN_HASH(x, chainPos);
        B[H2 & TABLE_MASK] |= FINGERPRINT(H);
    }

    // 3. Process the first Q qgrams.  They do not have a preceding qgram, but still need an entry in the table to
    // indicate their presence, if one does not already exist.  We don't have to set anything if there's already
    // something there after all other processing.  Doing this last ensures we only set additional bits if none
    // are already present for other reasons.
    int stop = m < (2*Q) ? m - 1 : (2*Q) - 1;
    for (int anchor = Q-1; anchor <= stop; anchor++) {
        H2 = ANCHOR_HASH(x, anchor);

        // If there's nothing in this hash table entry already, put in a fingerprint.
        //TODO: the XOR isn't necessarily helpful - work this out better.

        // reversing the bits of the actual hash of the data value it is, so we don't create a succession loop,
        // setting its fingerprint to be itself implies it follows itself.  The value itself doesn't really matter,
        // other than we want it to be randomized to some extent and not to point to itself.
        if (!(B[H2 & TABLE_MASK])) {
            B[H2 & TABLE_MASK] = FINGERPRINT(CHAIN_HASH(x, anchor) ^ 0xFF);
        }
    }

    return H;  // Return the hash value of processing the entire pattern (step 2)
}

/*
 * Searches for a pattern x of length m in a text y of length n and reports the number of occurences found.
 */
int search(unsigned char *x, int m, unsigned char *y, int n) {
    if (m < Q) return -1;  // have to be at least Q in length to search.

    const int MQ  = m - Q;
    const int MQ1 = MQ + 1;
    unsigned int H, Hm, V, B[ASIZE];

    /* Preprocessing */
    BEGIN_PREPROCESSING
    Hm = preprocessing(x, m, B); // Hm is the total hash of processing the entire pattern back from the end to the start.
    END_PREPROCESSING

    /* Searching */
    BEGIN_SEARCHING
    int count = 0;
    int pos = m - 1;
    while (pos < n) {
        H = ANCHOR_HASH(y, pos);
        V = B[H & TABLE_MASK];
        if (V) { // If the hash entry is not empty, we have a potential match for the anchor hash
            const int endFirstQgram = pos - MQ;
            chain:
            // As long as we're not before the end of the second qgram, we can move back Q and then read back Q bytes back from there safely.
            if (pos >= endFirstQgram + Q) {
                pos -= Q;
                H = (H << S2) + CHAIN_HASH(y, pos);
                if (V & FINGERPRINT(H)) {  // if next address fingerprint is in our current value,
                    V = B[H & TABLE_MASK]; // get the next value.
                    goto chain;            // and go round again.
                }
                // Did not find fingerprint of next hash - there is no match.  Drop through to pos += MQ1 and go round the main loop again.
            } else { // We read back as far as we can.  Check that the rolling hash equals Hm and if so, verify a match.
                if (H == Hm && memcmp(y + endFirstQgram - Q + 1, x, m) == 0) {
                    count++;
                }
                pos = endFirstQgram;
            }
        }
        pos += MQ1;
    }
    END_SEARCHING
    return count;
}
