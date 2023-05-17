The tuned sentinel version of Hash Chain uses the trick of copying the pattern just past the end of the text.
The copy of the pattern is called the sentinel.  You have to be able to write to the text buffer and not overflow it.

This lets us implement a fast loop during searching, where we don't have to check the position on each time around it,
because it is guaranteed we will run into the sentinel eventually, stopping the algorithm from running off into the rest of memory.
Most of the time though we will be within the text, so avoiding the position check on each loop iteration gives us a speed-up.

While this is a clever trick we can use to get a speed-up, in many practical scenarios, it is not possible to control the
memory allocation of the text buffers that are being searched, and to write additional data into them.  So this can
only be used where all of those factors can be fully controlled.

It can be faster than straight hash chain.
