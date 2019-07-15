To construct an index and query it for each input kmer:

```./seedlib index -d indexname test.fa```

TODOs:
 * Better API for Quentin
 * Serialization
 * Add the two references implementations in python, and a test workflow
 * Optimize position decoding in the last step of b3->b2_idx->(b2, pos)
 * Batching queries
