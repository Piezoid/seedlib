# Seedlib

A k-mer index for mapping reads with inexact seeds tolerating 1 error.

Given kmer sizes b1, b2, b3, it construct an index for searching simultaneously:
* 00: seeds of size b1+b2, exact matches,
* 010: seeds of size b1+b2(±1)+b3 with 1 error tolerated.

Overall, this gives the same sensibility as a search for b1+b2+b3 sized seeds within a levensthein edit distance of 1, thank to the pigeon hole principle.
The specificity is lowered but is partially recovered thant to filtering (LIS/chaining, removal of mapped region with few spurious seeds, low entropy kmers removal from the index).


This can be seen as a dumbed down, allowing only 1 error while being more performant, implementation of:

[Vroland, Christophe, et al. "Approximate search of short patterns with high error rates using the 01⁎0 lossless seeds." Journal of Discrete Algorithms 37 (2016): 3-16.](https://www.sciencedirect.com/science/article/pii/S1570866716300028)

### Compilation

```bash
git clone https://github.com/piezoid/seedlib
mkdir build-seedlib
cmake -S seedlib -B build-seedlib -DCMAKE_BUILD_TYPE=release
cmake --build build-seedlib
```


## Usage

### Index construction:

```./seedlib index -i indexname.idx genome.fa```

You can specify:
 * `--b1_len`, `--b2_len`, `--b3_len` for specifying the size of each part of the seeds. Giving only b1 will set the same size for b2 and b3. Otherwise, the results will differ from the statisticall model of 010 seeds, but might provide interesting properties: for example a wider b2 give more room to the seed's region tolerating an edit.
 * `--complexity`: A floating number between 0 and 4bits wich set the minimal entropy of 2-mers in the seeds selected for the index
 * `--downsample`: Fraction of kmer discarded, most common kmers are discarded first


### Query


```./seedlib query -i indexname.idx reads.fa```

Currently only querying with reads: all the kmer that pass the complexity filter are queried.

Output format is not really friendly at the time: each line is comprised of the <tab> separated values:
 * Query read number (in file orderer),
 * Query read length,
 * Mapped indexed sequence number,
 * Mapped indexed sequence length,
 * 0 for forward mapping, 1 for reverse mapping,
 * '-' A sperator/unused field,
 * A series for `query_position,target_position` pairs

The output is not stabilized and currently doesn't show self-mapped reads, or more insidiously, doesn't show mappings of reads having the same number than the indexed sequence they map to.

## Operating principle

[slides](https://piezoid.github.io/slides_seedlib/)
