# Read Mapper

## Algorithms overview

1. **FM-index construction.** We build a suffix array using the doubling/radix sort approach that repeatedly sorts suffixes by `2^k` prefixes, stops when all classes are unique, and then extracts the Burrows–Wheeler transform. The FM-index stores:
   * The full reference string (ending with `$`).
   * The suffix array (`sa`) and the BWT string (`bwt`).
   * `C`: cumulative counts of characters for LF-mapping.
   * `occ`: per-position checkpoints of nucleotide frequencies, so backward searches run in O(|pattern|).

2. **Reads processing.** We conditionally divide the reads into three classes:
    - reads without N (only ATGC) whose we try to map by standard technique (against FM-index). Noteworthy, not all the reads could be ideally mapped (sequencing errors).
    - reads with N. We divide them into non-intersecting fragments by splitting by N characters. Each of those fragments with length not lesser than `MIN_SEED = 20` - to decrease the calculations) we try to map against FM-index. Each one of the successful mapping we try to extend (in both directions) looking for exact matchings.
    - reads that weren't mapped on the previous steps. The main reason for that - sequencing errors. We use heuristic to map thos - as in the step two, we divide the read into fragments (`FALLBACK_FRAGMENT_LEN = 24`, if contains N, split by N) and try to map new fragments. For successful fragment mapping we extend the aligned sequence, looking not for exact matchings (it wouldn't result in anything new) but using DP Smith-Waterman algorithm for local alignments. To prevent too much computations, we restrict the number of initial seed hits and attempts of DP: `FRAGMENT_MAX_HITS = 64`, `SMITH_WATERMAN_MAX_ATTEMPTS = 400`.

Obviously, we accounted the possibilty of mapping to complementary DNA chain (producing mirror reads and mapping them as well).
To speed up the process of read mapping, we use parallelization and dynamic loading of reads in batches.

3. **Data quality.** To compute the quality of .fastq reads we use only information from .fastq file (this way, the whole mapping process may be scipped). We don't account the quality of separate reads in mapping - this would require too complex algorithms. So the only purpose of this part of the project is to provide overall sequencing quality - computed as mean of each nucleotide quality score (decoded by Phred quality score - ASCII symbols from '!' to 'I'). As some parameter of the quality of the mapping itself, we use the mean score of Smith-Waterman algorithm (for not ideally matching sequences). Base on this we could compute the overall meaning score (using the fact that ideally matching 100 bp sequences would have SW score 200).

## Program reports

./build/bioinf (mapper)
```
[INFO] Mapping finished
Total reads processed: 22720100
Mapping rate: 99.87% (22690573/22720100)
Unique mapping rate: 98.67% (22387677 reads)
Multi-mapping rate: 1.33% (302896 reads)
Alignment quality (mean score): 103.46
Genome coverage: 99.04% (4597061/4641651 bases)
Reads mapped by technique:
  - BWT full identity: 19619550
  - Smith-Waterman (reads with N): 4102
  - Smith-Waterman (reads without N): 3066921
Smith-Waterman mean score (non-zero alignments): 125.58
```

./build/fastq_quality (quality calculator)
```
[INFO] overall: 22720100 reads, 2272010000 bases, mean quality 36.16
```

./build/bioinf (mapper, old version, without Smith-Waterman - only exact matchings)
```
[INFO] Mapping finished
Total reads: 22720100
Mapped reads: 19619550
Mapping rate: 86.3533%
Unique mappings: 19256495
Multi-mappings: 363055
```

## Results
This way we can clearly present requested by the task values and properties:
1. Algorithms used: BWT, SA, FM-index, custom-made splitting the read into substrings, Smith-Waterman local alignment.
2. Total number of reads processed: 22720100 (100%)
3. Mapping rates: 86.4% (exact matches only), 99.0% (with errors allowed). Noteworthy, the second value depends on the threshold of allowed Smith-Waterman score (in our case, 50/200).
4. Percentage of unique/multi-mapped reads (out of mapped reads): 98.1%/1.9% (exact matches only), 98.7%/1.3% (with errors allowed).
5. Allignment quality: 125.6 (out of 200, mean Smith-Waterman score for non-exactly matched reads); 36.16 (Phred quality score, describes sequencing).
6. Genome read coverage: 99.0% (with errors allowed).

According to Phred quality score value, about 1 in 5000 nucleotides is sequenced incorrectly. This way, Smith-Waterman is almost completely redundant: if we increase it's threshold even up to 180/200 (Phred predicts the most accurate values with threshold about 196/200, the mapping rate will be almost the same as if we remove it completely. 

## Usage

Build:
```
mkdir -p build
cmake -S . -B build
# building the mapper
cmake --build build -j
# building the quality calculator
cmake --build build --target fastq_quality
```

Run (example):
```
# running the mapper
./build/bioinf
# running the quality calculator
./build/fastq_quality <fastq>
```

The executables print progress and the short report with results.
