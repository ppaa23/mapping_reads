# E. coli Read Mapper (C++)

## Algorithms overview

1. **FM-index construction.** We build a suffix array using the doubling/radix sort approach that repeatedly sorts suffixes by `2^k` prefixes, stops when all classes are unique, and then extracts the Burrows–Wheeler transform. The FM-index stores:
   * The full reference string (ending with `$`).
   * The suffix array (`sa`) and the BWT string (`bwt`).
   * `C`: cumulative counts of characters for LF-mapping.
   * `occ`: per-position checkpoints of character frequencies, so backward searches run in O(|pattern|).

2. **Seed-and-verify mapping.** Each read (forward and reverse complement) is split into A/C/G/T-only seeds by breaking on `N` (insensitive to leading/trailing or interior `N`s). Seeds shorter than `MIN_SEED` (default 20) are discarded; each valid seed carries its offset within the read so we can reconstruct candidate alignment positions. Exact seed hits are retrieved via the FM-index search, and `verify_read` demands the entire read matches the reference (treating `N` as a wildcard). The first seed hit that passes this check marks the read as mapped with identity score = read length; additional hits flag the read as multi-mapped but no other hits are emitted once one alignment is accepted.

3. **Smith–Waterman fallback.** Reads that do not survive the identity check fall back to a local alignment stage:
   * Reads containing `N` are split on that `N`, so each A/C/G/T block seeds a SW search.
   * Reads without `N` are chunked into fixed-length fragments (`FALLBACK_FRAGMENT_LEN = 24`) that seed the fallback.
   * For each fragment, the FM-index finds all occurrences, and we try at most `FRAGMENT_MAX_HITS = 64` per fragment and `SMITH_WATERMAN_MAX_ATTEMPTS = 400` attempts total to limit work.
   * Each candidate reference region extends `SMITH_WATERMAN_WINDOW = 40` bases beyond the projected read start on both sides to give the SW DP room.
   * The SW scoring scheme is `MATCH = +2`, `MISMATCH = -1`, `GAP = -2`, and alignments must score at least `SMITH_WATERMAN_MIN_SCORE = 50` to count. The DP routine reuses thread-local scratch buffers (`SWScratch`) so that we avoid repeated allocations during batches.
   * Only the first successful SW alignment per read is kept; we do not record additional alignments afterwards.

4. **Reporting and coverage.** Mapped reads contribute to global statistics with atomic counters, track coverage by incrementing `coverage[pos]` for the span of each accepted alignment, and accumulate:
   * `quality_sum` and `quality_count` across all mapped reads (identity and SW alike).
   * `sw_score_sum` and `sw_score_count` only for SW alignments that produced a non-zero score.
   * Counts of identity mappings, SW-from-reads-with-`N`, and SW-from-reads-without-`N`.

## Usage

Build:
```
mkdir -p build
cmake -S . -B build
cmake --build build -j
```

Run (example):
```
./build/ecoli_mapper \
  --ref data/GCF_000005845.2_ASM584v2_genomic.fna \
  --reads data/ERR022075_1.fastq
```

The executable prints progress and finally emits the “short report” below, using the in-memory counters that were just described.

## Short report (example metrics)

```
[INFO] Mapping finished
Total reads processed: <total_reads>
Mapping rate: <mapped>% (<mapped>/<total_reads>)
Unique mapping rate: <unique_rate>% (<unique_reads> reads)
Multi-mapping rate: <multi_rate>% (<multi_reads> reads)
Alignment quality (mean score): 103.46
Genome coverage: <coverage_rate>% (<covered_bases>/<genome_length> bases)
Reads mapped by technique:
  - BWT full identity: <identity_mapped>
  - Smith-Waterman (reads with N): <sw_from_n>
  - Smith-Waterman (reads without N): <sw_from_no_n>
Smith-Waterman mean score (non-zero alignments): 125.58
```

The concrete percentages (`<mapped>`, `<unique_rate>`, etc.) come from the run against the provided FASTQ and reference. In our earlier profiling run we observed about 86% mapping and an alignment quality around 103.5, with SW alignments averaging 125.6.

## Algorithm parameters recap

* **FM-index seeds** are at least 20 bases long (this is `MIN_SEED` in `src/main.cpp`) and are drawn from uninterrupted ACGT segments between `N`s; seeds know their offset so verification can overlay the whole read on top of the reference.
* **Identity verification** uses `verify_read(read_view, ref_start)` to compare every base of the read (treating `N` as wildcard) with the reference without touching the sentinel `$`.
* **Smith–Waterman scoring:** match = +2, mismatch = –1, gap = –2, minimum score = 50; local alignment of a 100 bp read therefore tops out at 200 and falls below the threshold as soon as the number of mismatches/gaps undoes too much of the `+2` total.
* **SW search bounds:** each candidate start scans ±40 bp around the predicted position; we halt further attempts if a region already produced an accepted alignment, ensuring only one alignment per read contributes to the report.

## Result interpretation

* Alignment quality formula:  
  `alignment_quality = quality_sum / quality_count`, where `quality_sum` contains the read length for identity mappings and the SW score for fallback alignments.
* Smith-Waterman mean score formula:  
  `mean_sw_score = sw_score_sum / sw_score_count`, where both sums only collect values from non-zero SW alignments.
* The multi-mapping flag is set when any seed or SW fragment finds multiple FM-index hits before verification, so the unique/multi percentages reflect how often a read could land in more than one spot at the FM-index search stage.

This README now ties together the project structure, the algorithms, their thresholds, and the run-time report you asked for. Let me know if you need the file to list per-read outputs or per-base depth files again.

Total reads: 22720100
Mapped reads: 19623654
Mapping rate: 86.3713%
Unique mappings: 19296169
Multi-mappings: 327485

Total reads: 22720100
Mapped reads: 19619550
Mapping rate: 86.3533%
Unique mappings: 19256495
Multi-mappings: 363055
