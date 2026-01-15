#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <atomic>
#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <iomanip>
#include "fm_index.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

std::string read_fasta(const std::string& path) {
    std::ifstream in(path);
    std::string line, s;
    while (std::getline(in, line))
        if (!line.empty() && line[0] != '>')
            s += line;
    return s;
}

inline char comp(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string revcomp(const std::string& s) {
    std::string r;
    r.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) r.push_back(comp(*it));
    return r;
}

constexpr size_t SMITH_WATERMAN_WINDOW = 40;
constexpr size_t SMITH_WATERMAN_MAX_ATTEMPTS = 400;
constexpr int SMITH_WATERMAN_MATCH = 2;
constexpr int SMITH_WATERMAN_MISMATCH = -1;
constexpr int SMITH_WATERMAN_GAP = -2;
constexpr int SMITH_WATERMAN_MIN_SCORE = 50;
constexpr size_t FALLBACK_FRAGMENT_LEN = 24;
constexpr int FRAGMENT_MAX_HITS = 64;

using SeedFragment = std::pair<size_t, std::string_view>;

std::vector<SeedFragment> split_seeds(std::string_view s, size_t min_len) {
    std::vector<SeedFragment> seeds;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == 'N') {
            if (i > start && i - start >= min_len)
                seeds.emplace_back(start, s.substr(start, i - start));
            start = i + 1;
        }
    }
    return seeds;
}

std::vector<SeedFragment> chunk_fragments(std::string_view s, size_t chunk_len) {
    std::vector<SeedFragment> chunks;
    size_t pos = 0;
    while (pos < s.size()) {
        size_t len = std::min(chunk_len, s.size() - pos);
        chunks.emplace_back(pos, s.substr(pos, len));
        pos += chunk_len;
    }
    return chunks;
}

struct SWScratch {
    std::vector<int> dp;
    std::vector<uint8_t> trace;
};

bool smith_waterman(std::string_view read, const std::string& reference, size_t segment_origin, size_t segment_len,
                    int min_score, size_t& out_ref_start, size_t& out_span, int& out_score) {
    if (read.empty() || segment_len == 0) return false;
    if (segment_origin + segment_len > reference.size()) return false;
    size_t rows = read.size() + 1;
    size_t cols = segment_len + 1;
    size_t cells = rows * cols;

    thread_local SWScratch scratch;
    if (scratch.dp.size() < cells) scratch.dp.resize(cells);
    if (scratch.trace.size() < cells) scratch.trace.resize(cells);
    std::fill_n(scratch.dp.data(), cells, 0);
    std::fill_n(scratch.trace.data(), cells, 0);

    int best_score = 0;
    size_t best_i = 0, best_j = 0;
    const char* read_data = read.data();

    for (size_t i = 1; i < rows; ++i) {
        size_t base = i * cols;
        size_t base_prev = base - cols;
        for (size_t j = 1; j < cols; ++j) {
            char a = read_data[i - 1];
            char b = reference[segment_origin + j - 1];
            int match = (a == 'N' || b == 'N' || a == b) ? SMITH_WATERMAN_MATCH : SMITH_WATERMAN_MISMATCH;
            int diag = scratch.dp[base_prev + j - 1] + match;
            int up = scratch.dp[base_prev + j] + SMITH_WATERMAN_GAP;
            int left = scratch.dp[base + j - 1] + SMITH_WATERMAN_GAP;
            int cell = diag;
            uint8_t dir = 0;
            if (up > cell) { cell = up; dir = 1; }
            if (left > cell) { cell = left; dir = 2; }
            if (cell <= 0) { cell = 0; dir = 3; }
            scratch.dp[base + j] = cell;
            scratch.trace[base + j] = dir;
            if (cell > best_score) {
                best_score = cell;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_score < min_score) return false;

    size_t ci = best_i;
    size_t cj = best_j;
    while (ci > 0 && cj > 0 && scratch.dp[ci * cols + cj] > 0) {
        uint8_t dir = scratch.trace[ci * cols + cj];
        if (dir == 0) {
            --ci;
            --cj;
        } else if (dir == 1) {
            --ci;
        } else if (dir == 2) {
            --cj;
        } else {
            break;
        }
    }

    size_t span = best_j - cj;
    if (span == 0) return false;
    out_ref_start = segment_origin + cj;
    out_span = span;
    out_score = best_score;
    return true;
}

struct MappingOutcome {
    bool mapped = false;
    bool multi = false;
    bool identity = false;
    bool sw_from_n = false;
    bool sw_from_no_n = false;
    size_t ref_start = 0;
    size_t span = 0;
    int score = 0;
};

static inline MappingOutcome map_one(const FMIndex& fm, const std::string& seq, size_t min_seed) {
    MappingOutcome outcome;
    auto rc = revcomp(seq);
    auto try_identity = [&](const std::string& r) -> bool {
        auto seeds = split_seeds(r, min_seed);
        for (const auto& seed : seeds) {
            auto iv = fm.search(seed.second);
            int hits = iv.second - iv.first;
            if (hits <= 0) continue;
            if (hits > 1) outcome.multi = true;
            for (int idx = iv.first; idx < iv.second; ++idx) {
                size_t seed_pos = fm.sa[idx];
                if (seed_pos < seed.first) continue;
                size_t read_start = seed_pos - seed.first;
                if (fm.verify_read(r, read_start)) {
                    outcome.mapped = true;
                    outcome.identity = true;
                    outcome.ref_start = read_start;
                    outcome.span = r.size();
                    outcome.score = static_cast<int>(r.size());
                    return true;
                }
            }
        }
        return false;
    };

    if (try_identity(seq) || try_identity(rc)) {
        return outcome;
    }

    bool has_n = (seq.find('N') != std::string::npos);
    auto try_sw = [&](const std::string& r) -> bool {
        auto fragments = has_n ? split_seeds(r, min_seed) : chunk_fragments(r, FALLBACK_FRAGMENT_LEN);
        if (fragments.empty()) return false;
        std::unordered_set<size_t> tried_positions;
        size_t attempts = 0;
        for (const auto& fragment : fragments) {
            if (attempts >= SMITH_WATERMAN_MAX_ATTEMPTS) break;
            auto iv = fm.search(fragment.second);
            int hits = iv.second - iv.first;
            if (hits <= 0) continue;
            int processed = 0;
            for (int idx = iv.first; idx < iv.second && attempts < SMITH_WATERMAN_MAX_ATTEMPTS && processed < FRAGMENT_MAX_HITS; ++idx, ++processed) {
                size_t seed_pos = fm.sa[idx];
                if (seed_pos < fragment.first) continue;
                size_t read_start = seed_pos - fragment.first;
                if (!tried_positions.insert(read_start).second) continue;
                size_t max_ref = fm.text.size() - 1;
                if (read_start >= max_ref) continue;
                size_t seg_start = (read_start > SMITH_WATERMAN_WINDOW) ? read_start - SMITH_WATERMAN_WINDOW : 0;
                size_t seg_end = std::min(max_ref, read_start + r.size() + SMITH_WATERMAN_WINDOW);
                if (seg_start >= seg_end) continue;
                size_t seg_len = seg_end - seg_start;
                size_t aligned_start = 0, aligned_span = 0;
                int sw_score = 0;
                if (smith_waterman(r, fm.text, seg_start, seg_len, SMITH_WATERMAN_MIN_SCORE, aligned_start, aligned_span, sw_score)) {
                    outcome.mapped = true;
                    outcome.ref_start = aligned_start;
                    outcome.span = aligned_span;
                    outcome.score = sw_score;
                    if (has_n) outcome.sw_from_n = true;
                    else outcome.sw_from_no_n = true;
                    return true;
                }
                ++attempts;
            }
        }
        return false;
    };

    for (const auto& r : {seq, rc}) {
        if (try_sw(r)) break;
    }

    return outcome;
}

struct AggregateStats {
    size_t identity_mapped = 0;
    size_t sw_from_n = 0;
    size_t sw_from_no_n = 0;
    uint64_t quality_sum = 0;
    size_t quality_count = 0;
    uint64_t sw_score_sum = 0;
    size_t sw_score_count = 0;
};

int main() {
    const size_t MIN_SEED = 20;
    const size_t LOG_STEP = 100000;
    const size_t BATCH = 200000;

    std::cout << "[INFO] Reading reference\n" << std::flush;
    std::string ref = read_fasta("data/GCF_000005845.2_ASM584v2_genomic.fna");

    std::cout << "[INFO] Building FM-index\n" << std::flush;
    FMIndex fm;
    fm.build(ref);
    std::cout << "[INFO] FM-index built\n" << std::flush;

    std::vector<uint32_t> coverage(ref.size(), 0);
    AggregateStats stats;

    std::ifstream fq("data/ERR022075_1.fastq");
    if (!fq) {
        std::cerr << "[ERROR] FASTQ not found (run from project root)\n";
        return 1;
    }

#ifdef _OPENMP
    std::cout << "[INFO] OpenMP threads: " << omp_get_max_threads() << "\n" << std::flush;
#else
    std::cout << "[INFO] OpenMP not enabled\n" << std::flush;
#endif

    std::atomic<size_t> total{0}, mapped{0}, unique{0}, multi{0};

    std::string id, seq, plus, qual;
    std::vector<std::string> batch;
    batch.reserve(BATCH);

    std::cout << "[INFO] FASTQ streaming started\n" << std::flush;

    while (true) {
        batch.clear();
        for (size_t i = 0; i < BATCH; ++i) {
            if (!std::getline(fq, id)) break;
            if (!std::getline(fq, seq)) break;
            if (!std::getline(fq, plus)) break;
            if (!std::getline(fq, qual)) break;
            batch.push_back(seq);
        }
        if (batch.empty()) break;

        size_t batch_n = batch.size();
        size_t start_total = total.fetch_add(batch_n, std::memory_order_relaxed) + batch_n;

        std::cout << "[INFO] Loaded reads: " << start_total << "\n" << std::flush;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (size_t i = 0; i < batch_n; ++i) {
            auto outcome = map_one(fm, batch[i], MIN_SEED);

            if (outcome.mapped) {
                mapped.fetch_add(1, std::memory_order_relaxed);
                if (outcome.multi) multi.fetch_add(1, std::memory_order_relaxed);
                else unique.fetch_add(1, std::memory_order_relaxed);
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    if (outcome.identity) stats.identity_mapped++;
                    if (outcome.sw_from_n) stats.sw_from_n++;
                    if (outcome.sw_from_no_n) stats.sw_from_no_n++;
                    stats.quality_sum += outcome.score;
                    stats.quality_count += 1;
                    if (outcome.sw_from_n || outcome.sw_from_no_n) {
                        stats.sw_score_sum += outcome.score;
                        stats.sw_score_count += 1;
                    }
                    if (outcome.span > 0 && outcome.ref_start < coverage.size()) {
                        size_t end = std::min(coverage.size(), outcome.ref_start + outcome.span);
                        for (size_t pos = outcome.ref_start; pos < end; ++pos) {
                            coverage[pos]++;
                        }
                    }
                }
            }

            size_t done = start_total - batch_n + i + 1;
            if (done % LOG_STEP == 0) {
#ifdef _OPENMP
                #pragma omp critical
#endif
                std::cout << "[INFO] Mapped reads processed: " << done << "\n" << std::flush;
            }
        }
    }

    size_t t = total.load();
    size_t m = mapped.load();
    size_t u = unique.load();
    size_t mm = multi.load();
    size_t genome_len = coverage.empty() ? 0 : coverage.size() - 1;
    size_t covered_bases = 0;
    for (size_t i = 0; i < genome_len; ++i)
        if (coverage[i] > 0) ++covered_bases;

    double mapping_rate = t ? (100.0 * m / t) : 0.0;
    double unique_rate = m ? (100.0 * u / m) : 0.0;
    double multi_rate = m ? (100.0 * mm / m) : 0.0;
    double coverage_rate = genome_len ? (100.0 * covered_bases / genome_len) : 0.0;
    double alignment_quality = stats.quality_count ? (double)stats.quality_sum / stats.quality_count : 0.0;
    double mean_sw_score = stats.sw_score_count ? (double)stats.sw_score_sum / stats.sw_score_count : 0.0;

    std::cout << "[INFO] Mapping finished\n" << std::flush;
    std::cout << "Total reads processed: " << t << "\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Mapping rate: " << mapping_rate << "% (" << m << "/" << t << ")\n";
    std::cout << "Unique mapping rate: " << unique_rate << "% (" << u << " reads)\n";
    std::cout << "Multi-mapping rate: " << multi_rate << "% (" << mm << " reads)\n";
    std::cout << "Alignment quality (mean score): " << alignment_quality << "\n";
    std::cout << "Genome coverage: " << coverage_rate << "% (" << covered_bases << "/" << genome_len << " bases)\n";
    std::cout << "Reads mapped by technique:\n"
              << "  - BWT full identity: " << stats.identity_mapped << "\n"
              << "  - Smith-Waterman (reads with N): " << stats.sw_from_n << "\n"
              << "  - Smith-Waterman (reads without N): " << stats.sw_from_no_n << "\n";
    std::cout << "Smith-Waterman mean score (non-zero alignments): " << mean_sw_score << "\n";
}
