#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <atomic>
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

std::vector<std::string> split_seeds(const std::string& s, size_t min_len) {
    std::vector<std::string> seeds;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == 'N') {
            if (i > start && i - start >= min_len) seeds.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return seeds;
}

static inline void map_one(const FMIndex& fm, const std::string& seq, size_t min_seed, bool& hit, bool& multi_hit) {
    hit = false;
    multi_hit = false;

    auto rc = revcomp(seq);

    for (const auto& r : {seq, rc}) {
        auto seeds = split_seeds(r, min_seed);
        for (const auto& seed : seeds) {
            auto iv = fm.search(seed);
            int hits = iv.second - iv.first;
            if (hits > 0) {
                hit = true;
                if (hits > 1) multi_hit = true;
            }
        }
    }
}

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
            bool hit = false, multi_hit = false;
            map_one(fm, batch[i], MIN_SEED, hit, multi_hit);

            if (hit) {
                mapped.fetch_add(1, std::memory_order_relaxed);
                if (multi_hit) multi.fetch_add(1, std::memory_order_relaxed);
                else unique.fetch_add(1, std::memory_order_relaxed);
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

    std::cout << "[INFO] Mapping finished\n" << std::flush;
    std::cout << "Total reads: " << t << "\n";
    std::cout << "Mapped reads: " << m << "\n";
    std::cout << "Mapping rate: " << (t ? (100.0 * (double)m / (double)t) : 0.0) << "%\n";
    std::cout << "Unique mappings: " << u << "\n";
    std::cout << "Multi-mappings: " << mm << "\n";
}
