#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

namespace fs = std::filesystem;

static constexpr size_t kDefaultProgressStep = 1'000'000;
static constexpr size_t kGzReadBuffer = 16 * 1024;

void strip_carriage_return(std::string &line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}

class FastqReader {
public:
    explicit FastqReader(const fs::path &path)
        : path_(path), is_gzip_(path.extension() == ".gz") {
        if (is_gzip_) {
            gz_handle_ = gzopen(path_.c_str(), "rb");
            if (!gz_handle_) {
                throw std::runtime_error("failed to open gzip file");
            }
        } else {
            file_.open(path_, std::ios::binary);
            if (!file_) {
                throw std::runtime_error("failed to open file");
            }
        }
    }

    ~FastqReader() {
        if (gz_handle_) {
            gzclose(gz_handle_);
        }
    }

    FastqReader(const FastqReader &) = delete;
    FastqReader &operator=(const FastqReader &) = delete;

    bool read_line(std::string &line) {
        if (is_gzip_) {
            return read_gz_line(line);
        }
        if (!std::getline(file_, line)) {
            return false;
        }
        strip_carriage_return(line);
        return true;
    }

private:
    bool read_gz_line(std::string &line) {
        line.clear();
        bool got_data = false;
        char buffer[kGzReadBuffer];
        while (true) {
            char *next = gzgets(gz_handle_, buffer, sizeof(buffer));
            if (!next) {
                return got_data;
            }
            got_data = true;
            size_t chunk_len = std::strlen(buffer);
            if (chunk_len == 0) {
                continue;
            }
            line.append(buffer, chunk_len);
            if (buffer[chunk_len - 1] == '\n') {
                if (!line.empty()) {
                    line.pop_back();
                }
                strip_carriage_return(line);
                return true;
            }
            if (gzeof(gz_handle_)) {
                strip_carriage_return(line);
                return true;
            }
        }
    }

    fs::path path_;
    bool is_gzip_;
    std::ifstream file_;
    gzFile gz_handle_ = nullptr;
};

struct FileStats {
    uint64_t reads = 0;
    uint64_t bases = 0;
    uint64_t quality_sum = 0;
};

uint64_t phred33_sum(const std::string &quality_line) {
    uint64_t sum = 0;
    for (unsigned char byte : quality_line) {
        if (byte >= 33) {
            sum += uint64_t(byte - 33);
        }
    }
    return sum;
}

std::optional<FileStats> process_fastq(const fs::path &path, size_t progress_step) {
    FastqReader reader(path);
    std::string header;
    std::string seq;
    std::string plus;
    std::string quality;
    FileStats stats;

    while (reader.read_line(header)) {
        if (!reader.read_line(seq) || !reader.read_line(plus) || !reader.read_line(quality)) {
            std::cerr << "[ERROR] " << path << ": truncated FASTQ record\n";
            return std::nullopt;
        }

        stats.reads += 1;
        stats.bases += quality.size();
        stats.quality_sum += phred33_sum(quality);

        if (progress_step > 0 && stats.reads % progress_step == 0) {
            std::cout << "[INFO] " << path << ": " << stats.reads << " reads processed\n";
            std::cout << std::flush;
        }
    }

    return stats;
}

struct Options {
    std::vector<fs::path> inputs;
    size_t progress_step = kDefaultProgressStep;
};

void print_usage(const char *program_name) {
    std::cout << "Usage: " << program_name << " [--progress-step N] FASTQ...\n";
    std::cout << "Compute mean Phred+33 base qualities from FASTQ/FASTQ.gz files.\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --progress-step N   Emit progress logs every N reads (default "
              << kDefaultProgressStep << ")\n";
}

Options parse_options(int argc, char **argv) {
    Options opts;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--progress-step") {
            if (++i >= argc) {
                throw std::runtime_error("--progress-step requires a value");
            }
            opts.progress_step = std::stoull(argv[i]);
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            opts.inputs.emplace_back(arg);
        }
    }
    return opts;
}

int main(int argc, char **argv) {
    Options opts;
    try {
        opts = parse_options(argc, argv);
    } catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }

    if (opts.inputs.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    uint64_t total_reads = 0;
    uint64_t total_bases = 0;
    uint64_t total_quality = 0;

    for (const auto &path : opts.inputs) {
        if (!fs::exists(path)) {
            std::cerr << "[ERROR] " << path << " does not exist\n";
            return 1;
        }

        auto stats = process_fastq(path, opts.progress_step);
        if (!stats) {
            return 1;
        }

        const auto &file_stats = stats.value();
        total_reads += file_stats.reads;
        total_bases += file_stats.bases;
        total_quality += file_stats.quality_sum;

        double mean_quality =
            file_stats.bases ? double(file_stats.quality_sum) / file_stats.bases : 0.0;
        std::cout << "[INFO] " << path << ": "
                  << file_stats.reads << " reads, "
                  << file_stats.bases << " bases, "
                  << "mean quality " << std::fixed << std::setprecision(2)
                  << mean_quality << "\n";
    }

    if (total_bases > 0) {
        double overall = double(total_quality) / total_bases;
        std::cout << "[INFO] overall: "
                  << total_reads << " reads, "
                  << total_bases << " bases, "
                  << "mean quality " << std::fixed << std::setprecision(2)
                  << overall << "\n";
    }

    return 0;
}
