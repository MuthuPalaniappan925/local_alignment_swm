// Importing Packages
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>

// Structure to hold alignment results
struct AlignmentResult {
    std::vector<std::vector<double>> result_matrix;
    std::vector<std::vector<std::string>> traceback_matrix;
    double score;
    std::pair<int, int> score_pos;
};

// Structure to hold aligned sequences
struct AlignedSequences {
    std::string seq_a_aligned;
    std::string seq_b_aligned;
    std::string matches;
};

// Structure to hold benchmark results
struct BenchmarkResult {
    double alignment_time;
    double traceback_time;
    double score;
    size_t sequence_length;
    AlignedSequences alignment;
};

// ScoringSystem class
class ScoringSystem {
private:
    int match;
    int mismatch;
    int gap;
    std::map<char, double> quality_weights;

    double fuzzy_similarity(char a, char b, char quality_a) {
        auto weight = quality_weights.find(quality_a);
        return (weight != quality_weights.end()) ? weight->second : 1.0;
    }

    double default_scoring(char a, char b, char quality_a) {
        if (a == '-' || b == '-') return gap;
        double fuzzy_score = fuzzy_similarity(a, b, quality_a);
        return (a == b) ? match * fuzzy_score : mismatch * fuzzy_score;
    }

public:
    ScoringSystem(int match = 2, int mismatch = -1, int gap = -2) 
        : match(match), mismatch(mismatch), gap(gap) {
        // Quality Weights Init
        char symbols[] = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
        double value = 0.0;
        for (char c : symbols) {
            quality_weights[c] = value;
            value += 0.01;
        }
    }

    double score(char a, char b, char quality_a) {
        return default_scoring(a, b, quality_a);
    }
};

class SequencesAnalyzer {
private:
    std::string seq_a;
    std::string seq_b;
    std::string quality_a;
    ScoringSystem scoring_sys;
    const std::map<int, std::string> traceback_symbols{{0, "↖"}, {1, "↑"}, {2, "←"}};
    // Traceback
    AlignedSequences perform_traceback(const std::vector<std::vector<double>>& result_matrix,
                                     const std::vector<std::vector<std::string>>& traceback_matrix,
                                     std::pair<int, int> start_pos) {
        AlignedSequences result;
        int row = start_pos.first, col = start_pos.second;

        while (result_matrix[row][col] > 0) {
            const std::string& symbol = traceback_matrix[row][col];

            if (symbol == "↖") {
                result.seq_a_aligned = seq_a[row - 1] + result.seq_a_aligned;
                result.seq_b_aligned = seq_b[col - 1] + result.seq_b_aligned;
                result.matches = (seq_a[row - 1] == seq_b[col - 1] ? "|" : " ") + result.matches;
                row--; col--;
            }
            else if (symbol == "↑") {
                result.seq_a_aligned = seq_a[row - 1] + result.seq_a_aligned;
                result.seq_b_aligned = "-" + result.seq_b_aligned;
                result.matches = " " + result.matches;
                row--;
            }
            else if (symbol == "←") {
                result.seq_a_aligned = "-" + result.seq_a_aligned;
                result.seq_b_aligned = seq_b[col - 1] + result.seq_b_aligned;
                result.matches = " " + result.matches;
                col--;
            }
            else break;
        }

        return result;
    }

public:
    SequencesAnalyzer(const std::string& seq_a, const std::string& seq_b, 
                      const std::string& quality_a, int match = 2, 
                      int mismatch = -1, int gap = -2)
        : seq_a(seq_a), seq_b(seq_b), quality_a(quality_a),
          scoring_sys(match, mismatch, gap) {}

    //SWM
    AlignmentResult smith_waterman_algorithm() {
        int rows = seq_a.length() + 1;
        int cols = seq_b.length() + 1;
        std::vector<std::vector<double>> H(rows, std::vector<double>(cols, 0.0));
        std::vector<std::vector<std::string>> traceback(rows, std::vector<std::string>(cols));
        double max_score = 0.0;
        std::pair<int, int> max_pos = {0, 0};

        for (int row = 1; row < rows; row++) {
            for (int col = 1; col < cols; col++) {
                char a = seq_a[row - 1];
                char b = seq_b[col - 1];
                char qa = (row - 1 < quality_a.length()) ? quality_a[row - 1] : '-';

                double score_diag = H[row - 1][col - 1] + scoring_sys.score(a, b, qa);
                double score_up = H[row - 1][col] + scoring_sys.score(a, '-', qa);
                double score_left = H[row][col - 1] + scoring_sys.score('-', b, qa);
                H[row][col] = std::max({0.0, score_diag, score_up, score_left});
                
                if (H[row][col] == score_diag) traceback[row][col] = traceback_symbols.at(0);
                else if (H[row][col] == score_up) traceback[row][col] = traceback_symbols.at(1);
                else if (H[row][col] == score_left) traceback[row][col] = traceback_symbols.at(2);
                
                if (H[row][col] > max_score) {
                    max_score = H[row][col];
                    max_pos = {row, col};
                }
            }
        }
        return {H, traceback, max_score, max_pos};
    }

    AlignedSequences get_alignment(const AlignmentResult& result) {
        return perform_traceback(result.result_matrix, result.traceback_matrix, result.score_pos);
    }

    // Utility - benchmarking
    static std::string generate_random_sequence(size_t length) {
        static const std::string nucleotides = "ACGT";
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<> dis(0, nucleotides.length() - 1);
        
        std::string sequence;
        sequence.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            sequence += nucleotides[dis(gen)];
        }
        return sequence;
    }

    static std::string generate_random_quality(size_t length) {
        static const std::string quality_chars = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<> dis(0, quality_chars.length() - 1);
        
        std::string quality;
        quality.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            quality += quality_chars[dis(gen)];
        }
        return quality;
    }

    // Benchmark runner
    static std::vector<BenchmarkResult> run_benchmarks(
        const std::vector<size_t>& lengths,
        size_t iterations,
        size_t num_references,
        const std::string& output_file
    ) {
        std::vector<BenchmarkResult> all_results;
        std::ofstream file(output_file);
        
        for (size_t length : lengths) {
            std::cout << "Running benchmarks for length " << length << std::endl;
            
            // Generate reference sequences
            std::vector<std::string> reference_sequences;
            for (size_t i = 0; i < num_references; ++i) {
                reference_sequences.push_back(generate_random_sequence(length));
            }
            
            // Run iterations
            for (size_t iter = 0; iter < iterations; ++iter) {
                std::string seq_a = generate_random_sequence(length);
                std::string quality_a = generate_random_quality(length);
                
                for (const auto& ref_seq : reference_sequences) {
                    // Create analyzer and measure alignment time
                    auto start_align = std::chrono::high_resolution_clock::now();
                    SequencesAnalyzer analyzer(seq_a, ref_seq, quality_a);
                    auto result = analyzer.smith_waterman_algorithm();
                    auto end_align = std::chrono::high_resolution_clock::now();
                    
                    // Measure traceback time
                    auto start_trace = std::chrono::high_resolution_clock::now();
                    auto alignment = analyzer.get_alignment(result);
                    auto end_trace = std::chrono::high_resolution_clock::now();
                    
                    // Calculate timings
                    double align_time = std::chrono::duration<double>(end_align - start_align).count();
                    double trace_time = std::chrono::duration<double>(end_trace - start_trace).count();
                    
                    // Store results
                    BenchmarkResult bench_result{
                        align_time,
                        trace_time,
                        result.score,
                        length,
                        alignment
                    };
                    all_results.push_back(bench_result);
                    
                    // Log results periodically
                    if ((iter + 1) % 10 == 0) {
                        write_benchmark_results(file, bench_result, iter + 1, iterations);
                    }
                }
            }
            
            // Write summary statistics for this length
            write_length_summary(file, all_results, length, iterations * num_references);
        }
        
        return all_results;
    }

private:
    static void write_benchmark_results(
        std::ofstream& file,
        const BenchmarkResult& result,
        size_t current_iter,
        size_t total_iter
    ) {
        file << "Iteration " << current_iter << "/" << total_iter 
             << " - Length " << result.sequence_length << ":\n";
        file << "Alignment Time: " << std::fixed << std::setprecision(6) 
             << result.alignment_time << " sec\n";
        file << "Traceback Time: " << result.traceback_time << " sec\n";
        file << "Score: " << result.score << "\n\n";
        file << "Alignment:\n";
        file << result.alignment.seq_a_aligned << "\n";
        file << result.alignment.matches << "\n";
        file << result.alignment.seq_b_aligned << "\n\n";
        file << std::string(80, '-') << "\n";
    }

    static void write_length_summary(
        std::ofstream& file,
        const std::vector<BenchmarkResult>& results,
        size_t length,
        size_t total_iterations
    ) {
        std::vector<double> align_times, trace_times;
        for (const auto& result : results) {
            if (result.sequence_length == length) {
                align_times.push_back(result.alignment_time);
                trace_times.push_back(result.traceback_time);
            }
        }
        
        double avg_align = std::accumulate(align_times.begin(), align_times.end(), 0.0) / align_times.size();
        double avg_trace = std::accumulate(trace_times.begin(), trace_times.end(), 0.0) / trace_times.size();
        
        file << "\n[SUMMARY] Length " << length << " Statistics:\n";
        file << "Average Alignment Time: " << avg_align << " sec\n";
        file << "Average Traceback Time: " << avg_trace << " sec\n";
        file << "Total Iterations: " << total_iterations << "\n\n";
    }
};


int main(int argc, char* argv[]) {
    // Run automated benchmarks
    std::cout << "Starting benchmark runs...\n";
    
    std::vector<size_t> lengths = {32, 64, 128, 256};
    size_t iterations = 70;
    size_t num_references = 5;
    std::string output_file = "cpp_alignment_benchmark_serial_set.txt";
    
    for (size_t length : lengths) {
        std::cout << "Running benchmark for sequence length: " << length << "...\n";

        auto result = SequencesAnalyzer::run_benchmarks(
            {length}, iterations, num_references, output_file);

        std::cout << "Completed benchmark for sequence length: " << length << ".\n";
    }
    
    std::cout << "All benchmarks completed. Results written to: " << output_file << std::endl;
    
    return 0;
}