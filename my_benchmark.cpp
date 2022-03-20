// This file is part of la_vector <https://github.com/aboffa/Learned-Compressed-Rank-Select-TALG22>.
// Copyright (c) 2022 Antonio Boffa.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <iostream>
#include <cstdint>
#include <vector>
#include <random>

#include "benchmark_utils.hpp"
#include "index_reader.hpp"
#include "arguments_parser.hpp"

int main(int argc, char *argv[]) {
    std::cout << std::scientific;
    std::cout << std::setprecision(3);

    bool mean_variance = false;
    bool array_bench = false;
    bool ef_bench = false;
    bool rrr_bench = false;
    bool rle_bench = false;
    bool la_vector_epsilon = false;
    bool la_vector_opt = false;
    bool elias_codes = false;
    bool ds2i = false;
    bool entropy_corrections = false;
    bool s18_bench = false;

    using type = uint32_t;

    int num_opt = parse_argument(argc, argv, mean_variance, array_bench, ef_bench, rrr_bench, rle_bench,
                                 la_vector_epsilon, la_vector_opt, elias_codes, ds2i, entropy_corrections, s18_bench);

    print_csv_header(mean_variance, array_bench, ef_bench, rrr_bench, rle_bench, la_vector_epsilon, la_vector_opt,
                     elias_codes, ds2i, entropy_corrections, s18_bench);

    for (int h = num_opt; h < argc; ++h) {
        std::string filename_path(argv[h]);
        std::vector<type> data(read_data_binary<type, type>(filename_path, true));
        assert(std::is_sorted(data.begin(), data.end()));
        // universe
        type u = data.back() + 1;
        // number of integers
        size_t n = data.size();
        auto ratio = static_cast<double>(n) / static_cast<double>(u);
        std::size_t found = filename_path.find_last_of('/');
        std::string filename = filename_path.substr(found + 1);
        std::cout << "'" << filename << "'," << n << "," << u << "," << ratio << std::flush;
        if (u == 0) {
            std::cout << ",error_while_reading_dataset" << std::endl << std::flush;
            continue;
        }
        if (mean_variance) {
            // Printing mean and variance of the gaps of the integers in the datasets
            std::vector<type> data_gaps;
            data_gaps.reserve(n - 1);
            std::vector<size_t> data_runs;
            data_runs.reserve(n);
            size_t this_run_size = 1;
            for (auto j = 0; j < n - 1; ++j) {
                data_gaps.push_back(data[j + 1] - data[j]);
                if (data[j] == (data[j + 1] - 1)) {
                    this_run_size++;
                } else {
                    data_runs.push_back(this_run_size);
                    this_run_size = 1;
                }
            }
            data_runs.shrink_to_fit();
            // mean and standard deviation fo the gaps
            print_mean_variance(data_gaps);
            // mean and standard deviation fo the length of the runs of consecutive elements
            print_mean_variance(data_runs);
        }
        // Generating datasets of integers for the queries
        size_t TIMES_TEST = data.size() / 5; // query 20% of the dataset
        std::mt19937 mt1(2323);
        // select query
        std::uniform_int_distribution<size_t> dist1(1, n - 1);
        std::vector<size_t> rands1(TIMES_TEST);
        for (auto i = 0; i < TIMES_TEST; ++i) {
            rands1[i] = (dist1(mt1));
        }
        std::mt19937 mt2(4242);
        // rank query
        std::uniform_int_distribution<type> dist2(data.front(), data.back() - 1);
        std::vector<type> rands2(TIMES_TEST);
        for (auto i = 0; i < TIMES_TEST; ++i) {
            rands2[i] = (dist2(mt2));
        }
        // Start testing data structures
        if (array_bench) {
            array<type> arr(data);
            test_select_bpk_rank(arr, rands1, rands2);
        }
        if (ef_bench) {
            wrapper_sdsl<type, sdsl::sd_vector<>> ef(data);
            test_select_bpk_rank(ef, rands1, rands2);
        }
        if (rrr_bench) {
            static_for<start_rrr, end_rrr, step_rrr>([&data, &rands1, &rands2](auto b_size_exp) {
                constexpr auto block_size = (1u << b_size_exp) - 1;
                wrapper_sdsl<type, sdsl::rrr_vector<block_size>> rrr(data);
                test_select_bpk_rank(rrr, rands1, rands2);
            });
        }
        if (rle_bench) {
            static_for<start_rle, end_rle, step_rle>([&data, &rands1, &rands2](auto block_size) {
                wrapper_sdsl<type, sdsl::rle_vector<block_size>> rle(data);
                test_select_bpk_rank(rle, rands1, rands2);
            });
        }
        if (la_vector_epsilon) {
            static_for<start_bpc, end_bpc, step_bpc>([&](auto bits_per_correction) {
                la_vector<type, bits_per_correction> la_vec(data);
                test_select_bpk_rank(la_vec, rands1, rands2);
                std::cout << "," << la_vec.segments_count();
            });
        }
        if (la_vector_opt) {
            la_vector<type> la_vec_opt(data);
            test_select_bpk_rank(la_vec_opt, rands1, rands2);
            std::cout << "," << la_vec_opt.segments_count();
        }
        if (elias_codes) {
            static_for<start_elias, end_elias, step_elias>([&data, &rands1, &rands2](auto e) {
                constexpr auto dens = 1u << e;
                {
                    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<type>> v(
                            data);
                    test_select_bpk_rank(v, rands1, rands2);
                }
                {
                    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<type>> v(
                            data);
                    test_select_bpk_rank(v, rands1, rands2);
                }
            });
        }
        if (ds2i) {
            // test_ds2i<ds2i::compact_elias_fano, T>(a, rands1, rands2);
            // Template parameter is the level of hybridness
            test_ds2i<ds2i::uniform_partitioned_sequence<0, ds2i::indexed_sequence<0>>, type>(data, rands1, rands2);
            test_ds2i<ds2i::uniform_partitioned_sequence<1, ds2i::indexed_sequence<1>>, type>(data, rands1, rands2);
            test_ds2i<ds2i::partitioned_sequence<0, ds2i::indexed_sequence<0>>, type>(data, rands1, rands2);
            test_ds2i<ds2i::partitioned_sequence<1, ds2i::indexed_sequence<1>>, type>(data, rands1, rands2);
            test_ds2i<ds2i::partitioned_sequence<2, ds2i::indexed_sequence<2>>, type>(data, rands1, rands2);
        }
        if (entropy_corrections) {
            static_for<start_bpc_entropy, end_bpc_entropy, step_bpc_entropy>([&](auto bits_per_correction) {
                // MMAX SLOPE
                la_vector_slope<type, bits_per_correction, 0> la_vec(data);
                test_one_entropy_corrections<bits_per_correction, true, true>(u_cvec.get());
            });
            static_for<start_bpc_entropy, end_bpc_entropy, step_bpc_entropy>([&](auto bits_per_correction) {
                // MMID SLOPE
                la_vector_slope<type, bits_per_correction, 4> la_vec(data);
                test_one_entropy_corrections<bits_per_correction, true, true>(u_cvec.get());
            });
            static_for<start_bpc_entropy, end_bpc_entropy, step_bpc_entropy>([&](auto bits_per_correction) {
                // BEST SLOPE
                la_vector_bs<type, bits_per_correction> la_vec(data);
                test_one_entropy_corrections<bits_per_correction, true, true>(u_cvec.get());
            });
        }
        if (s18_bench) {
            static_for<start_bs18, end_bs18, step_bs18>([&](auto b_size_exp) {
                test_select_bpk_rank_s18<b_size_exp>(data, rands1, rands2);
            });
        }
        std::cout << std::endl;
    }
    return 0;
}
