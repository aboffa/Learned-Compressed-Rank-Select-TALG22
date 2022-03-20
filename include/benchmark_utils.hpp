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

#pragma once

#include <chrono>
#include <iostream>
#include "climits"

#include "array.hpp"
#include <sdsl/enc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <index_types.hpp>
#include <succinct/mapper.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/io.hpp>
#include "wrapper_sdsl.hpp"
#include "la_vector.hpp"
#include "la_vector_best_slope.hpp"
#include "la_vector_slope.hpp"
#include <sdsl/s18_vector.hpp>
#include <sdsl/rle_vector.hpp>

const int start_bpc = 6;
const int end_bpc = 14;
const int step_bpc = 1;

const int start_rrr = 4;
const int end_rrr = 7;
const int step_rrr = 1;

const int start_elias = 4;
const int end_elias = 7;
const int step_elias = 1;

const int start_rle = 32;
const int end_rle = 160;
const int step_rle = 32;

const int start_bs18 = 1;
const int end_bs18 = 5;
const int step_bs18 = 1;

const int start_bpc_entropy = 5;
const int end_bpc_entropy = 7;
const int step_bpc_entropy = 1;

using timer = std::chrono::high_resolution_clock;

template<class T>
void do_not_optimize(T const &value) {
    asm volatile("" : : "r,m"(value) : "memory");
}

template<typename D>
void print_mean_variance(D data) {
    uint64_t sum = std::accumulate(data.begin(), data.end(), (uint64_t) 0);
    double mean = static_cast<double>(sum) / data.size();
    std::vector<double> diff(data.size());
    std::transform(data.begin(), data.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / data.size());
    std::cout << "," << mean << "," << stdev;
}

template<typename V, typename Q>
void test_select_bpk(const V &v, const Q &queries) {
    auto start = timer::now();
    auto cnt = 0;
    for (auto i = 0; i < queries.size(); ++i) {
        cnt += v.select(queries[i]);
    }
    do_not_optimize(cnt);

    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
    std::cout << "," << static_cast<double>(elapsed) / queries.size()
              << "," << v.bits_per_element() << std::flush;
}

template<typename V, typename Q>
void test_rank(const V &v, const Q &queries) {
    auto start = timer::now();
    auto cnt = 0;
    for (auto i = 0; i < queries.size(); ++i) {
        auto res_lower_bound = v.rank(queries[i]);
        if constexpr (std::is_scalar<decltype(res_lower_bound)>::value)
            cnt += res_lower_bound;
        else
            cnt += *res_lower_bound;
    }

    do_not_optimize(cnt);
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
    std::cout << "," << static_cast<double>(elapsed) / queries.size() << std::flush;
}

template<int First, int Last, int step, typename Lambda>
inline void static_for(Lambda const &f) {
    if constexpr (First <= Last) {
        f(std::integral_constant<size_t, First>{});
        static_for<First + step, Last, step>(f);
    }
}

const int rrr_block_size = 63;
const int SA_sample = 32;
const int ISA_sample = 32;
using bit_vector_compression = sdsl::rrr_vector<rrr_block_size>;

template<uint8_t bits_per_correction, bool _csa, bool _wt>
void test_one_entropy_corrections(std::vector<uint16_t> corrections) {
    constexpr uint8_t bit_alignment = (bits_per_correction < 8) ? 8 : 16;
    corrections.shrink_to_fit();
    // + 1 just because sdsl::wt_huff uses 0 as a special value, so we use one bit more and we add 1 to each correction
    sdsl::int_vector<bits_per_correction + 1> iv(corrections.size());
    sdsl::int_vector<bit_alignment> byte_aligned_corrections(corrections.size());
    for (auto i = 0; i < corrections.size(); ++i) {
        corrections[i] += 1;
        byte_aligned_corrections[i] = corrections[i];
        assert(byte_aligned_corrections[i] != 0);
        iv[i] = corrections[i];
    }
    std::cout << "," << (sdsl::size_in_bytes(iv) * 8.0) / iv.size();
    if constexpr (_wt) {
        sdsl::wt_huff<bit_vector_compression> wt;
        sdsl::construct_im(wt, byte_aligned_corrections);
        std::cout << "," << (sdsl::size_in_bytes(wt) * 8.0) / static_cast<double>(wt.size());
    }
    if constexpr (_csa) {
        sdsl::csa_wt<sdsl::wt_huff<bit_vector_compression>, SA_sample, ISA_sample> fm_index;
        sdsl::construct_im(fm_index, byte_aligned_corrections);
        std::cout << "," << (sdsl::size_in_bytes(fm_index) * 8.0) / static_cast<double>(fm_index.size()) << std::flush;
    }
}


template<typename Sequence, class T>
void test_ds2i(const std::vector<T> &dataset, const std::vector<size_t> &rands1, const std::vector<T> &rands2) {
    ds2i::global_parameters params;
    succinct::bit_vector_builder docs_bits;
    T universe = dataset.back() + 1;
    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> tuple1;
    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> tuple2;
    if constexpr (std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<0, ds2i::indexed_sequence<0>>> ||
                  std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<1, ds2i::indexed_sequence<1>>>)
        tuple2 = Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    else
        tuple1 = Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    succinct::bit_vector bit_vector(&docs_bits);
    typename Sequence::enumerator enumerator(bit_vector, 0, universe, dataset.size(), params);
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto val = 0;
        for (size_t i = 0; i < rands1.size(); ++i)
            val += enumerator.move(rands1[i]).second;
        do_not_optimize(val);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed =
                static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) /
                rands1.size();
        double bpk = static_cast<double>(((succinct::mapper::size_tree_of(bit_vector)->size) * 8)) / dataset.size();
        std::cout << "," << elapsed << "," << bpk << std::flush;
    }
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto val = 0;
        for (size_t i = 0; i < rands2.size(); ++i)
            val += enumerator.next_geq(rands2[i]).second;

        do_not_optimize(val);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed =
                static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) /
                rands2.size();
        std::cout << "," << elapsed << std::flush;
        if constexpr (std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<0, ds2i::indexed_sequence<0>>> ||
                      std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<1, ds2i::indexed_sequence<1>>>) {
            // static_cast<double> just to have scientific notation
            std::cout << "," << static_cast<double>(std::get<0>(tuple2)) / dataset.size();
            std::cout << "," << std::get<1>(tuple2);
            std::cout << "," << std::get<2>(tuple2);
            std::cout << "," << std::get<3>(tuple2) << std::flush;
        } else {
            std::cout << "," << static_cast<double>(std::get<0>(tuple1)) / dataset.size();
            std::cout << "," << std::get<1>(tuple1);
            std::cout << "," << std::get<2>(tuple1);
            std::cout << "," << std::get<3>(tuple1);
            std::cout << "," << std::get<4>(tuple1);
            std::cout << "," << std::get<5>(tuple1);
            std::cout << "," << std::get<6>(tuple1);
            std::cout << "," << std::get<7>(tuple1);
            std::cout << "," << std::get<8>(tuple1) << std::flush;
        }
    }
}

template<typename Sequence, class T>
void test_ds2i(const std::vector<T> &dataset) {
    ds2i::global_parameters params;
    succinct::bit_vector_builder docs_bits;
    T universe = dataset.back() + 1;
    if constexpr (std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<0, ds2i::indexed_sequence<0>>> ||
                  std::is_same_v<Sequence, ds2i::uniform_partitioned_sequence<1, ds2i::indexed_sequence<1>>>)
        Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    else
        Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    succinct::bit_vector bit_vector(&docs_bits);
}

template<typename D, typename Qsize, typename Qvalue>
void test_select_bpk_rank(D &d, Qsize &rands1, Qvalue &rands2) {
    test_select_bpk(d, rands1);
    test_rank(d, rands2);
}

template<uint8_t b_size_exp, typename D, typename Qsize, typename Qvalue>
void test_select_bpk_rank_s18(D data, Qsize &rands1, Qvalue &rands2) {
    sdsl::bit_vector bv(data.back() + 1, 0);
    for (auto i = 0; i < data.size(); ++i) {
        bv[data[i]] = 1;
    }
    constexpr auto b_size = 1u << b_size_exp;
    sdsl::s18::vector<b_size> s18v(bv);
    sdsl::s18::rank_support<1, b_size> s18rank(s18v);
    sdsl::s18::select_support<1, b_size> s18select(s18v);
    {
        auto start = timer::now();
        auto cnt = 0;
        for (auto i = 0; i < rands1.size(); ++i) {
            cnt += s18select(rands1[i]);
        }
        do_not_optimize(cnt);

        double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
        std::cout << "," << static_cast<double>(elapsed) / rands1.size()
                  << "," << (sdsl::size_in_bytes(s18v) * 8.0) / data.size() << std::flush;
    }
    {
        auto start = timer::now();
        auto cnt = 0;
        for (auto i = 0; i < rands2.size(); ++i) {
            cnt += s18rank(rands2[i]);
        }

        do_not_optimize(cnt);
        double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
        std::cout << "," << static_cast<double>(elapsed) / rands2.size() << std::flush;
    }
}

template<uint8_t b_size_exp, typename D>
void test_s18(D data) {
    sdsl::bit_vector bv(data.back() + 1, 0);
    for (auto i = 0; i < data.size(); ++i) {
        bv[data[i]] = 1;
    }
    constexpr auto b_size = 1u << b_size_exp;
    sdsl::s18::vector<b_size> s18v(bv);
    sdsl::s18::rank_support<1, b_size> s18rank(s18v);
    sdsl::s18::select_support<1, b_size> s18select(s18v);
}

void print_csv_header(bool mean_variance, bool array, bool ef, bool rrr, bool rle, bool la_vector_epsilon,
                      bool la_vector_opt,
                      bool elias_codes, bool ds2i, bool entropy_corrections, bool s18_bench) {
    std::cout << "filename,n,u,ratio";
    if (mean_variance) {
        std::cout << ",mean_gap,stddev_gap,mean_run,stddev_run";
    }

    if (array) {
        std::cout << ",array_time_select,array_bpk"
                     ",array_time_rank";
    }
    if (ef) {
        std::cout << ",ef_sd_time_select,ef_sd_bpk"
                     ",ef_sd_time_rank";
    }
    if (rrr) {
        static_for<start_rrr, end_rrr, step_rrr>([](auto e) {
            constexpr auto block_size = (1u << e) - 1;
            std::string block_size_str = std::to_string(block_size);
            std::cout << ",rrr_" + block_size_str + "_time_select,rrr_" + block_size_str + "_bpk"
                                                                                           ",rrr_" + block_size_str +
                         "_time_rank";
        });
    }
    if (rle) {
        static_for<start_rle, end_rle, step_rle>([](auto block_size) {
            std::string block_size_str = std::to_string(block_size);
            std::cout << ",rle_" + block_size_str + "_time_select,rle_" + block_size_str + "_bpk"
                                                                                           ",rle_" + block_size_str +
                         "_time_rank";
        });
    }
    if (la_vector_epsilon) {
        static_for<start_bpc, end_bpc, step_bpc>([](uint8_t bits_per_correction) {
            std::string bpc = std::to_string(bits_per_correction);
            std::cout << ",la_vector_" + bpc + "_time_select,la_vector_" + bpc + "_bpk"
                                                                                 ",la_vector_" + bpc + "_time_rank"
                                                                                                       ",la_vector_" +
                         bpc + "_segments";
        });
    }
    if (la_vector_opt) {
        std::cout << ",la_vector_opt_time_select,la_vector_opt_bpk"
                     ",la_vector_opt_time_rank"
                     ",la_vector_opt_segments";
    }
    if (elias_codes) {
        static_for<start_elias, end_elias, step_elias>([](auto e) {
            auto e_str = std::to_string(e);
            std::cout << ",elias_delta_" + e_str + "_time_select,elias_delta_" + e_str + "_bpk,"
                                                                                         "elias_delta_" + e_str +
                         "_time_rank";
            std::cout << ",elias_gamma_" + e_str + "_time_select,elias_gamma_" + e_str + "_bpk,"
                                                                                         "elias_gamma_" + e_str +
                         "_time_rank";
        });
    }
    if (ds2i) {
        std::cout << ",uniform_hyb0_select"
                     ",uniform_hyb0_sequence_bpk"
                     ",uniform_hyb0_rank"
                     ",uniform_hyb0_estimate_bpk"
                     ",uniform_hyb0_ef"
                     ",uniform_hyb0_rcbv"
                     ",uniform_hyb0_allones";

        std::cout << ",uniform_hyb1_select"
                     ",uniform_hyb1_sequence_bpk"
                     ",uniform_hyb1_rank"
                     ",uniform_hyb1_estimate_bpk"
                     ",uniform_hyb1_ef"
                     ",uniform_hyb1_rcbv"
                     ",uniform_hyb1_allones";

        std::cout << ",partitioned_hyb0_select"
                     ",partitioned_hyb0_bpk"
                     ",partitioned_hyb0_rank"
                     ",partitioned_hyb0_estimate_bpk"
                     ",partitioned_hyb0_n_ef"
                     ",partitioned_hyb0_ef"
                     ",partitioned_hyb0_n_rcbv"
                     ",partitioned_hyb0_rcbv"
                     ",partitioned_hyb0_n_allones"
                     ",partitioned_hyb0_allones"
                     ",partitioned_hyb0_n_segments"
                     ",partitioned_hyb0_segments";

        std::cout << ",partitioned_hyb1_select"
                     ",partitioned_hyb1_bpk"
                     ",partitioned_hyb1_rank"
                     ",partitioned_hyb1_estimate_bpk"
                     ",partitioned_hyb1_n_ef"
                     ",partitioned_hyb1_ef"
                     ",partitioned_hyb1_n_rcbv"
                     ",partitioned_hyb1_rcbv"
                     ",partitioned_hyb1_n_allones"
                     ",partitioned_hyb1_allones"
                     ",partitioned_hyb1_n_segments"
                     ",partitioned_hyb1_segments";

        std::cout << ",partitioned_hyb2_select"
                     ",partitioned_hyb2_bpk"
                     ",partitioned_hyb2_rank"
                     ",partitions_hyb2_estimate_bpk"
                     ",partitioned_hyb2_n_ef"
                     ",partitioned_hyb2_ef"
                     ",partitioned_hyb2_n_rcbv"
                     ",partitioned_hyb2_rcbv"
                     ",partitioned_hyb2_n_allones"
                     ",partitioned_hyb2_allones"
                     ",partitioned_hyb2_n_segments"
                     ",partitioned_hyb2_segments";
    }
    if (entropy_corrections) {
        auto slopes = {"s_max_", "s_mid_", "s_best_"};
        auto measures = {"plain_bits_", "wt_bits_", "csa_bits_"};
        std::vector<std::string> bits_per_corrections;
        static_for<start_bpc_entropy, end_bpc_entropy, step_bpc_entropy>([&](auto bits_per_correction) {
            bits_per_corrections.push_back(std::to_string(bits_per_correction));
        });
        for (auto &slope: slopes) {
            for (auto &bpc: bits_per_corrections) {
                for (auto &measure: measures) {
                    std::cout << "," << slope << measure << bpc;
                }
            }
        }
    }
    if (s18_bench) {
        static_for<start_bs18, end_bs18, step_bs18>([](auto e) {
            auto e_str = std::to_string(e);
            std::cout << ",s18_" + e_str + "_select,s18_" + e_str + "_bpk,s18_" + e_str + "_rank";
        });
    }
    std::cout << std::endl;
}
