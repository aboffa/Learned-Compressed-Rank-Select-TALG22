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

#ifndef MY_BENCHMARK_ARGUMENTS_PARSER_HPP
#define MY_BENCHMARK_ARGUMENTS_PARSER_HPP

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int parse_argument(int argc, char *argv[], bool &mean_variance, bool &array_bench, bool &ef_bench, bool &rrr_bench, bool &rle_bench,
                    bool &la_vector_epsilon, bool &la_vector_opt, bool &elias_codes, bool &ds2i,
                    bool &entropy_corrections, bool &s18_bench) {
    // Declare the supported options.
    po::options_description desc(
            "This program compares the select time, rank time and bit per integer of different data structures. Usage: \n"
            + std::string(argv[0]) + " <options> <file_1> <file_2> ... <file_n>. \n"
                                     "Files must be generated using the function 'write_data_binary<uint32_t, first_is_size = true>' in the file generate_datasets.cpp in this repo (https://github.com/aboffa/Learned-Rank-Select-ALENEX21). \n"
                                     "Allowed options");

    desc.add_options()
            ("help", "Produce help message")
            ("meva", "Print also mean and variance of the gaps/runs-length")
            ("plai", "Plain 32-bit vector (std::vector) + binary search")
            ("elia", "Elias-Fano (sdsl::sd_vector)")
            ("rrrv", "rrr_vector varying the block size (sdsl::rrr_vector)")
            ("rlev", "rle_vector varying the block size (sdsl::rle_vector)")
            ("lave", "la_vector varying epsilon (la_vector<c>)")
            ("lavo", "la_vector space optimized (la_vector_opt)")
            ("encv", "Gap+encoding + elias delta/elias gamma varying the sample density (sdsl::enc_vector)")
            ("ds2i", "Partitioned Elias Fano (uniform partition ds2i::uniform_partitioned_sequence) No hybrid\n"
                     "Partitioned Elias Fano (uniform partition ds2i::uniform_partitioned_sequence) hybrid\n"
                     "Partitioned Elias Fano (optimal partition ds2i::partitioned_sequence) No hybrid\n"
                     "Partitioned Elias Fano (optimal partition ds2i::partitioned_sequence) hybrid\n"
                     "Partitioned Elias Fano (optimal partition ds2i::partitioned_sequence) hybrid with segments")
            ("s18v", "s18_vector varying block size")
            ("entr",
             "Checks the entropy of the corrections and how sdsl::wt_huff and sdsl::csa_wt are good in exploiting it varying the slope of the segments");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    int num_opt = 1;

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(EXIT_SUCCESS);
    }
    if(argc < 2){
        std::cout << desc << "\n";
        exit(EXIT_FAILURE);
    }
    if (vm.count("meva")) {
        mean_variance = true;
        num_opt++;
    }
    if (vm.count("plai")) {
        array_bench = true;
        num_opt++;
    }
    if (vm.count("elia")) {
        ef_bench = true;
        num_opt++;
    }
    if (vm.count("rrrv")) {
        rrr_bench = true;
        num_opt++;
    }
    if (vm.count("rlev")) {
        rle_bench = true;
        num_opt++;
    }
    if (vm.count("lave")) {
        la_vector_epsilon = true;
        num_opt++;
    }
    if (vm.count("lavo")) {
        la_vector_opt = true;
        num_opt++;
    }
    if (vm.count("encv")) {
        elias_codes = true;
        num_opt++;
    }
    if (vm.count("ds2i")) {
        ds2i = true;
        num_opt++;
    }
    if (vm.count("s18v")) {
        s18_bench = true;
        num_opt++;
    }
    if (vm.count("entr")) {
        entropy_corrections = true;
        num_opt++;
    }
    if(num_opt == 1){
        std::cout << desc << "\n";
        exit(EXIT_FAILURE);
    }
    return num_opt;
}
#endif //MY_BENCHMARK_ARGUMENTS_PARSER_HPP
