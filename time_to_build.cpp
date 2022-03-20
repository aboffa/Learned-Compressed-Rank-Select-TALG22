#include <iostream>
#include <cstdint>
#include <vector>

#include "benchmark_utils.hpp"
#include "index_reader.hpp"
#include "arguments_parser.hpp"

using type = uint32_t;

using timer = std::chrono::high_resolution_clock;

void print_csv_header_time_to_build() {
    std::cout << "filename,n,u,ratio";
    std::cout << ",array_construction_time";
    std::cout << ",ef_sd_construction_time";
    static_for<start_rrr, end_rrr, step_rrr>([](auto e) {
        constexpr auto block_size = (1u << e) - 1;
        std::string block_size_str = std::to_string(block_size);
        std::cout << ",rrr_" + block_size_str + "_construction_time";
    });
    static_for<start_rle, end_rle, step_rle>([](auto block_size) {
        std::string block_size_str = std::to_string(block_size);
        std::cout << ",rle_" + block_size_str + "_construction_time";
    });

    static_for<start_bpc, end_bpc, step_bpc>([](uint8_t bits_per_correction) {
        std::string bpc = std::to_string(bits_per_correction);
        std::cout << ",la_vector_" + bpc + "_construction_time";
    });

    std::cout << ",la_vector_opt_construction_time";
    static_for<start_elias, end_elias, step_elias>([](auto e) {
        auto e_str = std::to_string(e);
        std::cout << ",elias_delta_" + e_str + "_construction_time";
        std::cout << ",elias_gamma_" + e_str + "_construction_time";
    });

    std::cout << ",uniform_hyb0_construction_time";
    std::cout << ",uniform_hyb1_construction_time";
    std::cout << ",partitioned_hyb0_construction_time";
    std::cout << ",partitioned_hyb1_construction_time";
    std::cout << ",partitioned_hyb2_construction_time";

    static_for<start_bs18, end_bs18, step_bs18>([](auto e) {
        auto e_str = std::to_string(e);
        std::cout << ",s18_" + e_str + "_construction_time";
    });
    std::cout << std::endl;
}


class timer_class {
    timer::time_point start;

public:
    timer_class() : start(timer::now()) {}

    void print_time() {
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start);
        std::cout << "," << elapsed.count() << std::flush;
    };
};

int main(int argc, char *argv[]) {
    std::cout << std::scientific;
    std::cout << std::setprecision(3);

    print_csv_header_time_to_build();
    for (int h = 1; h < argc; ++h) {
        std::string filename_path(argv[h]);
        std::vector<type> data(read_data_binary<type, type>(filename_path, true));
        assert(std::is_sorted(data.begin(), data.end()));
        type u = data.back() + 1;
        size_t n = data.size();
        auto ratio = static_cast<double>(n) / static_cast<double>(u);
        std::size_t found = filename_path.find_last_of('/');
        std::string filename = filename_path.substr(found + 1);
        std::cout << "'" << filename << "'," << n << "," << u << "," << ratio;
        if (u == 0) {
            std::cout << ",error_while_reading_dataset" << std::endl << std::flush;
            continue;
        }
        {
            timer_class t;
            array<type> arr(data);
            t.print_time();
        }
        {
            timer_class t;
            wrapper_sdsl<type, sdsl::sd_vector<>> ef(data);
            t.print_time();
        }
        static_for<start_rrr, end_rrr, step_rrr>([&data](auto b_size_exp) {
            constexpr auto block_size = (1u << b_size_exp) - 1;
            timer_class t;
            wrapper_sdsl<type, sdsl::rrr_vector<block_size>> rrr(data);
            t.print_time();
        });
        static_for<start_rle, end_rle, step_rle>([&data](auto block_size) {
            timer_class t;
            wrapper_sdsl<type, sdsl::rle_vector<block_size>> rle(data);
            t.print_time();
        });


        static_for<start_bpc, end_bpc, step_bpc>([&](auto bits_per_correction) {
            timer_class t;
            la_vector<type, bits_per_correction> la_vec(data);
            t.print_time();
        });
        {
            timer_class t;
            la_vector<type> la_vec_opt(data);
            t.print_time();
        }
        static_for<start_elias, end_elias, step_elias>([&data](auto e) {
            constexpr auto dens = 1u << e;
            {
                timer_class t;
                wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<type>> v(
                        data);
                t.print_time();
            }
            {
                timer_class t;
                wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<type>> v(
                        data);
                t.print_time();
            }
        });
        // test_ds2i<ds2i::compact_elias_fano, T>(a, rands1, rands2);
        // Template parameter is the level of hybridness
        {
            timer_class t;
            test_ds2i<ds2i::uniform_partitioned_sequence<0, ds2i::indexed_sequence<0>>, type>(data);
            t.print_time();
        }
        {
            timer_class t;
            test_ds2i<ds2i::uniform_partitioned_sequence<1, ds2i::indexed_sequence<1>>, type>(data);
            t.print_time();
        }
        {
            timer_class t;
            test_ds2i<ds2i::partitioned_sequence<0, ds2i::indexed_sequence<0>>, type>(data);
            t.print_time();
        }
        {
            timer_class t;
            test_ds2i<ds2i::partitioned_sequence<1, ds2i::indexed_sequence<1>>, type>(data);
            t.print_time();
        }
        {
            timer_class t;
            test_ds2i<ds2i::partitioned_sequence<2, ds2i::indexed_sequence<2>>, type>(data);
            t.print_time();
        }
        static_for<start_bs18, end_bs18, step_bs18>([&](auto b_size_exp) {
            timer_class t;
            test_s18<b_size_exp>(data);
            t.print_time();
        });
        std::cout << std::endl;
    }
    return 0;
}
