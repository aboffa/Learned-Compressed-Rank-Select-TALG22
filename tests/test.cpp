#include <random>
#include <datasets_generation_utils.hpp>
#include <algorithm>
#include <cstdlib>

#include "gtest/gtest.h"
#include "array.hpp"
#include "la_vector.hpp"
#include "benchmark_utils.hpp"
#include "synthetic.hpp"

template<typename Sequence>
void test_ds2i_select_rank(std::vector <uint32_t> &dataset, array <uint32_t> &arr, uint32_t n, uint32_t u_pow) {
    ds2i::global_parameters params;
    succinct::bit_vector_builder docs_bits;
    uint32_t universe = dataset.back() + 1;
    Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    succinct::bit_vector bit_vector(&docs_bits);
    typename Sequence::enumerator enumerator(bit_vector, 0, universe, dataset.size(), params);
    for (uint64_t i = 1; i < n; ++i) {
        ASSERT_EQ(arr.select(i), enumerator.move(i - 1).second) << i;
    }
    for (uint64_t i = 1; i < arr.select(n - 1); ++i) {
        ASSERT_EQ(arr.rank(i), enumerator.next_geq(i).first) << i;
    }
}

const uint32_t u_pow = 20;
const uint32_t n = 100000;

TEST(Rank_Select_Test_sdsl, uniform) {
    std::srand(42);
    std::mt19937 gen(42);
    std::uniform_int_distribution<uint32_t> distribution(0, (1ul << u_pow) - 1);
    auto data = generate_unique(distribution, gen, n, true);
    array<uint32_t> arr(data);
    wrapper_sdsl<uint32_t, sdsl::sd_vector<>> ef_sd(data);
    wrapper_sdsl<uint32_t, sdsl::rrr_vector<>> rrr(data);
    la_vector<uint32_t, 7> lav(data);
    la_vector<uint32_t> lav_opt(data);
    const auto dens = 8;
    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<uint32_t>> v_delta(data);
    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<uint32_t>> v_gamma(data);
    wrapper_sdsl<uint32_t, sdsl::rle_vector<32>> v_rle(data);
    for (int i = 1; i < n; ++i) {
        auto arr_result = arr.select(i);
        ASSERT_EQ(arr_result, ef_sd.select(i)) << i;
        ASSERT_EQ(arr_result, rrr.select(i)) << i;
        ASSERT_EQ(arr_result, v_gamma.select(i)) << i;
        ASSERT_EQ(arr_result, v_delta.select(i)) << i;
        ASSERT_EQ(arr_result, v_rle.select(i)) << i;
        ASSERT_EQ(arr_result, lav.select(i)) << i;
        ASSERT_EQ(arr_result,  lav_opt.select(i)) << i;
    }

    for (int i = 0; i < arr.select(n - 1); ++i) {
        auto arr_result_rank = arr.rank(i);
        ASSERT_EQ(arr_result_rank, ef_sd.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_gamma.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_delta.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, rrr.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_rle.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, lav.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, lav_opt.rank(i)) << i;
    }
}

TEST(Rank_Select_Test_sdsl, clustered) {
    std::vector<uint32_t> data(n);
    ClusteredDataGenerator cg(42);
    cg.fillClustered(data.begin(), data.end(), 0, (1ul << u_pow) - 1);
    array <uint32_t> arr(data);
    wrapper_sdsl<uint32_t, sdsl::sd_vector<>> ef_sd(data);
    wrapper_sdsl<uint32_t, sdsl::rrr_vector<>> rrr(data);
    la_vector<uint32_t, 7> lav(data);
    la_vector<uint32_t> lav_opt(data);
    const auto dens = 8;
    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<uint32_t>> v_delta(data);
    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<uint32_t>> v_gamma(data);
    wrapper_sdsl<uint32_t, sdsl::rle_vector<32>> v_rle(data);
    for (int i = 1; i < n; ++i) {
        auto arr_result = arr.select(i);
        ASSERT_EQ(arr_result, ef_sd.select(i)) << i;
        ASSERT_EQ(arr_result, rrr.select(i)) << i;
        ASSERT_EQ(arr_result, v_gamma.select(i)) << i;
        ASSERT_EQ(arr_result, v_delta.select(i)) << i;
        ASSERT_EQ(arr_result, v_rle.select(i)) << i;
        ASSERT_EQ(arr_result, lav.select(i)) << i;
        ASSERT_EQ(arr_result,  lav_opt.select(i)) << i;
    }

    for (int i = 0; i < arr.select(n - 1); ++i) {
        auto arr_result_rank = arr.rank(i);
        ASSERT_EQ(arr_result_rank, ef_sd.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_gamma.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_delta.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, rrr.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, v_rle.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, lav.rank(i)) << i;
        ASSERT_EQ(arr_result_rank, lav_opt.rank(i)) << i;
    }
}

TEST(Rank_Select_Test_ds2i, uniform) {
    std::srand(42);
    std::mt19937 gen(42);
    std::uniform_int_distribution <uint32_t> distribution(0, (1ul << u_pow) - 1);
    auto data = generate_unique(distribution, gen, n, true);
    array <uint32_t> arr(data);

    test_ds2i_select_rank<ds2i::compact_ranked_bitvector>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::compact_elias_fano>(data, arr, n, u_pow);

    test_ds2i_select_rank<ds2i::uniform_partitioned_sequence <0>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::uniform_partitioned_sequence <1>>(data, arr, n, u_pow);

    test_ds2i_select_rank<ds2i::partitioned_sequence <0>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::partitioned_sequence <1>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::partitioned_sequence <2>>(data, arr, n, u_pow);
}

TEST(Rank_Select_Test_ds2i, clustered) {
    std::vector<uint32_t> data(n);
    ClusteredDataGenerator cg(42);
    cg.fillClustered(data.begin(), data.end(), 0, (1ul << u_pow) - 1);
    array <uint32_t> arr(data);
    test_ds2i_select_rank<ds2i::compact_ranked_bitvector>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::compact_elias_fano>(data, arr, n, u_pow);

    test_ds2i_select_rank<ds2i::uniform_partitioned_sequence <0>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::uniform_partitioned_sequence <1>>(data, arr, n, u_pow);

    test_ds2i_select_rank<ds2i::partitioned_sequence <0>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::partitioned_sequence <1>>(data, arr, n, u_pow);
    test_ds2i_select_rank<ds2i::partitioned_sequence <2>>(data, arr, n, u_pow);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}