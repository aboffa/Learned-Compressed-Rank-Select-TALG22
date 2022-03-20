#pragma once

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

namespace ds2i {

    struct global_parameters {
        global_parameters()
            : ef_log_sampling0(9)
            , ef_log_sampling1(8)
            , rb_log_rank1_sampling(9)
            , rb_log_sampling1(8)
            , log_partition_size(7)
            , max_bpc(11)
            , bit_width_max_bpc(BIT_WIDTH(max_bpc))
        {}

        template <typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (ef_log_sampling0, "ef_log_sampling0")
                (ef_log_sampling1, "ef_log_sampling1")
                (rb_log_rank1_sampling, "rb_log_rank1_sampling")
                (rb_log_sampling1, "rb_log_sampling1")
                (log_partition_size, "log_partition_size")
                ;
        }

        uint8_t ef_log_sampling0;
        uint8_t ef_log_sampling1;
        uint8_t rb_log_rank1_sampling;
        uint8_t rb_log_sampling1;
        uint8_t log_partition_size;
        uint8_t max_bpc;
        uint8_t bit_width_max_bpc;
    };

}
