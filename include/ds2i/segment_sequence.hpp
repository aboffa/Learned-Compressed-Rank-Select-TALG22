#pragma once

#include "global_parameters.hpp"
#include "util.hpp"

#include "piecewise_linear_model.hpp"

/** Computes (bits_per_correction > 0 ? 2^(bits_per_correction-1) - 1 : 0) without the conditional operator. */
#define BPC_TO_EPSILON(bits_per_correction) (((1ul << (bits_per_correction)) + 1) / 2 - 1)
#define CACHE_LINE_BITS (64 * CHAR_BIT)

namespace ds2i {

    struct segment_sequence {

        struct offsets {
            offsets(uint64_t base_offset,
                    uint64_t universe,
                    uint64_t n,
                    global_parameters const& params,
                    uint8_t c)
                    : universe(universe)
                    , n(n)
                    , c_offset(base_offset + 32)
                    , intercept_offset(c_offset + params.bit_width_max_bpc)
                    , corrections_offset(intercept_offset + c)
                    , end(corrections_offset + (c * n))
            {};

            offsets(){};

            uint64_t universe;
            uint64_t n;

            uint64_t c_offset;
            uint64_t intercept_offset;
            uint64_t corrections_offset;
            uint64_t end;
        };

        /* I need c here too, it is a problem. But I think I don't need this method
         * inline static uint64_t
        bitsize(global_parameters const& params, uint64_t universe, uint64_t n)
        {
            return offsets(0, universe, n, params).end;
        }*/

        template <typename Iterator>
        static void write(succinct::bit_vector_builder& bvb,
                          Iterator begin,
                          uint64_t universe,
                          uint64_t n,
                          global_parameters const& params,
                          uint8_t c)
        {
            assert(c >= 0 && c != 1);
            auto in_fun = [begin](auto i) { return std::pair<uint32_t , uint32_t>(i, begin[i]); };
            auto p = in_fun(0);

            OptimalPiecewiseLinearModel<uint32_t , uint32_t> opt(BPC_TO_EPSILON(c));
            opt.add_point(p.first, p.second);

            //should be <=n?
            for (auto i = 1; i < n; ++i) {
                auto next_p = in_fun(i);
                if (i != 1 && next_p.first == p.first)
                    continue;
                p = next_p;
                bool result = opt.add_point(p.first, p.second);
                assert(result);
            }

            auto[significand, exponent, intercept]  = opt.get_segment().get_fixed_point_segment_PEF(0, n);
            bvb.append_bits(significand, 27);
            bvb.append_bits(exponent, 5);
            bvb.append_bits(c, params.bit_width_max_bpc);
            if (c != 0) {
                bvb.append_bits(std::abs(intercept), c);
                for (auto k = 0; k < n; k++) {
                    int64_t approximation = ((int64_t(k) * significand) >> exponent) + intercept;
                    auto error = begin[k] - approximation;
                    auto correction = uint64_t(error + BPC_TO_EPSILON(c));
                    if (BIT_WIDTH(correction) > c)
                        throw std::overflow_error("Segment correction too large");
                    bvb.append_bits(correction, c);
                }
            }
        }

        class enumerator {
        public:

            typedef std::pair<uint64_t, uint64_t> value_type; // (position, value)

            enumerator(succinct::bit_vector const &bv,
                       uint64_t offset,
                       uint64_t universe,
                       uint64_t n,
                       global_parameters const &params) {
                auto tmp_offset = offset;
                m_significand = bv.get_bits(tmp_offset, 27);
                tmp_offset += 27;
                m_exponent = bv.get_bits(tmp_offset, 5);
                tmp_offset += 5;
                m_c = bv.get_bits(tmp_offset, params.bit_width_max_bpc);
                m_epsilon = BPC_TO_EPSILON(m_c);
                //m_linear_threshold = 2 * CACHE_LINE_BITS / (m_c);
                tmp_offset += params.bit_width_max_bpc;
                // this - may be the problem
                m_intercept = - (int64_t (bv.get_bits(tmp_offset, m_c)));
                tmp_offset += m_c;
                m_corrections = &tmp_offset;
                m_of = *new offsets(offset, universe, n, params, m_c);
                m_bv = &bv;
                m_position = size();
                m_value = m_of.universe;
            }

            value_type move(uint64_t position)
            {
                assert(position <= size());
                m_value = decompress(position);
                m_position = position;
                return value_type(m_position, m_value + 1);
            }

            inline value_type value() const
            {
                return value_type(m_position, m_value);
            }

            value_type next_geq(uint64_t lower_bound)
            {
                assert(lower_bound <= m_of.universe);
                auto[pos, bound] = approximate_position(lower_bound);
                pos = std::clamp(pos, 0ul, size());

                auto lo = pos <= bound ? 0 : pos - bound;
                auto hi = std::min(pos + bound + 1, size());

                /*uint64_t val;
                while ((hi - lo) > m_linear_threshold && lo < hi) {
                    auto mid = lo + (hi - lo) / 2;
                    if ((val = decompress(mid)) < lower_bound)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                while (lo < hi && val < lower_bound)
                    val = decompress(++lo);
                */

                while (lo < hi) {
                    auto mid = lo + (hi - lo) / 2;
                    if (decompress(mid) < lower_bound)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                lo = (lo >= 1) ? --lo : lo;
                m_value = decompress(lo);
                m_position = lo;
                return value_type(m_position, m_value + 1);
            }

            value_type next()
            {
                m_position++;
                m_value = decompress(m_position);
                return value();
            }

            uint64_t size() const
            {
                return m_of.n;
            }

            uint64_t prev_value() const
            {
                if (m_position == 0) {
                    return 0;
                }
                return decompress(m_position - 1);
            }

        private:

            inline int64_t approximate(int64_t i) const {
                return ((int64_t(m_significand) * i) >> m_exponent) + m_intercept;
            }

            std::pair<uint64_t , uint64_t> approximate_position(uint32_t value) const {
                auto significand = int64_t(m_significand == 0 ? 1 : m_significand);
                auto p = ((int64_t(value) - m_intercept) << m_exponent) / significand;
                auto bound = 1 + (m_epsilon << m_exponent) / int64_t(significand);
                return {std::max<int64_t>(0, p), bound};
            }

            uint32_t decompress(size_t i) const {
                uint64_t correction = m_bv->get_bits(m_of.corrections_offset + (i * m_c), m_c);
                return approximate(i) + correction - m_epsilon;
            }

            succinct::bit_vector const* m_bv;
            uint64_t m_position;
            uint64_t m_value;

            offsets m_of;
            int64_t m_significand;
            uint8_t m_exponent;
            int64_t m_intercept;
            uint8_t m_c;
            uint64_t m_epsilon;
            uint64_t * m_corrections;
            uint64_t m_linear_threshold;

        };
    };
}
