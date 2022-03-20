#pragma once

#include <stdexcept>

#include "compact_elias_fano.hpp"
#include "compact_ranked_bitvector.hpp"
#include "all_ones_sequence.hpp"
#include "global_parameters.hpp"
#include "util.hpp"

#include "segment_sequence.hpp"

namespace ds2i {

    // hybrid means the type of hybridness
    // 0 -> no hybridness (just Elias-Fano)
    // 1 -> original hybridness (Elias-Fano + Ranked Characteristic Bit Vector + Run Length)

    template<uint8_t hybrid = 1>
    struct indexed_sequence {

        //static const uint64_t type_bits = 1; // all_ones is implicit
        static const uint64_t type_bits = (hybrid < 2) ? 1 : 2;

        static DS2I_FLATTEN_FUNC uint64_t
        bitsize(global_parameters const &params, uint64_t universe, uint64_t n) {
            if constexpr (hybrid == 0)
                return compact_elias_fano::bitsize(params, universe, n) + type_bits;
            if constexpr (hybrid > 0) {
                uint64_t best_cost = all_ones_sequence::bitsize(params, universe, n);

                uint64_t ef_cost = compact_elias_fano::bitsize(params, universe, n) + type_bits;
                if (ef_cost < best_cost) {
                    best_cost = ef_cost;
                }

                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, universe, n) + type_bits;
                if (rb_cost < best_cost) {
                    best_cost = rb_cost;
                }

                return best_cost;
            }
        }

        static DS2I_FLATTEN_FUNC std::pair<index_type, uint64_t>
        bitsize_type(global_parameters const &params, uint64_t u, uint64_t n) {
            if constexpr (hybrid == 0)
                return std::pair(elias_fano, compact_elias_fano::bitsize(params, u, n) + type_bits);
            if constexpr (hybrid > 0) {
                if (all_ones_sequence::bitsize(params, u, n) == 0)
                    return std::pair<index_type, uint64_t>(all_ones, 0);
                uint64_t ef_cost = compact_elias_fano::bitsize(params, u, n) + type_bits;
                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, u, n) + type_bits;
                if (ef_cost < rb_cost)
                    return std::pair<index_type, uint64_t>(elias_fano, ef_cost);
                else
                    return std::pair<index_type, uint64_t>(ranked_bitvector, rb_cost);
            }
        }

        template<typename Iterator>
        static void write(succinct::bit_vector_builder &bvb,
                          Iterator begin,
                          uint64_t universe, uint64_t n,
                          global_parameters const &params,
                          int8_t c = -1) {
            uint64_t best_cost = 0;
            int best_type = 0;

            if constexpr (hybrid == 0) {
                best_type = elias_fano;
                bvb.append_bits(best_type, type_bits);
            }

            if constexpr (hybrid == 1) {
                best_cost = all_ones_sequence::bitsize(params, universe, n);
                best_type = all_ones;

                if (best_cost) {
                    uint64_t ef_cost = compact_elias_fano::bitsize(params, universe, n) + type_bits;
                    if (ef_cost < best_cost) {
                        best_cost = ef_cost;
                        best_type = elias_fano;
                    }

                    uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, universe, n) + type_bits;
                    if (rb_cost < best_cost) {
                        best_cost = rb_cost;
                        best_type = ranked_bitvector;
                    }

                    bvb.append_bits(best_type, type_bits);
                }
            }

            if constexpr (hybrid == 2) {
                if (c != -1) {
                    best_type = segment;
                    bvb.append_bits(best_type, type_bits);
                } else {
                    best_cost = all_ones_sequence::bitsize(params, universe, n);
                    best_type = all_ones;

                    if (best_cost) {
                        uint64_t ef_cost = compact_elias_fano::bitsize(params, universe, n) + type_bits;
                        if (ef_cost < best_cost) {
                            best_cost = ef_cost;
                            best_type = elias_fano;
                        }

                        uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, universe, n) + type_bits;
                        if (rb_cost < best_cost) {
                            best_cost = rb_cost;
                            best_type = ranked_bitvector;
                        }

                        bvb.append_bits(best_type, type_bits);
                    }
                }
            }

            switch (best_type) {
                case elias_fano:
                    compact_elias_fano::write(bvb, begin,
                                              universe, n,
                                              params);
                    break;
                case ranked_bitvector:
                    compact_ranked_bitvector::write(bvb, begin,
                                                    universe, n,
                                                    params);
                    break;
                case all_ones:
                    all_ones_sequence::write(bvb, begin,
                                             universe, n,
                                             params);
                    break;
                case segment:
                    segment_sequence::write(bvb, begin,
                                            universe, n,
                                            params, c);
                    break;
                default:
                    assert(false);
            }
        }

        class enumerator {
        public:

            typedef std::pair<uint64_t, uint64_t> value_type; // (position, value)

            enumerator() {}

            enumerator(succinct::bit_vector const &bv, uint64_t offset,
                       uint64_t universe, uint64_t n,
                       global_parameters const &params) {
                if constexpr (hybrid == 0) {
                    m_type = elias_fano;
                }
                if constexpr (hybrid > 0) {
                    if (all_ones_sequence::bitsize(params, universe, n) == 0) {
                        m_type = all_ones;
                    } else {
                        m_type = index_type(bv.get_word56(offset)
                                            & ((uint64_t(1) << type_bits) - 1));
                    }
                }

                switch (m_type) {
                    case elias_fano:
                        m_ef_enumerator = compact_elias_fano::enumerator(bv, offset + type_bits,
                                                                         universe, n,
                                                                         params);
                        break;
                    case ranked_bitvector:
                        m_rb_enumerator = compact_ranked_bitvector::enumerator(bv, offset + type_bits,
                                                                               universe, n,
                                                                               params);
                        break;
                    case all_ones:
                        m_ao_enumerator = all_ones_sequence::enumerator(bv, offset + type_bits,
                                                                        universe, n,
                                                                        params);
                        break;
                    case segment:
                        m_se_enumerator = segment_sequence::enumerator(bv, offset + type_bits,
                                                                       universe, n,
                                                                       params);
                        break;
                    default:
                        throw std::invalid_argument("Unsupported type");
                }
            }

#define ENUMERATOR_METHOD(RETURN_TYPE, METHOD, FORMALS, ACTUALS)    \
            RETURN_TYPE DS2I_FLATTEN_FUNC METHOD FORMALS              \
            {                                                       \
                switch (__builtin_expect(m_type, elias_fano)) {     \
                case elias_fano:                                    \
                    return m_ef_enumerator.METHOD ACTUALS;          \
                case ranked_bitvector:                              \
                    return m_rb_enumerator.METHOD ACTUALS;          \
                case all_ones:                                      \
                    return m_ao_enumerator.METHOD ACTUALS;          \
                case segment:                                       \
                    return m_se_enumerator.METHOD ACTUALS;          \
                default:                                            \
                    assert(false);                                  \
                    __builtin_unreachable();                        \
                }                                                   \
            }                                                       \
            /**/

            // semicolons are redundant but they are needed to get emacs to
            // align the lines properly
            ENUMERATOR_METHOD(value_type, move, (uint64_t position), (position));

            ENUMERATOR_METHOD(value_type, next_geq, (uint64_t lower_bound), (lower_bound));

            ENUMERATOR_METHOD(value_type, next, (), ());

            ENUMERATOR_METHOD(uint64_t, size, () const, ());

            ENUMERATOR_METHOD(uint64_t, prev_value, () const, ());

#undef ENUMERATOR_METHOD
#undef ENUMERATOR_VOID_METHOD

        private:
            index_type m_type;
            union {
                compact_elias_fano::enumerator m_ef_enumerator;
                compact_ranked_bitvector::enumerator m_rb_enumerator;
                all_ones_sequence::enumerator m_ao_enumerator;
                segment_sequence::enumerator m_se_enumerator;
            };
        };
    };
}
