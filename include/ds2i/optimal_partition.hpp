#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include "util.hpp"

#include "piecewise_linear_model.hpp"
#include "la_vector.hpp"

namespace ds2i {

    typedef uint32_t posting_t;
    typedef uint64_t cost_t;

    using canonical_segment = typename OptimalPiecewiseLinearModel<uint32_t, uint32_t>::CanonicalSegment;

    std::vector<std::vector<canonical_segment>> segmentations(0);
    auto const &conf = configuration::get();

    template<uint8_t hybrid = 1>
    struct optimal_partition {

        std::vector<posting_t> partition;
        cost_t cost_opt = 0; // the costs are in bits!

        template<typename ForwardIterator>
        struct cost_window {
            // a window represents the cost of the interval [start, end)

            ForwardIterator start_it;
            ForwardIterator end_it;
            // starting and ending position of the window
            posting_t start = 0;
            posting_t end = 0; // end-th position is not in the current window
            posting_t min_p = 0; // element that precedes the first element of the window
            posting_t max_p = 0;

            cost_t cost_upper_bound; // The maximum cost for this window

            cost_window(ForwardIterator begin, cost_t cost_upper_bound)
                    : start_it(begin), end_it(begin), min_p(*begin), max_p(0), cost_upper_bound(cost_upper_bound) {}

            uint64_t universe() const {
                return max_p - min_p + 1;
            }

            uint64_t size() const {
                return end - start;
            }

            void advance_start() {
                min_p = *start_it + 1;
                ++start;
                ++start_it;
            }

            void advance_end() {
                max_p = *end_it;
                ++end;
                ++end_it;
            }

        };

        template<typename Iterator>
        void build_segmentation(Iterator begin, uint64_t n, global_parameters const &params) {
            // Precompute all segmentations
            segmentations = std::vector<std::vector<canonical_segment>>(params.max_bpc);
            auto end = begin;
            std::advance(end, n);
            #pragma omp parallel for default(none) firstprivate(n) shared(begin, end, segmentations, params)
            for (uint8_t bpc = 0; bpc < params.max_bpc; ++bpc) {
                if (bpc == 1)
                    continue;
                auto eps = BPC_TO_EPSILON(bpc);
                std::vector<canonical_segment> out;
                out.reserve(eps > 0 ? n / (eps * eps) : n / 8);
                auto in_fun = [begin](auto i) { return std::pair<uint32_t, uint32_t>(i, begin[i]); };
                auto out_fun = [&out](auto cs) { out.push_back(cs); };
                make_segmentation_par(n, eps, in_fun, out_fun);
                segmentations[bpc] = out;
                //segmentations[bpc].emplace_back(n);
                //segmentations[bpc].shrink_to_fit();
            }
        }

        //if we include segments as encoding scheme, we have to add an extra bit
        static const uint64_t type_bits = (hybrid < 2) ? 1 : 2;

        /*template<bool is_vec = false, typename Iterator>
        std::pair<index_type, uint64_t>
        cost_fun(Iterator begin, uint64_t i, uint64_t j, global_parameters const &params) {
            uint64_t _u = 0;
            // is_vec means the encoding start from the 0 or from the first value begin[i]
            if constexpr (is_vec)
                _u = begin[j] + 1;
            else
                _u = begin[j] - begin[i] + 1;
            uint64_t _n = j - i + 1;
            if constexpr (hybrid == 0)
                return std::pair<index_type, uint64_t>(elias_fano,
                                                       compact_elias_fano::bitsize(params, _u, _n) + type_bits +
                                                       conf.fix_cost);
            if constexpr (hybrid == 1) {
                uint64_t best_cost = all_ones_sequence::bitsize(params, _u, _n);
                index_type idx_type = all_ones;

                uint64_t ef_cost = compact_elias_fano::bitsize(params, _u, _n) + type_bits;
                if (ef_cost < best_cost) {
                    best_cost = ef_cost;
                    idx_type = elias_fano;
                }

                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, _u, _n) + type_bits;
                if (rb_cost < best_cost) {
                    best_cost = rb_cost;
                    idx_type = ranked_bitvector;
                }

                return std::pair<index_type, uint64_t>(idx_type, best_cost + conf.fix_cost);

            }
            //Look for the segment (with lower c) covering S[i,j]
            if constexpr (hybrid == 2) {
                assert(!segmentations.empty());
                if (all_ones_sequence::bitsize(params, _u, _n) == 0)
                    return std::pair<index_type, uint64_t>(all_ones, conf.fix_cost);
                uint64_t best_cost = 0;
                uint64_t ef_cost = compact_elias_fano::bitsize(params, _u, _n) + type_bits;
                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, _u, _n) + type_bits;
                if (ef_cost < rb_cost)
                    best_cost = ef_cost;
                else
                    best_cost = rb_cost;
                uint64_t segment_cost = 0;
                for (auto c = 0; c < max_bpc; ++c) {
                    if (c == 1)
                        continue;
                    auto comparator = [](const canonical_segment &seg1, uint64_t i)
                            -> bool { return seg1.get_first_x() <= i; };
                    // finding the first segment which starts just before i
                    auto cur_segment = std::prev(std::lower_bound(segmentations[c].begin(),
                                                                  segmentations[c].end(), i, comparator));
                    auto seg_i = cur_segment->get_first_x();
                    auto seg_j = (cur_segment == std::prev(segmentations[c].end())) ? _n : std::next(cur_segment)->get_first_x();
                    if (seg_i <= i && seg_j > j) {
                        // compute bits for the segment
                        // vector of corrections + slope + intercept (just the offset respect the first value)
                        segment_cost = (c * _n) + (32 + c + BIT_WIDTH(max_bpc)) + type_bits;
                        if (segment_cost < best_cost)
                            return std::pair<index_type, uint64_t>(segment, segment_cost + conf.fix_cost);
                    }
                }
                if (ef_cost < rb_cost) {
                    return std::pair<index_type, uint64_t>(elias_fano, ef_cost + conf.fix_cost);
                }
                return std::pair<index_type, uint64_t>(ranked_bitvector, rb_cost + conf.fix_cost);
            }
        }*/

        std::tuple<index_type, uint64_t, int8_t>
        cost_fun(uint64_t universe, uint64_t n, uint64_t i, global_parameters const &params) {
            if constexpr (hybrid == 0)
                return std::tuple<index_type, uint64_t, uint8_t>(elias_fano,
                                                       compact_elias_fano::bitsize(params, universe, n) + type_bits +
                                                       conf.fix_cost, -1);
            if constexpr (hybrid == 1) {
                uint64_t best_cost = all_ones_sequence::bitsize(params, universe, n);
                index_type idx_type = all_ones;
                if (best_cost == 0)
                    return std::tuple<index_type, uint64_t, uint8_t>(idx_type, best_cost + conf.fix_cost, -1);

                uint64_t ef_cost = compact_elias_fano::bitsize(params, universe, n) + type_bits;
                if (ef_cost < best_cost) {
                    best_cost = ef_cost;
                    idx_type = elias_fano;
                }

                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, universe, n) + type_bits;
                if (rb_cost < best_cost) {
                    best_cost = rb_cost;
                    idx_type = ranked_bitvector;
                }

                return std::tuple<index_type, uint64_t, uint8_t>(idx_type, best_cost + conf.fix_cost, -1);

            }
            //Look for the segment (with lower c) covering S[i,i+n]
            if constexpr (hybrid == 2) {
                assert(!segmentations.empty());
                uint64_t best_cost = all_ones_sequence::bitsize(params, universe, n);
                index_type idx_type = all_ones;
                if (best_cost == 0)
                    return std::tuple<index_type, uint64_t, uint8_t>(idx_type, best_cost + conf.fix_cost, -1);

                uint64_t ef_cost = compact_elias_fano::bitsize(params, universe, n) + type_bits;
                if (ef_cost < best_cost) {
                    best_cost = ef_cost;
                    idx_type = elias_fano;
                }

                uint64_t rb_cost = compact_ranked_bitvector::bitsize(params, universe, n) + type_bits;
                if (rb_cost < best_cost) {
                    best_cost = rb_cost;
                    idx_type = ranked_bitvector;
                }

                uint64_t segment_cost = 0;
                uint64_t j = i + n - 1;
                for (auto c = 0; c < params.max_bpc; ++c) {
                    if (c == 1)
                        continue;
                    auto comparator = [](const canonical_segment &seg1, uint64_t i)
                            -> bool { return seg1.get_first_x() <= i; };
                    // finding the first segment which starts just before i
                    auto cur_segment = std::prev(std::lower_bound(segmentations[c].begin(),
                                                                  segmentations[c].end(), i, comparator));
                    auto seg_i = cur_segment->get_first_x();
                    auto seg_j = (cur_segment == std::prev(segmentations[c].end())) ? n : std::next(
                            cur_segment)->get_first_x();
                    if (seg_i <= i && seg_j > j) {
                        // compute bits for the segment
                        // vector of corrections + slope + intercept (just the offset respect the first value)
                        segment_cost = (c * n) + (32 + c + params.bit_width_max_bpc) + type_bits;
                        if (segment_cost < best_cost)
                            return std::tuple<index_type, uint64_t, uint8_t>(segment, segment_cost + conf.fix_cost, c);
                    }
                }
                return std::tuple<index_type, uint64_t, uint8_t>(idx_type, best_cost + conf.fix_cost, -1);
            }
        }


        optimal_partition() = default;

        template<typename ForwardIterator>
        optimal_partition(ForwardIterator begin, uint64_t universe, uint64_t size, global_parameters const &params,
                          double eps1, double eps2) {
            if constexpr (hybrid == 2) {
                build_segmentation(begin, size, params);
            }
            cost_t single_block_cost = std::get<1>(cost_fun(universe, size, 0, params));
            std::vector<cost_t> min_cost(size + 1, single_block_cost);
            min_cost[0] = 0;

            // create the required window: one for each power of approx_factor
            std::vector<cost_window<ForwardIterator>> windows;
            cost_t cost_lb = std::get<1>(cost_fun(1, 1, 0, params)); // minimum cost
            cost_t cost_bound = cost_lb;
            while (eps1 == 0 || cost_bound < cost_lb / eps1) {
                windows.emplace_back(begin, cost_bound);
                if (cost_bound >= single_block_cost) break;
                cost_bound = cost_bound * (1 + eps2);
            }

            std::vector<posting_t> path(size + 1, 0);
            for (posting_t i = 0; i < size; i++) {
                size_t last_end = i + 1;
                for (auto &window: windows) {

                    assert(window.start == i);
                    while (window.end < last_end) {
                        window.advance_end();
                    }

                    cost_t window_cost;
                    while (true) {
                        window_cost = std::get<1>(cost_fun(window.universe(), window.size(), window.start, params));
                        if ((min_cost[i] + window_cost < min_cost[window.end])) {
                            min_cost[window.end] = min_cost[i] + window_cost;
                            path[window.end] = i;
                        }
                        last_end = window.end;
                        if (window.end == size) break;
                        if (window_cost >= window.cost_upper_bound) break;
                        window.advance_end();
                    }

                    window.advance_start();
                }
            }

            posting_t curr_pos = size;
            while (curr_pos != 0) {
                partition.push_back(curr_pos);
                curr_pos = path[curr_pos];
            }
            std::reverse(partition.begin(), partition.end());
            cost_opt = min_cost[size];
        }
    };

}
