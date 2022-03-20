#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <climits>

template<typename Dist, typename Gen>
std::vector<typename Dist::result_type> generate_unique(Dist &distribution, Gen &generator, size_t n, bool sorted) {
    using T = typename Dist::result_type;

    if constexpr (std::is_same<Dist, std::uniform_int_distribution<T>>::value) {
        std::vector<T> out(n);
        size_t i = 0;
        size_t u = distribution.max() - distribution.min();
        for (auto k = 0; k < u && i < n; ++k)
            if (generator() % (u - k) < n - i)
                out[i++] = k + distribution.min();

        if (!sorted)
            std::random_shuffle(out.begin(), out.end());

        return out;
    }

    std::unordered_set<T> set;
    set.reserve(n);

    while (set.size() < n)
        set.insert(distribution(generator));
    std::vector<T> out(set.begin(), set.end());

    if (sorted)
        std::sort(out.begin(), out.end());
    return out;
}

struct datasetStats {
    float averageLCP = 0;
    float averageLength = 0;
    int maxSize = 0;
    int minSize = INT_MAX;
    size_t size = 0;
    std::vector<char> chars;
};
