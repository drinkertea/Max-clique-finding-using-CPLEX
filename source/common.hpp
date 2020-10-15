#pragma once

#include "solver.h"
#include "graph.h"

#define NOMINMAX
#include <Windows.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <array>
#include <deque>
#include <functional>
#include <mutex>
#include <atomic>

constexpr size_t g_invalid_index = std::numeric_limits<size_t>::max();

constexpr auto g_strategies = {
    Graph::ColorizationType::default,
    Graph::ColorizationType::maxdegree,
    Graph::ColorizationType::mindegree,
    Graph::ColorizationType::random,
    Graph::ColorizationType::maxdegree_random,
    Graph::ColorizationType::mindegree_random,
};

static double EpsValue(double x)
{
    constexpr double eps = 1e8;
    return std::round(x * eps) / eps;
}

static bool IsInteger(double x)
{
    return EpsValue(x) == std::round(x);
}

static double DiffToInteger(double x)
{
    return std::abs(std::round(x) - x);
}

static std::set<uint32_t> ExtractClique(const std::vector<double>& variables)
{
    std::set<uint32_t> clique;
    for (uint64_t i = 0; i < variables.size(); ++i)
        if (EpsValue(variables[i]) == 1.0)
            clique.emplace(i);
    return clique;
}

struct Solution
{
    double              upper_bound = 0.0;
    std::vector<double> variables;
    size_t              branching_index = g_invalid_index;
    uint64_t            integer_count = 0;
};
