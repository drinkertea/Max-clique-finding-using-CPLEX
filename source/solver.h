#pragma once

#include <cstdint>
#include <vector>
#include <atomic>

#ifndef IL_STD
#define IL_STD
#endif
#include <ilcplex/ilocplex.h>

struct Graph;

struct CliqueFinder
{
    std::atomic_bool stop = false;
    uint64_t         branch_count = 0;
    double           average_heuristic_time = 0;
    double           average_loop_time = 0;
    double           average_solve_time = 0;

    std::vector<uint32_t> FindMaxCliqueInteger(const Graph& graph);
    std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph);
};

