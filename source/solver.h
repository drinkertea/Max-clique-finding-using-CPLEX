#pragma once

#include <cstdint>
#include <vector>

#ifndef IL_STD
#define IL_STD
#endif
#include <ilcplex/ilocplex.h>

struct Graph;

std::vector<uint32_t> FindMaxClique(const Graph& graph);
