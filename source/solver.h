#pragma once

#include <cstdint>
#include <vector>

#ifndef IL_STD
#define IL_STD
#endif
#include <ilcplex/ilocplex.h>

struct Graph;

std::vector<uint32_t> FindMaxCliqueInteger(const Graph& graph);
std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph);
