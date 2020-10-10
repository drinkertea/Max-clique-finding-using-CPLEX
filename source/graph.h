#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct Graph
{
    explicit Graph(const std::string& path);

    bool   HasEdge(uint32_t i, uint32_t j) const { return m_graph.at(i).at(j); }
    size_t Size() const { return m_graph.size(); }
    uint64_t GetDegree(uint32_t i) const { return m_degrees[i]; }

private:
    std::vector<std::vector<bool>> m_graph;
    std::vector<uint64_t> m_degrees;
};
