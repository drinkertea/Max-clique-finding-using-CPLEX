#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct Graph
{
    explicit Graph(const std::string& path);

    bool   HasEdge(uint32_t i, uint32_t j) const { return m_graph.at(i).at(j); }
    size_t Size() const { return m_graph.size(); }

private:
    std::vector<std::vector<bool>> m_graph;
};
