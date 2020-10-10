#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "ngraph.hpp"

struct Graph
{
    using ColorToVerts = std::map<uint32_t, std::vector<uint32_t>>;

    explicit Graph(const std::string& path);

    explicit Graph(NGraph::tGraph<uint32_t> g) : m_graph(std::move(g)) {}

    const NGraph::tGraph<uint32_t>& get() const { return m_graph; }

    enum class ColorizationType
    {
        default = 0,
        mindegree,
        maxdegree,
    };

    ColorToVerts Colorize(ColorizationType type) const;

private:
    NGraph::tGraph<uint32_t> m_graph;
};
