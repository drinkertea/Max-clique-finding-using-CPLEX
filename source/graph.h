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

    enum class ColorizationType
    {
        default = 0,
        mindegree,
        maxdegree,
        random,
        mindegree_random,
        maxdegree_random,
    };

    ColorToVerts Colorize(ColorizationType type) const;

    bool IsClique(const std::set<uint32_t>& verts) const
    {
        for (auto v : verts)
        {
            auto neights = m_graph.in_neighbors(v) + m_graph.out_neighbors(v);
            for (auto v1 : neights)
            {
                if (v != v1 && !neights.count(v1))
                    return false;
            }
        }
        return true;
    }

    size_t GetSize() const
    {
        return m_graph.num_vertices();
    }

    uint32_t GetDegree(uint32_t vertex) const
    {
        return m_graph.in_neighbors(vertex).size() + m_graph.out_neighbors(vertex).size();
    }

    Graph GetSubGraph(uint32_t vertex) const
    {
        return Graph(m_graph.subgraph(m_graph.in_neighbors(vertex) + m_graph.out_neighbors(vertex)));
    }

    bool HasEdge(uint32_t i, uint32_t j) const
    {
        return m_graph.out_neighbors(i).count(j);
    }

private:
    NGraph::tGraph<uint32_t> m_graph;
};
