#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

template <class  T>
std::set<T> operator*(const std::set<T>& A, const std::set<T>& B)
{
    std::set<T> res;

    std::set_intersection(A.begin(), A.end(), B.begin(), B.end(),
        inserter(res, res.begin()));

    return res;
}

struct Graph
{
    using ColorToVerts = std::map<uint32_t, std::vector<uint32_t>>;

    explicit Graph(const std::string& path);

    explicit Graph(size_t size)
    {
        m_graph.resize(size, std::vector<bool>(size, false));
        m_adj.resize(size);
    }

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
            for (auto v1 : verts)
            {
                if (v != v1 && !HasEdge(v, v1))
                    return false;
            }
        }
        return true;
    }

    size_t GetSize() const
    {
        return m_graph.size();
    }

    uint32_t GetDegree(uint32_t vertex) const
    {
        return m_adj[vertex].size();
    }

    Graph GetSubGraph(uint32_t vertex) const
    {
        Graph sub(m_graph.size());
        for (auto a : m_adj[vertex])
        {
            for (auto b : (m_adj[vertex] * m_adj[a]))
            {
                if (a == b)
                    continue;

                sub.m_graph[a][b] = true;
                sub.m_graph[b][a] = true;

                sub.m_adj[a].emplace(b);
                sub.m_adj[b].emplace(a);
            }
        }

        return sub;
    }

    bool HasEdge(uint32_t i, uint32_t j) const
    {
        return m_graph[i][j];
    }

private:
    std::vector<std::vector<bool>> m_graph;
    std::vector<std::set<uint32_t>> m_adj;
};
