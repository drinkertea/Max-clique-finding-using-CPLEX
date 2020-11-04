#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <functional>

template <class  T>
std::set<T> operator*(const std::set<T>& A, const std::set<T>& B)
{
    std::set<T> res;

    std::set_intersection(A.begin(), A.end(), B.begin(), B.end(),
        inserter(res, res.begin()));

    return res;
}

template <class T>
std::set<T>& operator+=(std::set<T>& A, const std::set<T>& B)
{
    // A.insert(B.begin(), B.end());
    for (typename std::set<T>::const_iterator p = B.begin(); p != B.end(); p++)
        A.insert(*p);

    return A;
}

template <class T>
std::set<T>& operator*=(std::set<T>& A, const std::set<T>& B)
{
    std::erase_if(A, [&B](const auto& x) {
        return !B.count(x);
    });
    return A;
}

struct WeightNode
{
    WeightNode(uint32_t l, double w, uint32_t ro)
        : label(l)
        , weight(w)
        , eps100w(uint32_t(w * 100))
        , random_order(ro)
    {
    }

    bool operator<(const WeightNode& r) const
    {
        return std::tie(eps100w, random_order, weight, label) >
               std::tie(r.eps100w, r.random_order, r.weight, r.label);
    }

    uint32_t label = 0;
    double   weight = 0.0;

private:
    uint32_t eps100w = 0;
    uint32_t random_order = 0;
};

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

    std::set<std::vector<uint32_t>> CheckSolution(const std::set<uint32_t>& verts) const
    {
        std::set<std::vector<uint32_t>> res;
        for (auto v : verts)
        {
            for (auto v1 : verts)
            {
                if (v != v1 && !HasEdge(v, v1))
                {
                    res.emplace(std::vector<uint32_t>{ v, v1 });
                }
            }
        }
        return res;
    }

    size_t GetSize() const
    {
        return m_graph.size();
    }

    uint32_t GetDegree(uint32_t vertex) const
    {
        return m_adj[vertex].size();
    }

    const std::set<uint32_t>& GetNeighbors(uint32_t vertex) const
    {
        return m_adj[vertex];
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

    std::vector<uint32_t> GetOrderedNodes(ColorizationType type) const;

    std::set<std::set<uint32_t>> GetHeuristicConstr(ColorizationType type) const
    {
        return GetHeuristicConstr(GetOrderedNodes(type));
    }

    std::set<std::set<uint32_t>> GetHeuristicConstr(const std::vector<uint32_t>& ordered_nodes) const
    {
        std::vector<uint32_t> nodes_order(m_graph.size());
        uint32_t order = 0;
        for (auto node : ordered_nodes)
            nodes_order[node] = order++;

        struct Node
        {
            uint32_t val = 0;
            const std::vector<uint32_t>& order;

            bool operator<(const Node& r) const
            {
                return std::tie(order[val], val) < std::tie(order[r.val], r.val);
            }
        };

        std::vector<std::set<Node>> m_non_adj;
        m_non_adj.resize(m_graph.size());
        for (uint32_t i = 0; i < m_graph.size(); ++i)
        {
            for (uint32_t j = 0; j < m_graph.size(); ++j)
            {
                if (i == j || m_graph[i][j])
                    continue;

                m_non_adj[i].emplace(Node{ j, nodes_order });
                m_non_adj[j].emplace(Node{ i, nodes_order });
            }
        }

        std::set<Node> nodes;
        for (uint32_t i = 0; i < m_graph.size(); ++i)
            nodes.emplace(Node{ i, nodes_order });

        std::set<std::set<uint32_t>> res;
        while (!nodes.empty())
        {
            std::set<Node> constr;
            constr.emplace(Node{ nodes.begin()->val, nodes_order });
            auto t = m_non_adj[nodes.begin()->val];
            nodes.erase(nodes.begin());
            while (!t.empty())
            {
                auto first = t.begin()->val;
                t = t * m_non_adj[first];
                constr.emplace(Node{ first, nodes_order });
            }

            if (constr.size() < 4)
                continue;

            std::set<uint32_t> converted;
            for (const auto& node : constr)
                converted.emplace(node.val);

            res.emplace(converted);
        }
        return res;
    }

    std::vector<std::set<WeightNode>> GetWeightlyNonAdj(uint32_t start, const std::vector<double>& weights) const;
    void GetWeightHeuristicConstrFor(
        uint32_t start,
        const std::vector<double>& weights,
        const std::function<void(std::vector<uint32_t>&&)>& callback
    ) const;

private:
    std::vector<std::vector<bool>> m_graph;
    std::vector<std::set<uint32_t>> m_adj;
    std::vector<std::set<uint32_t>> m_non_adj;

    std::vector<std::vector<uint32_t>> m_random_metrics;
    double density = 0.0;
};
