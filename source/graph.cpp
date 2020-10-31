#include <sstream>
#include <fstream>
#include <algorithm>
#include <random>
#include <numeric>
#include "graph.h"
#include "common.hpp"

Graph::Graph(const std::string& path)
{
    std::ifstream infile(path);
    if (!infile)
        throw std::runtime_error("Invalid file path!");

    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        char type = 0;
        if (!(iss >> type))
            continue;

        if (type == 'c')
        {
            continue;
        }
        else if (type == 'p')
        {
            std::string name;
            uint32_t vertex_count = 0u;
            uint32_t edge_count   = 0u;
            if (!(iss >> name >> vertex_count >> edge_count))
                continue;

            m_graph.resize(vertex_count, std::vector<bool>(vertex_count, false));
            m_adj_vec.resize(vertex_count);
            m_adj.resize(vertex_count);
            m_random_metrics.resize(vertex_count);
        }
        else if (type == 'e')
        {
            uint32_t a = 0u;
            uint32_t b = 0u;
            if (!(iss >> a >> b))
                continue;

            if (!a || !b)
                throw std::runtime_error("Unexpected vertex index!");

            m_graph[a - 1u][b - 1u] = true;
            m_graph[b - 1u][a - 1u] = true;

            m_adj[a - 1u].emplace(b - 1u);
            m_adj[b - 1u].emplace(a - 1u);

            m_adj_vec[a - 1u].emplace_back(b - 1u);
            m_adj_vec[b - 1u].emplace_back(a - 1u);
        }
        else
        {
            throw std::runtime_error("Invalid file format!");
        }
    }

    m_non_adj.resize(m_graph.size());
    for (uint32_t i = 0; i < m_graph.size(); ++i)
    {
        for (uint32_t j = 0; j < m_graph.size(); ++j)
        {
            if (i == j || m_graph[i][j])
                continue;

            m_non_adj[i].emplace(j);
            m_non_adj[j].emplace(i);
        }
    }

    for (auto& rm : m_random_metrics)
        rm = GetOrderedNodes(ColorizationType::random);
}

Graph::ColorToVerts Graph::Colorize(ColorizationType type) const
{
    auto n = m_graph.size();
    constexpr uint32_t c_no_color = std::numeric_limits<uint32_t>::max();

    ColorToVerts result;
    std::vector<bool> available(n, false);

    auto nodes = GetOrderedNodes(type);

    std::map<uint32_t, uint32_t> colors = { { nodes.front(), 0 } };

    for (const auto& node : nodes)
    {
        for (const auto& adj : m_adj[node])
            if (colors.count(adj))
                available[colors[adj]] = true;

        uint32_t cr = 0;
        for (; cr < n; cr++)
            if (available[cr] == false)
                break;

        colors[node] = cr;
        result[cr].emplace_back(node);

        for (const auto& adj : m_adj[node])
            if (colors.count(adj))
                available[colors[adj]] = false;
    }

    return result;
}

std::vector<uint32_t> Graph::GetOrderedNodes(ColorizationType type) const
{
    std::vector<uint32_t> nodes(m_graph.size(), 0);
    std::iota(nodes.begin(), nodes.end(), 0);

    std::vector<uint32_t> random_metric(m_graph.size(), 0);
    std::iota(random_metric.begin(), random_metric.end(), 0);
    std::mt19937 g(42);
    std::shuffle(random_metric.begin(), random_metric.end(), g);

    if (type == ColorizationType::maxdegree)
    {
        std::sort(nodes.begin(), nodes.end(), [&](const auto& l, const auto& r) {
            return GetDegree(l) > GetDegree(r);
        });
    }
    else if (type == ColorizationType::mindegree)
    {
        std::sort(nodes.begin(), nodes.end(), [&](const auto& l, const auto& r) {
            return GetDegree(l) < GetDegree(r);
        });
    }
    else if (type == ColorizationType::random)
    {
        std::shuffle(nodes.begin(), nodes.end(), g);
    }
    else if (type == ColorizationType::maxdegree_random)
    {
        std::sort(nodes.begin(), nodes.end(), [&](const auto& l, const auto& r) -> bool {
            if (GetDegree(l) == GetDegree(r))
                return random_metric[l] > random_metric[r];
            return GetDegree(l) > GetDegree(r);
        });
    }
    else if (type == ColorizationType::mindegree_random)
    {
        std::sort(nodes.begin(), nodes.end(), [&](const auto& l, const auto& r) -> bool {
            if (GetDegree(l) == GetDegree(r))
                return random_metric[l] > random_metric[r];
            return  GetDegree(l) < GetDegree(r);
        });
    }

    return nodes;
}

struct LocalSearchHelper
{
    const std::vector<std::set<WeightNode>>& non_adj_sorted;
    std::mt19937 g{42};

    void Search(const std::vector<uint32_t>& constr, const std::vector<double>& weights, uint32_t n, const std::function<void(std::vector<uint32_t>&&)>& callback)
    {
        constexpr uint32_t k = 3;
        if (constr.size() <= k * 2 + 1)
            return;

        auto nodes = constr;
        std::array<uint32_t, k> banned;
        for (uint32_t i = 0; i < k; ++i)
        {
            std::uniform_int_distribution<size_t> dist(0, nodes.size() - 1);
            auto val = dist(g);
            std::swap(nodes[val], nodes.back());
            banned[i] = nodes.back();
            nodes.pop_back();
        }

        auto t = non_adj_sorted[nodes.front()];
        double weight = weights[nodes.front()];
        for (uint32_t i = 1; i < nodes.size(); ++i)
        {
            t *= non_adj_sorted[nodes[i]];
            weight += weights[nodes[i]];
        }
        std::erase_if(t, [&banned](const auto& v) {
            return banned[0] == v.label || banned[1] == v.label || banned[2] == v.label;
        });

        if (t.size() <= k)
            return;

        for (const auto& start_node : t)
        {
            auto start = start_node.label;
            auto tt = t * non_adj_sorted[start];
            std::vector<uint32_t> new_nodes = { start };
            double add_weight = weights[start];

            while (!tt.empty())
            {
                auto node = tt.begin()->label;
                add_weight += tt.begin()->weight;
                new_nodes.emplace_back(node);

                tt *= non_adj_sorted[node];
            }

            if (EpsValue(weight + add_weight) <= 1.0)
                continue;

            if (new_nodes.size() + nodes.size() < constr.size())
                continue;

            new_nodes.insert(new_nodes.end(), nodes.begin(), nodes.end());
            callback(std::move(new_nodes));
        }
    }
};

std::vector<std::set<WeightNode>> Graph::GetWeightlyNonAdj(uint32_t start, const std::vector<double>& weights) const
{
    std::vector<std::set<WeightNode>> non_adj_sorted(m_graph.size());
    for (uint32_t i = 0; i < m_graph.size(); ++i)
    {
        for (auto j : m_non_adj[i])
        {
            non_adj_sorted[i].emplace(j, weights[j], m_random_metrics[start][j]);
            non_adj_sorted[j].emplace(i, weights[i], m_random_metrics[start][i]);
        }
    }
    return non_adj_sorted;
}

void Graph::GetWeightHeuristicConstr(
    const std::vector<double>& weights,
    const std::function<void(std::vector<uint32_t>&&)>& callback
) const
{
    std::set<WeightNode> consider;
    uint32_t max = 0;
    for (uint32_t i = 0; i < m_graph.size(); ++i)
    {
        auto w = weights[i];
        if (w < 0.1)
            continue;

        consider.emplace(i, w, m_random_metrics[i][i]);

        if (weights[max] >= w)
            continue;

        max = i;
    }
    if (consider.empty())
        return;

    auto non_adj_sorted = GetWeightlyNonAdj(max, weights);
    int k = 0;
    for (const auto& v : consider)
    {
        if (k++ == 10)
            break;

        GetWeightHeuristicConstrFor(v.label, weights, callback, non_adj_sorted, non_adj_sorted[v.label]);
    }
}

void Graph::GetWeightHeuristicConstrFor(uint32_t start, const std::vector<double>& weights, const std::function<void(std::vector<uint32_t>&&)>& callback, bool alternative) const
{
    auto non_adj = GetWeightlyNonAdj(start, weights);
    if (alternative)
    {
        std::set<WeightNode> consider;
        for (auto v : m_adj[start])
            consider.emplace(v, weights[v], m_random_metrics[start][v]);

        return GetWeightHeuristicConstrFor(start, weights, callback, non_adj, consider);
    }
    GetWeightHeuristicConstrFor(start, weights, callback, non_adj, non_adj[start]);
}

void Graph::GetWeightHeuristicConstrFor(
    uint32_t start,
    const std::vector<double>& weights,
    const std::function<void(std::vector<uint32_t>&&)>& callback,
    const std::vector<std::set<WeightNode>>& non_adj_sorted,
    const std::set<WeightNode>& nodes_queue
) const
{
    std::set<WeightNode> nodes = nodes_queue;

    std::vector<uint32_t> constr;
    constr.reserve(m_graph.size());
    LocalSearchHelper lsh{ non_adj_sorted };
    while (!nodes.empty())
    {
        auto t = non_adj_sorted[start] * non_adj_sorted[nodes.begin()->label];
        nodes.erase(nodes.begin());

        constr.clear();
        constr.emplace_back(start);
        double weight = weights[start];
        while (!t.empty())
        {
            auto node = t.begin()->label;
            weight += t.begin()->weight;

            t *= non_adj_sorted[node];

            constr.emplace_back(node);
        }

        if (EpsValue(weight) <= 1.0)
            continue;

        lsh.Search(constr, weights, 10, callback);
        callback(std::move(constr));
    }
}
