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

HeuristicConstrain Extract(const std::vector<int>& is, const std::vector<double>& weights)
{
    HeuristicConstrain res{};
    for (auto x : is)
    {
        res.nodes.emplace(x);
        res.weight += weights[x];
    }
    return res;
}

std::vector<uint32_t> Graph::GetOrderedNodes(ColorizationType type) const
{
    std::vector<uint32_t> nodes(m_graph.size(), 0);
    std::iota(nodes.begin(), nodes.end(), 0);

    std::vector<uint32_t> random_metric(m_graph.size(), 0);
    std::iota(random_metric.begin(), random_metric.end(), 0);
    static uint64_t seed = 42;
    std::mt19937 g(seed++);
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

void Graph::GetWeightHeuristicConstrFor(uint32_t start, const std::vector<double>& weights, const std::function<void(const HeuristicConstrain&)>& callback) const
{
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

    std::vector<std::set<WeightNode>> non_adj_sorted(m_graph.size());
    for (uint32_t i = 0; i < m_graph.size(); ++i)
    {
        for (auto j : m_non_adj[i])
        {
            non_adj_sorted[i].emplace(j, weights[j], m_random_metrics[start][j]);
            non_adj_sorted[j].emplace(i, weights[i], m_random_metrics[start][i]);
        }
    }
    std::set<WeightNode> nodes = non_adj_sorted[start];

    HeuristicConstrain base_constr{ weights[start], { start } };

    while (!nodes.empty())
    {
        auto t = non_adj_sorted[start] * non_adj_sorted[nodes.begin()->label];
        nodes.erase(nodes.begin());

        HeuristicConstrain constr = base_constr;
        while (!t.empty())
        {
            auto node = t.begin()->label;
            auto weight = t.begin()->weight;

            t *= non_adj_sorted[node];

            constr.weight += weight;
            constr.nodes.emplace(node);
        }

        if (EpsValue(constr.weight) <= 1.0)
            continue;

        callback(constr);
    }
}
