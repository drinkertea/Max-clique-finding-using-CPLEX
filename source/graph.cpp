#include <sstream>
#include <fstream>
#include <algorithm>
#include <random>
#include <numeric>
#include "graph.h"

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
            m_adj.resize(vertex_count);
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
        }
        else
        {
            throw std::runtime_error("Invalid file format!");
        }
    }
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
