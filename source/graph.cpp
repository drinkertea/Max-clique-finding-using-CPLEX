#include <sstream>
#include <fstream>

#include "graph.h"

Graph::Graph(const std::string& path)
{
    m_graph.set_undirected();

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
        }
        else if (type == 'e')
        {
            uint32_t a = 0u;
            uint32_t b = 0u;
            if (!(iss >> a >> b))
                continue;

            if (!a || !b)
                throw std::runtime_error("Unexpected vertex index!");

            m_graph.insert_edge(a - 1u, b - 1u);
        }
        else
        {
            throw std::runtime_error("Invalid file format!");
        }
    }

    Colorize();
}

Graph::ColorToVerts Graph::Colorize() const
{
    auto n = m_graph.rbegin()->first + 1;
    constexpr uint32_t c_no_color = std::numeric_limits<uint32_t>::max();

    ColorToVerts result;
    std::map<uint32_t, uint32_t> colors = { { m_graph.begin()->first, 0 } };
    result[0].emplace_back(m_graph.begin()->first);

    std::vector<bool> available(n, false);

    for (const auto& node : m_graph)
    {
        auto all = node.second.second + node.second.first;
        for (const auto& adj : all)
            if (colors.count(adj))
                available[colors[adj]] = true;

        uint32_t cr = 0;
        for (; cr < n; cr++)
            if (available[cr] == false)
                break;

        colors[node.first] = cr;
        result[cr].emplace_back(node.first);

        for (const auto& adj : all)
            if (colors.count(adj))
                available[colors[adj]] = false;
    }

    return result;
}
