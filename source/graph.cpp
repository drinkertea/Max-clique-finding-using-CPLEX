#include <sstream>
#include <fstream>

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
            m_degrees.resize(vertex_count, 0);
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

            m_degrees[a - 1u]++;
            m_degrees[b - 1u]++;
        }
        else
        {
            throw std::runtime_error("Invalid file format!");
        }
    }
}
