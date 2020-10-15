#pragma once
#include "common.hpp"

bool Heuristic(const Graph& first_graph, const Graph& graph, std::set<uint32_t>& curr_clique, Graph::ColorizationType strategy)
{
    auto color_data = graph.Colorize(strategy);
    auto higher_color = *color_data.rbegin();
    if (color_data.size() == 1)
        return false;

    auto min_degree_vert = higher_color.second[0];
    auto min_degree = graph.GetDegree(higher_color.second[0]);
    for (auto vert : higher_color.second)
    {
        auto degree = graph.GetDegree(vert);
        if (min_degree <= degree)
            continue;

        min_degree = degree;
        min_degree_vert = vert;
    }

    curr_clique.emplace(min_degree_vert);
    if (!Heuristic(first_graph, graph.GetSubGraph(min_degree_vert), curr_clique, strategy))
    {
        for (auto neigh : graph.GetNeighbors(min_degree_vert))
        {
            curr_clique.emplace(neigh);
            if (first_graph.IsClique(curr_clique))
                return true;

            curr_clique.erase(neigh);
        }
    }
    return true;
}

std::set<uint32_t> Heuristic(const Graph& graph)
{
    std::set<uint32_t> best_clique;
    for (auto str : g_strategies)
    {
        std::set<uint32_t> curr_clique;
        Heuristic(graph, graph, curr_clique, str);
        if (curr_clique.size() > best_clique.size())
            best_clique = std::move(curr_clique);
    }
    return best_clique;
}
