#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>

std::vector<uint32_t> FindMaxClique(const Graph& graph)
{
    IloEnv env{};
    IloModel model(env);

    auto n = graph.Size();

    // y[i] == 1 if y[1] is a max clique vertex
    IloNumVarArray y(env, n);

    for (uint32_t i = 0; i < n; ++i)
    {
        std::string name = "y[" + std::to_string(i + 1) + "]";
        y[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.c_str());
    }

    std::vector<std::pair<uint32_t, uint32_t>> non_edge_pairs;
    for (uint32_t i = 0; i < n; ++i)
    {
        for (uint32_t j = i + 1; j < n; ++j)
        {
            if (graph.HasEdge(i, j))
                continue;

            non_edge_pairs.emplace_back(i, j);
        }
    }

    // if y[i] and y[j] has not edge - y[i] + y[j] <= 1
    IloRangeArray constrains(env, non_edge_pairs.size());

    uint32_t k = 0u;
    for (auto [i, j] : non_edge_pairs)
    {
        IloExpr expr(env);
        expr += y[i];
        expr += y[j];
        std::string name = "constr[" + std::to_string(i) + "][" + std::to_string(j) + "]";
        constrains[k++] = IloRange(env, 0, expr, 1, name.c_str());
    }

    IloExpr obj_expr(env);
    for (uint32_t i = 0; i < n; ++i)
    {
        obj_expr += y[i];
    }
    IloObjective obj(env, obj_expr, IloObjective::Maximize);

    model.add(constrains);
    model.add(obj);

    IloCplex cplex(model);
    cplex.exportModel("model.lp");

    bool solved = cplex.solve();
    std::vector<uint32_t> res;
    if (solved) {
        for (uint32_t i = 0; i < n; ++i)
        {
            if (cplex.getValue(y[i]) > 0.)
                res.emplace_back(i);
        }
    }
    else
    {
        throw std::runtime_error("Cplex error");
    }
    return res;
}
