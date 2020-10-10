#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>

ModelData::ModelData(const Graph& graph, IloNumVar::Type type)
    : m_model(m_env)
    , m_variables(m_env, graph.get().num_vertices())
    , m_graph(graph)
{
    auto n = graph.get().num_vertices();

    for (uint32_t i = 0; i < n; ++i)
    {
        std::string name = "y[" + std::to_string(i + 1) + "]";
        m_variables[i] = IloNumVar(m_env, 0, 1, type, name.c_str());
    }

    std::vector<std::pair<uint32_t, uint32_t>> non_edge_pairs;
    for (uint32_t i = 0; i < n; ++i)
    {
        for (uint32_t j = i + 1; j < n; ++j)
        {
            if (graph.get().out_neighbors(i).count(j))
                continue;

            non_edge_pairs.emplace_back(i, j);
        }
    }

    // if y[i] and y[j] has not edge - y[i] + y[j] <= 1
    IloRangeArray constrains(m_env, non_edge_pairs.size());

    uint32_t k = 0u;
    for (auto [i, j] : non_edge_pairs)
    {
        IloExpr expr(m_env);
        expr += m_variables[i];
        expr += m_variables[j];
        std::string name = "constr[" + std::to_string(i) + "][" + std::to_string(j) + "]";
        constrains[k++] = IloRange(m_env, 0, expr, 1, name.c_str());
    }

    IloExpr obj_expr(m_env);
    for (uint32_t i = 0; i < n; ++i)
    {
        obj_expr += m_variables[i];
    }
    IloObjective obj(m_env, obj_expr, IloObjective::Maximize);

    m_model.add(constrains);
    m_model.add(obj);
}

IloCplex ModelData::CreateSolver() const
{
    return { m_model };
}

ConstrainsGuard ModelData::AddScopedConstrains(uint32_t variable_index, IloNum lowerBound, IloNum upperBound)
{
    IloExpr expr(m_env);
    std::string name = std::to_string(lowerBound) + " <= y[" + std::to_string(variable_index) + "] <= " + std::to_string(upperBound);
    expr += m_variables[variable_index];
    IloRange new_constrain(m_env, lowerBound, expr, upperBound, name.c_str());
    return ConstrainsGuard(m_model, m_model.add(new_constrain));
}

double ModelData::ExtractValue(const IloCplex& cplex, uint32_t variable_index) const
{
    return cplex.getValue(m_variables[variable_index]);
}

ConstrainsGuard::ConstrainsGuard(IloModel& model, const IloExtractable& constrains)
    : m_model(model)
    , m_constrains(constrains)
{
}

ConstrainsGuard::~ConstrainsGuard()
{
    m_model.remove(m_constrains);
}

std::vector<uint32_t> FindMaxCliqueInteger(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Int);

    IloCplex cplex = model.CreateSolver();
    cplex.exportModel("model.lp");

    bool solved = cplex.solve();
    std::vector<uint32_t> res;
    std::vector<double> vars;

    auto n = graph.get().num_vertices();
    auto obj_val = cplex.getObjValue();
    if (solved) {
        for (uint32_t i = 0; i < n; ++i)
        {
            auto val = model.ExtractValue(cplex, i);
            vars.push_back(val);
            if (val)
                res.emplace_back(i);
        }
    }
    else
    {
        throw std::runtime_error("Cplex error");
    }
    return res;
}

double EpsValue(double x)
{
    constexpr double eps = 1e8;
    return std::round(x * eps) / eps;
}

bool IsInteger(double x)
{
    return EpsValue(x) == std::round(x);
}

bool IsIntegerSolution(const IloCplex& cplex, const ModelData& model)
{
    for (size_t i = 0; i < model.GetSize(); ++i)
    {
        if (!IsInteger(model.ExtractValue(cplex, i)))
            return false;
    }
    return true;
}

double DiffToInteger(double x)
{
    return std::abs(std::round(x) - x);
}

void BnB(ModelData& model, uint64_t& best_obj, std::set<uint32_t>& banned, bool is_first = false, size_t iff = -1)
{
    IloCplex cplex = model.CreateSolver();
    if (!cplex.solve())
        return;

    auto upper_bound = cplex.getObjValue();

    if (EpsValue(upper_bound) <= static_cast<double>(best_obj))
        return;

    if (IsInteger(upper_bound) && IsIntegerSolution(cplex, model))
        best_obj = static_cast<uint64_t>(upper_bound);

    auto i = -1;
    if (i == size_t(-1))
        return;

    banned.emplace(i);

    double var_value = model.ExtractValue(cplex, i);

    double left_res = 0;
    double right_res = 0;
    {
        auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
        IloCplex cplex = model.CreateSolver();
        if (cplex.solve())
            left_res = cplex.getObjValue();
    }
    {
        auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
        IloCplex cplex = model.CreateSolver();
        if (cplex.solve())
            right_res = cplex.getObjValue();
    }

    double ldiff = DiffToInteger(left_res);
    double rdiff = DiffToInteger(right_res);
    if (ldiff == rdiff ? right_res > left_res : rdiff < ldiff)
    {
        {
            auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
            BnB(model, best_obj, banned);
        }
        {
            auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
            BnB(model, best_obj, banned);
        }
    }
    else
    {
        {
            auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
            BnB(model, best_obj, banned);
        }
        {
            auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
            BnB(model, best_obj, banned);
        }
    }

    banned.erase(i);
}

void BnBC(ModelData& model, const Graph& graph, uint64_t& best_obj, uint32_t n)
{
    IloCplex cplex = model.CreateSolver();
    if (!cplex.solve())
        return;

    auto upper_bound = cplex.getObjValue();
    if (EpsValue(upper_bound) <= static_cast<double>(best_obj))
        return;

    std::vector<double> values(n, 0);
    for (uint32_t i = 0; i < n; ++i)
    {
        values[i] = model.ExtractValue(cplex, i);
    }

    if (IsInteger(upper_bound) && std::all_of(values.begin(), values.end(), [](double x) { return IsInteger(x); }))
        best_obj = static_cast<uint64_t>(upper_bound);

    auto verts_by_color = graph.Colorize();
    if (verts_by_color.size() <= best_obj)
        return;

    while (!verts_by_color.empty())
    {
        auto& color_data = *verts_by_color.rbegin();
        uint32_t color_value = color_data.first;
        uint32_t vertex_index = color_data.second.back();
        if (EpsValue(upper_bound) + static_cast<double>(color_value) <= best_obj)
            return;

        double var_value = values[vertex_index];
        {
            auto guard = model.AddScopedConstrains(vertex_index, 1.0, IloInfinity);
            BnBC(model, Graph(graph.get().subgraph(graph.get().in_neighbors(vertex_index) + graph.get().out_neighbors(vertex_index))), best_obj, n);
        }

        color_data.second.pop_back();
        if (color_data.second.empty())
            verts_by_color.erase(color_value);
    }
}

std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.get().num_vertices();

    uint64_t best_obj = 30;

    std::set<uint32_t> banned;
    std::vector<ConstrainsGuard> forever_constr;

    BnBC(model, graph, best_obj, n);
    return{};

}
