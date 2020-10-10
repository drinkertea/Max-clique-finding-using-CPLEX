#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>

void AddIndependetConst(const Graph& graph, const IloNumVarArray& m_variables, IloEnv& m_env, IloModel& m_model, Graph::ColorizationType type)
{
    auto verts_by_color = graph.Colorize(type);
    uint32_t k = 0;
    for (const auto& independed_nodes : verts_by_color)
    {
        if (independed_nodes.second.size() <= 2)
            continue;

        IloExpr expr(m_env);
        for (const auto& i : independed_nodes.second)
            expr += m_variables[i];
        std::string name = "constr for color " + std::to_string(independed_nodes.first);
        m_model.add(IloRange(m_env, 0, expr, 1, name.c_str()));
    }
}

ModelData::ModelData(const Graph& graph, IloNumVar::Type type)
    : m_model(m_env)
    , m_variables(m_env, graph.get().num_vertices())
    , m_graph(graph)
    , m_size(graph.get().num_vertices())
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

    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::default);
    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::maxdegree);
    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::mindegree);

    m_model.add(constrains);
    m_model.add(obj);
}

IloCplex ModelData::CreateSolver() const
{
    return { m_model };
}

ConstrainsGuard ModelData::AddScopedConstrains(uint32_t variable_index, IloNum lowerBound, IloNum upperBound)
{
    return ConstrainsGuard(*this, variable_index, lowerBound, upperBound);
}

double ModelData::ExtractValue(const IloCplex& cplex, uint32_t variable_index) const
{
    return cplex.getValue(m_variables[variable_index]);
}

ConstrainsGuard::ConstrainsGuard(ModelData& model, uint32_t variable_index, IloNum lowerBound, IloNum upperBound)
    : m_model(model)
    , m_expr(m_model.m_env)
{
    std::string name = std::to_string(lowerBound) + " <= y[" + std::to_string(variable_index) + "] <= " + std::to_string(upperBound);
    m_expr += m_model.m_variables[variable_index];
    m_constrain = IloRange(m_model.m_env, lowerBound, m_expr, upperBound, name.c_str());
    m_constrains = m_model.m_model.add(m_constrain);
}

ConstrainsGuard::~ConstrainsGuard()
{
    m_model.m_model.remove(m_constrains);
    m_constrain.end();
    m_expr.end();
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

double DiffToInteger(double x)
{
    return std::abs(std::round(x) - x);
}

size_t ClosestToIntIndex(const std::vector<double>& vars)
{
    double max = std::numeric_limits<double>::min();
    size_t res = -1;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        auto val = vars[i];
        if (IsInteger(val))
            continue;

        if (val <= max)
            continue;

        max = val;
        res = i;
    }
    return res;
}

std::tuple<double, std::vector<double>> Solve(ModelData& model)
{
    IloCplex cplex = model.CreateSolver();
    if (!cplex.solve())
        return { 0.0, {} };

    std::vector<double> vars(model.GetSize(), 0.0);
    for (size_t i = 0; i < model.GetSize(); ++i)
        vars[i] = model.ExtractValue(cplex, i);

    auto result = std::make_tuple(cplex.getObjValue(), std::move(vars));
    cplex.end();
    return result;
}

void BnB(ModelData& model, uint64_t& best_obj)
{
    auto [upper_bound, vars] = Solve(model);

    if (EpsValue(upper_bound) <= static_cast<double>(best_obj))
        return;

    if (IsInteger(upper_bound))
    {
        auto int_count = std::count_if(vars.begin(), vars.end(), [](double x) { return EpsValue(x) == 1.0; });
        if (int_count > best_obj)
            best_obj = static_cast<uint64_t>(int_count);
        //best_obj = static_cast<uint64_t>(upper_bound);
    }
    auto i = ClosestToIntIndex(vars);

    if (i == size_t(-1))
        return;

    double var_value = vars[i];
    {
        auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
        BnB(model, best_obj);
    }
    {
        auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
        BnB(model, best_obj);
    }
}

std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.get().num_vertices();

    uint64_t best_obj = 8;

    std::set<uint32_t> banned;
    std::vector<ConstrainsGuard> forever_constr;

    BnB(model, best_obj);
    return{};

}
