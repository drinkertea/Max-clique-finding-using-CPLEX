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
    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::random);
    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::maxdegree_random);
    AddIndependetConst(graph, m_variables, m_env, m_model, Graph::ColorizationType::mindegree_random);

    m_model.add(constrains);
    m_model.add(obj);
}

IloCplex ModelData::CreateSolver() const
{
    IloCplex res(m_model);
    res.setOut(m_env.getNullStream());
    return res;
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

uint64_t bc = 0;

struct BnBhelper
{
    ModelData& model;

    BnBhelper(ModelData& m) : model(m) {};

    uint64_t best_obj = 1;

    struct Solution
    {
        double upper_bound = 0.0;
        std::vector<double> vars;
        size_t max_non_int_index = -1;
        uint64_t int_count = 0;
    };

    Solution Solve()
    {
        Solution res{};

        IloCplex cplex = model.CreateSolver();
        if (!cplex.solve())
            return res;

        res.vars.resize(model.GetSize(), 0.0);
        double max = std::numeric_limits<double>::min();
        for (size_t i = 0; i < model.GetSize(); ++i)
        {
            double val = model.ExtractValue(cplex, i);
            res.vars[i] = val;
            if (IsInteger(val))
            {
                res.int_count += uint64_t(EpsValue(val) == 1.0);
                continue;
            }

            if (val <= max)
                continue;

            max = val;
            res.max_non_int_index = i;
        }

        res.upper_bound = cplex.getObjValue();
        cplex.end();
        return res;
    }

    size_t SelectBranch(const std::vector<double>& vars)
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

    void BnB(const Solution& solution)
    {
        bc++;

        if (EpsValue(solution.upper_bound) <= static_cast<double>(best_obj))
            return;

        if (IsInteger(solution.upper_bound))
        {
            if (solution.int_count >= best_obj)
            {
                best_obj = static_cast<uint64_t>(solution.int_count);
                std::cout << "Found " << best_obj << " " << bc << std::endl;
                for (uint64_t i = 0; i < solution.vars.size(); ++i)
                    if (solution.vars[i] > 0.999999)
                        std::cout << i << " ";
                std::cout << std::endl;
                std::cout << std::endl;
            }
        }

        size_t i = solution.max_non_int_index;
        if (i == size_t(-1))
            return;

        Solution right_solution{};
        {
            auto guard = model.AddScopedConstrains(i, 1.0, 1.0);
            right_solution = Solve();
        }
        Solution left_solution{};
        {
            auto guard = model.AddScopedConstrains(i, 0.0, 0.0);
            left_solution = Solve();
        }

        double rdiff = DiffToInteger(right_solution.upper_bound);
        double ldiff = DiffToInteger(left_solution.upper_bound);

        bool go_right = rdiff < ldiff;
        if (EpsValue(rdiff) == EpsValue(ldiff))
            go_right = right_solution.upper_bound > left_solution.upper_bound;

        if (go_right)
        {
            {
                auto guard = model.AddScopedConstrains(i, 1.0, 1.0);
                BnB(right_solution);
            }
            {
                auto guard = model.AddScopedConstrains(i, 0.0, 0.0);
                BnB(left_solution);
            }
        }
        else
        {
            {
                auto guard = model.AddScopedConstrains(i, 0.0, 0.0);
                BnB(left_solution);
            }
            {
                auto guard = model.AddScopedConstrains(i, 1.0, 1.0);
                BnB(right_solution);
            }
        }
    }
};


std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.get().num_vertices();

    uint64_t best_obj = 8;

    std::set<uint32_t> banned;
    std::vector<ConstrainsGuard> forever_constr;

    BnBhelper bnbh(model);
    auto initial_solution = bnbh.Solve();
    bnbh.BnB(initial_solution);
    return{};

}
