#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <array>


struct ModelData;
struct ConstrainsGuard
{
    ConstrainsGuard(ModelData& model, uint32_t variable_index, IloNum lowerBound, IloNum upperBound);
    ~ConstrainsGuard();

private:
    ModelData& m_model;
    IloExpr m_expr;
    IloRange m_constrain;
    IloExtractable m_constrains;
};

struct ModelData
{
    ModelData(const Graph& graph, IloNumVar::Type type)
        : m_graph(graph)
        , m_size(graph.get().num_vertices())
        , m_type(type)
        , m_model(m_env)
        , m_variables(m_env, m_size)
    {
        AddNonEdgePairs();
        AddNonEdgeTrios();
        AddIndependetConst(Graph::ColorizationType::default);
        AddIndependetConst(Graph::ColorizationType::maxdegree);
        AddIndependetConst(Graph::ColorizationType::mindegree);
        AddIndependetConst(Graph::ColorizationType::random);
        AddIndependetConst(Graph::ColorizationType::maxdegree_random);
        AddIndependetConst(Graph::ColorizationType::mindegree_random);
        InitModel();
    }

    ModelData(const ModelData& r)
        : m_graph(r.m_graph)
        , m_size(r.m_size)
        , m_type(r.m_type)
        , m_model(m_env)
        , m_variables(m_env, m_size)
        , m_sum_vert_less_one(r.m_sum_vert_less_one)
    {
        InitModel();
    }

    void InitModel()
    {
        for (uint32_t i = 0; i < m_size; ++i)
        {
            std::string name = "y[" + std::to_string(i + 1) + "]";
            m_variables[i] = IloNumVar(m_env, 0, 1, m_type, name.c_str());
        }

        IloRangeArray constrains(m_env, m_sum_vert_less_one.size());
        uint32_t k = 0u;
        for (const auto& constr : m_sum_vert_less_one)
        {
            IloExpr expr(m_env);
            for (auto vert : constr)
                expr += m_variables[vert];
            std::string name = "constr[" + std::to_string(k) + "]";
            constrains[k++] = IloRange(m_env, 0, expr, 1, name.c_str());
        }

        IloExpr obj_expr(m_env);
        for (uint32_t i = 0; i < m_size; ++i)
        {
            obj_expr += m_variables[i];
        }
        IloObjective obj(m_env, obj_expr, IloObjective::Maximize);

        m_model.add(constrains);
        m_model.add(obj);
    }

    IloCplex CreateSolver() const
    {
        IloCplex res(m_model);
        res.setOut(m_env.getNullStream());
        return res;
    }

    ConstrainsGuard AddScopedConstrains(uint32_t variable_index, IloNum lowerBound, IloNum upperBound)
    {
        return ConstrainsGuard(*this, variable_index, lowerBound, upperBound);
    }

    double ExtractValue(const IloCplex& cplex, uint32_t variable_index) const
    {
        return cplex.getValue(m_variables[variable_index]);
    }

    size_t GetSize() const
    { 
        return m_size; 
    }

private:
    void AddIndependetConst(Graph::ColorizationType type)
    {
        auto verts_by_color = m_graph.Colorize(type);
        uint32_t k = 0;
        for (const auto& independed_nodes : verts_by_color)
        {
            if (independed_nodes.second.size() <= 2)
                continue;

            m_sum_vert_less_one.emplace(std::set<uint32_t>{
                independed_nodes.second.begin(),
                    independed_nodes.second.end()
            });
        }
    }

    void AddNonEdgePairs()
    {
        for (uint32_t i = 0; i < m_size; ++i)
        {
            for (uint32_t j = i + 1; j < m_size; ++j)
            {
                if (m_graph.get().out_neighbors(i).count(j))
                    continue;

                m_sum_vert_less_one.emplace(std::set<uint32_t>{ i, j });
            }
        }
    }

    void AddNonEdgeTrios()
    {
        for (uint32_t i = 0; i < m_size; ++i)
        {
            for (uint32_t j = i + 1; j < m_size; ++j)
            {
                for (uint32_t k = j + 1; k < m_size; ++k)
                {
                    if (m_graph.get().out_neighbors(i).count(j))
                        continue;

                    if (m_graph.get().out_neighbors(i).count(k))
                        continue;

                    if (m_graph.get().out_neighbors(j).count(k))
                        continue;

                    m_sum_vert_less_one.emplace(std::set<uint32_t>{ i, j, k });
                }
            }
        }
    }

    friend struct ConstrainsGuard;

    std::set<std::set<uint32_t>> m_sum_vert_less_one;

    const Graph&    m_graph;
    IloEnv          m_env{};
    IloModel        m_model;
    IloNumVarArray  m_variables;
    size_t          m_size = 0;
    IloNumVar::Type m_type;
};

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

struct BnBhelper
{
    ModelData& model;

    uint64_t bc = 0;
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

void Heuristic(const Graph& graph, std::set<uint32_t>& curr_clique)
{
    if (graph.get().num_vertices() == 0)
        return;

    auto color_data = graph.Colorize(Graph::ColorizationType::mindegree_random);
    auto higher_color = *color_data.rbegin();

    static const auto get_degree = [](const Graph& graph, uint32_t vert) {
        return graph.get().in_neighbors(vert).size() + graph.get().out_neighbors(vert).size();
    };
    auto min_degree_vert = higher_color.second[0];
    auto min_degree = get_degree(graph, higher_color.second[0]);

    for (auto vert : higher_color.second)
    {
        auto degree = get_degree(graph, vert);
        if (min_degree <= degree)
            continue;

        min_degree = degree;
        min_degree_vert = vert;
    }

    curr_clique.emplace(min_degree_vert);
    Heuristic(Graph(graph.get().subgraph(graph.get().in_neighbors(min_degree_vert) + graph.get().out_neighbors(min_degree_vert))), curr_clique);
}

std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.get().num_vertices();

    std::set<uint32_t> curr_clique;
    Heuristic(graph, curr_clique);

    BnBhelper bnbh(model);
    auto initial_solution = bnbh.Solve();
    std::cout << "Initial solution: " << initial_solution.upper_bound << std::endl << std::endl;

    for (auto vert : curr_clique)
    {
        initial_solution.vars[vert] = 1.0;
    }
    initial_solution.max_non_int_index = bnbh.SelectBranch(initial_solution.vars);
    initial_solution.upper_bound = static_cast<double>(curr_clique.size());
    initial_solution.int_count = curr_clique.size();
    bnbh.BnB(initial_solution);

    std::cout << "Total branches: " << bnbh.bc << std::endl << std::endl;

    return{};

}
