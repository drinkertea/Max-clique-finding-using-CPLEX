#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>

ModelData::ModelData(const Graph& graph, IloNumVar::Type type)
    : m_model(m_env)
    , m_variables(m_env, graph.Size())
    , m_size(graph.Size())
    , m_graph(graph)
{
    auto n = graph.Size();

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
            if (graph.HasEdge(i, j))
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

    auto n = graph.Size();
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

template<typename Float>
constexpr auto cround(Float x) -> decltype(std::round(x))
{
    return (x >= 0.0) ? int(x + 0.5) : int(x - 0.5);
}
//static_assert(EpsValue(0.999999999) == 1.0);
//static_assert(EpsValue(0.000000001) == 0.0);
//static_assert(IsInteger(0.000000001));
//static_assert(IsInteger(0.999999999));
//static_assert(!IsInteger(0.99999999));
//static_assert(!IsInteger(0.00000001));

// work only for positive
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
//
//size_t ClosestToIntIndex(const IloCplex& cplex, const ModelData& model, const std::set<uint32_t>& banned)
//{
//    double max = std::numeric_limits<double>::min();
//    size_t res = -1;
//    for (size_t i = 0; i < model.GetSize(); ++i)
//    {
//        auto val = model.ExtractValue(cplex, i);
//        if (EpsValue(val) == 1.0)
//            continue;
//
//        if (val <= max)
//            continue;
//
//        if (banned.count(i))
//            continue;
//
//        max = val;
//        res = i;
//    }
//    return res;
//}
//
//std::map<std::vector<double>, uint64_t> results;
//
//void BnB(ModelData& model, uint64_t& best_obj, std::set<uint32_t>& banned, std::vector<ConstrainsGuard>& forever_constr)
//{
//    IloCplex cplex = model.CreateSolver();
//    if (!cplex.solve())
//        return;// throw std::runtime_error("Unsolved model!");
//
//    auto upper_bound = cplex.getObjValue();
//
//    std::vector<double> vars_vals;
//    vars_vals.reserve(model.GetSize());
//    for (size_t i = 0; i < model.GetSize(); ++i)
//    {
//        vars_vals.push_back(model.ExtractValue(cplex, i));
//        //std::cout << vars_vals.back() << " ";
//    }
//    //std::cout << std::endl;
//
//    if (++results[vars_vals] > 2)
//    {
//        volatile int k = 0;
//        ++k; --k;
//    }
//
//    if (EpsValue(upper_bound) <= static_cast<double>(best_obj))
//        return;
//
//    if (IsInteger(upper_bound) && IsIntegerSolution(cplex, model))
//    {
//        best_obj = static_cast<uint64_t>(upper_bound);
//        //for (size_t i = 0; i < model.GetSize(); ++i)
//        //{
//        //    auto val = model.ExtractValue(cplex, i);
//        //    if (EpsValue(val) == 1.0)
//        //        continue;
//        //
//        //    if (model.GetGraph().GetDegree(i) >= best_obj)
//        //        continue;
//        //
//        //    if (banned.count(i))
//        //        continue;
//        //
//        //    banned.emplace(i);
//        //    forever_constr.push_back(model.AddScopedConstrains(i, 0.0, 0.0));
//        //}
//    }
//
//    auto i = ClosestToIntIndex(cplex, model, banned);
//    if (i == size_t(-1))
//        return;
//
//    double var_value = model.ExtractValue(cplex, i);
//
//    double left_res = 0;
//    double right_res = 0;
//    {
//        auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
//        IloCplex cplex = model.CreateSolver();
//        if (cplex.solve())
//            left_res = cplex.getObjValue();
//    }
//    {
//        auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
//        IloCplex cplex = model.CreateSolver();
//        if (cplex.solve())
//            right_res = cplex.getObjValue();
//    }
//
//    double ldiff = DiffToInteger(left_res);
//    double rdiff = DiffToInteger(right_res);
//    if (ldiff == rdiff ? right_res > left_res : rdiff < ldiff)
//    {
//        {
//            auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
//            BnB(model, best_obj, banned, forever_constr);
//        }
//        {
//            auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
//            BnB(model, best_obj, banned, forever_constr);
//        }
//    }
//    else
//    {
//        {
//            auto guard = model.AddScopedConstrains(i, 0.0, std::floor(var_value));
//            BnB(model, best_obj, banned, forever_constr);
//        }
//        {
//            auto guard = model.AddScopedConstrains(i, std::ceil(var_value), IloInfinity);
//            BnB(model, best_obj, banned, forever_constr);
//        }
//    }
//}

size_t BranchingIndex(const IloCplex& cplex, const ModelData& model, const std::set<uint32_t>& banned)
{
    uint64_t max = 0;
    size_t res = -1;
    for (size_t i = 0; i < model.GetSize(); ++i)
    {
        auto val = model.GetGraph().GetDegree(i);    
        if (val <= max)
            continue;
    
        if (banned.count(i))
            continue;
    
        max = val;
        res = i;
    }
    return res;
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

    auto i = is_first ? iff : BranchingIndex(cplex, model, banned);
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


std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.Size();

    uint64_t best_obj = 1;

    std::set<uint32_t> banned;
    std::vector<ConstrainsGuard> forever_constr;

    for (size_t i = 0; i < model.GetSize(); ++i)
        BnB(model, best_obj, banned, true, i);

    return{};

}
