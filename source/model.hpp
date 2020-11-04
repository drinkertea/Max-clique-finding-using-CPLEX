#pragma once
#include "common.hpp"
#include "utils.hpp"

struct ModelData;
struct ConstrainsGuard
{
    ConstrainsGuard(ModelData& model, size_t variable_index, IloNum lowerBound, IloNum upperBound);
    ConstrainsGuard(ModelData& model, const std::vector<uint32_t>& constr);
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
        , m_size(graph.GetSize())
        , m_type(type)
        , m_model(m_env)
        , m_variables(m_env, m_size)
    {
        if (type == IloNumVar::Type::Float)
        {
            for (auto type : g_strategies)
                AddPerColorConstrains(type);

            if (!AddRecursiveConstrains())
            {
                for (auto type : g_strategies)
                    m_sum_vert_less_one += m_graph.GetHeuristicConstr(type);
            }
        }
        AddNonEdgePairs();

        InitModel(m_sum_vert_less_one);
    }

    ModelData(const ModelData& r)
        : m_graph(r.m_graph)
        , m_size(r.m_size)
        , m_type(r.m_type)
        , m_model(r.m_model.getClone(m_env))
        , m_variables(m_env, m_size)
    {
        auto it = IloModel::Iterator(m_model);
        while (it.ok())
        {
            auto extr = *it;
            if (extr.isVariable())
                m_variables[std::stoi(extr.getName())] = extr.asVariable();
            ++it;
        }
    }

    ~ModelData()
    {
        m_env.end();
    }

    IloCplex CreateSolver() const
    {
        IloCplex res(m_model);
        res.setOut(m_env.getNullStream());
        return res;
    }

    ConstrainsGuard AddScopedConstrains(size_t variable_index, IloNum lowerBound, IloNum upperBound)
    {
        return ConstrainsGuard(*this, variable_index, lowerBound, upperBound);
    }

    std::unique_ptr<ConstrainsGuard> AddScopedConstrainsPtr(size_t variable_index, IloNum lowerBound, IloNum upperBound)
    {
        return std::make_unique<ConstrainsGuard>(*this, variable_index, lowerBound, upperBound);
    }

    std::unique_ptr<ConstrainsGuard> AddScopedConstrain(const std::vector<uint32_t>& constr)
    {
        return std::make_unique<ConstrainsGuard>(*this, constr);
    }

    double ExtractValue(const IloCplex& cplex, size_t variable_index) const
    {
        return cplex.getValue(m_variables[variable_index]);
    }

    size_t GetSize() const
    {
        return m_size;
    }

private:
    void AddConstrains(const std::set<std::set<uint32_t>>& m_sum_vert_less_one)
    {
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
        m_model.add(constrains);
    }

    void InitModel(const std::set<std::set<uint32_t>>& m_sum_vert_less_one)
    {
        for (uint32_t i = 0; i < m_size; ++i)
        {
            std::string name = std::to_string(i);
            m_variables[i] = IloNumVar(m_env, 0, 1, m_type, name.c_str());
        }

        IloExpr obj_expr(m_env);
        for (uint32_t i = 0; i < m_size; ++i)
            obj_expr += m_variables[i];
        IloObjective obj(m_env, obj_expr, IloObjective::Maximize);

        AddConstrains(m_sum_vert_less_one);
        m_model.add(obj);
        m_model.add(m_variables);
    }

    void AddPerColorConstrains(Graph::ColorizationType type)
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
                if (m_graph.HasEdge(i, j))
                    continue;

                m_sum_vert_less_one.emplace(std::set<uint32_t>{ i, j });
            }
        }
    }

    bool AddRecursiveConstrains()
    {
        // If v1 .. vk all has not edge -> v1 + v2 + ... + vk <= 1
        uint32_t      break_count = uint32_t(double(m_size) * std::exp(3.14159));
        Timer         timer;
        TimeoutThread break_timer(std::chrono::milliseconds(500));

        // searching for optimal "k" < "m_break_count" to no add too much constrains
        NonEdgeKHelper helper{ m_size, break_count, m_size, m_graph, break_timer.stopped };

        helper.AddNonEdgeRec();

        if (break_timer.stopped)
        {
            std::cout << "AddHeuristicConstrains interrupted" << std::endl;
            return false;
        }

        for (uint32_t i = 3; i < m_size; ++i)
        {
            if (helper.res[i].size() >= break_count)
                continue;

            for (const auto& constr : helper.res[i])
                m_sum_vert_less_one.emplace(std::set<uint32_t>{constr.begin(), constr.end()});
        }

        std::cout << "AddHeuristicConstrains " << timer.Stop() << " ms" << std::endl;
        return true;
    }

    friend struct ConstrainsGuard;

    std::set<std::set<uint32_t>> m_sum_vert_less_one;

    const Graph&    m_graph;
    size_t          m_size = 0;
    IloNumVar::Type m_type;
    IloEnv          m_env{};
    IloModel        m_model;
    IloNumVarArray  m_variables;
};

ConstrainsGuard::ConstrainsGuard(ModelData& model, size_t variable_index, IloNum lowerBound, IloNum upperBound)
    : m_model(model)
    , m_expr(m_model.m_env)
{
    std::string name = std::to_string(lowerBound) + " <= y[" + std::to_string(variable_index) + "] <= " + std::to_string(upperBound);
    m_expr += m_model.m_variables[variable_index];
    m_constrain = IloRange(m_model.m_env, lowerBound, m_expr, upperBound, name.c_str());
    m_constrains = m_model.m_model.add(m_constrain);
}

ConstrainsGuard::ConstrainsGuard(ModelData& model, const std::vector<uint32_t>& constr)
    : m_model(model)
    , m_expr(m_model.m_env)
{
    for (auto index : constr)
        m_expr += m_model.m_variables[index];
    m_constrain = IloRange(m_model.m_env, 0.0, m_expr, 1.0);
    m_constrains = m_model.m_model.add(m_constrain);
}

ConstrainsGuard::~ConstrainsGuard()
{
    m_model.m_model.remove(m_constrains);
    m_constrain.end();
    m_expr.end();
}
