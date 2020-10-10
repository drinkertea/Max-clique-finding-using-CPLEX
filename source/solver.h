#pragma once

#include <cstdint>
#include <vector>

#ifndef IL_STD
#define IL_STD
#endif
#include <ilcplex/ilocplex.h>

struct Graph;

struct ConstrainsGuard
{
    ConstrainsGuard(IloModel& model, const IloExtractable& constrains);
    ~ConstrainsGuard();

private:
    IloModel&            m_model;
    const IloExtractable m_constrains;
};

struct ModelData
{
    ModelData(const Graph& graph, IloNumVar::Type type);

    IloCplex        CreateSolver() const;
    ConstrainsGuard AddScopedConstrains(uint32_t variable_index, IloNum lowerBound, IloNum upperBound = IloInfinity);
    double          ExtractValue(const IloCplex& cplex, uint32_t variable_index) const;
    size_t          GetSize() const { return m_size; }
    const Graph&    GetGraph() const { return m_graph; }

private:
    const Graph& m_graph;
    IloEnv          m_env{};
    IloModel        m_model;
    IloNumVarArray  m_variables;
    size_t          m_size = 0;
};

std::vector<uint32_t> FindMaxCliqueInteger(const Graph& graph);
std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph);
