#include "solver.h"
#include "graph.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <array>
#include <deque>
#include <functional>
#include <mutex>
#include <atomic>

struct ModelData;
struct ConstrainsGuard
{
    ConstrainsGuard(ModelData& model, size_t variable_index, IloNum lowerBound, IloNum upperBound);
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
        //AddNonEdgeTrios();
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

    ~ModelData()
    {
        m_env.end();
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

    ConstrainsGuard AddScopedConstrains(size_t variable_index, IloNum lowerBound, IloNum upperBound)
    {
        return ConstrainsGuard(*this, variable_index, lowerBound, upperBound);
    }

    std::unique_ptr<ConstrainsGuard> AddScopedConstrainsPtr(size_t variable_index, IloNum lowerBound, IloNum upperBound)
    {
        return std::make_unique<ConstrainsGuard>(*this, variable_index, lowerBound, upperBound);
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

constexpr size_t g_invalid_index = std::numeric_limits<size_t>::max();

struct Solution
{
    double upper_bound = 0.0;
    std::vector<double> vars;
    size_t max_non_int_index = g_invalid_index;
    uint64_t int_count = 0;
};

using Path = std::vector<std::pair<size_t, double>>;

struct DepthGuard
{
    bool stop = false;

    DepthGuard(Path& p, bool s)
        : path(p)
        , stop(s)
    {
    }

    ~DepthGuard()
    {
        path.pop_back();
    }

private:
    Path& path;
};

struct BranchingCallback
{
    BranchingCallback(size_t md)
        : max_depth(md)
    {
    }

    DepthGuard OnBranch(const Solution& solution, size_t vertex, double constr)
    {
        curr_path.emplace_back(vertex, constr);
        if (curr_path.size() > max_depth)
        {
            branches.emplace(curr_path, solution);
            return DepthGuard(curr_path, true);
        }
        return DepthGuard(curr_path, false);
    }

    const std::map<Path, Solution>& GetBranches()
    {
        return branches;
    }

private:
    std::map<Path, Solution> branches;
    size_t                   max_depth = 0;
    Path                     curr_path;
};

class BnBHelper;

struct Printer
{
    Printer(const Graph& g)
        : graph(g)
    {
    }

    void AddCallback(const std::shared_ptr<BnBHelper>& callback)
    {
        std::lock_guard<std::mutex> lg(print_mutex);
        best_upd_callbacks.emplace_back(callback);
    }

    uint64_t GetBest() const
    {
        return best_obj;
    }

    void PrintAndValidate(const Solution& solution, uint64_t bc);

private:
    uint64_t best_obj = 0;
    std::vector<std::weak_ptr<BnBHelper>> best_upd_callbacks;
    std::mutex print_mutex;
    const Graph& graph;
};

class BnBHelper
{
    uint64_t best_obj = 1;
    uint64_t bc = 0;

    ModelData&         model;
    Printer&           printer;
    BranchingCallback* branch_callback = nullptr;

    struct BranchingData
    {
        using CRef = std::reference_wrapper<const BranchingData>;

        double way = 0.0;
        Solution res;
    };

public:
    BnBHelper(ModelData& m, Printer& p, uint64_t best, BranchingCallback* callback = nullptr)
        : model(m)
        , printer(p)
        , branch_callback(callback)
        , best_obj(best)
    {
    };

    uint64_t GetBranchCount() const
    {
        return bc;
    }

    void NewBestObj(uint64_t newval)
    {
        best_obj = std::max(newval, best_obj);
    }

    Solution Solve()
    {
        Solution res{};

        IloCplex cplex = model.CreateSolver();
        cplex.setParam(IloCplex::Param::Threads, 1);
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

    size_t SelectBranch(const std::vector<double>& vars) const
    {
        double max = std::numeric_limits<double>::min();
        size_t res = g_invalid_index;

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

    void Branching(const Solution& solution, size_t vertex, double constr)
    {
        auto guard = model.AddScopedConstrains(vertex, constr, constr);
        BnB(solution);
    }

    void Branching(BranchingCallback& callback, const Solution& solution, size_t vertex, double constr)
    {
        auto guard = callback.OnBranch(solution, vertex, constr);
        if (guard.stop)
            return;
        Branching(solution, vertex, constr);
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
                printer.PrintAndValidate(solution, bc);
            }
        }

        size_t i = solution.max_non_int_index;
        if (i == g_invalid_index)
            return;

        BranchingData right{ 1.0 };
        {
            auto guard = model.AddScopedConstrains(i, right.way, right.way);
            right.res = Solve();
        }

        BranchingData left{ 0.0 };
        {
            auto guard = model.AddScopedConstrains(i, left.way, left.way);
            left.res = Solve();
        }

        bool go_right = right.res.int_count > left.res.int_count;
        if (right.res.int_count == left.res.int_count)
            go_right = right.res.upper_bound > left.res.upper_bound;

        BranchingData::CRef right_ref = right;
        BranchingData::CRef left_ref = left;
        if (!go_right)
            std::swap(right_ref, left_ref);

        if (branch_callback)
        {
            Branching(*branch_callback, right_ref.get().res, i, right_ref.get().way);
            Branching(*branch_callback, left_ref.get().res, i, left_ref.get().way);
        }
        else
        {
            Branching(right_ref.get().res, i, right_ref.get().way);
            Branching(left_ref.get().res, i, left_ref.get().way);
        }
    }
};

void Printer::PrintAndValidate(const Solution& solution, uint64_t bc)
{
    std::lock_guard<std::mutex> lg(print_mutex);
    if (solution.int_count > best_obj)
    {
        best_obj = solution.int_count;
        for (auto& callback : best_upd_callbacks)
        {
            auto catch_obj = callback.lock();
            if (!catch_obj)
                continue;

            catch_obj->NewBestObj(best_obj);
        }
    }

    std::cout << "Found " << solution.int_count << " " << bc << std::endl;

    std::set<uint32_t> clique;
    for (uint64_t i = 0; i < solution.vars.size(); ++i)
        if (EpsValue(solution.vars[i]) == 1.0)
        {
            std::cout << i << " ";
            clique.emplace(i);
        }

    bool valid = true;
    for (auto v : clique)
    {
        auto neights = graph.get().in_neighbors(v) + graph.get().out_neighbors(v);
        for (auto v1 : clique)
        {
            if (v != v1 && !neights.count(v1))
                valid = false;
        }
    }
    std::cout << std::endl;
    std::cout << "Valid " << (valid ? "True" : "False") << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
}

struct ParallelExecturor
{
    using Task = std::function<void()>;
    using Tasks = std::deque<Task>;
    ParallelExecturor(Tasks&& tasks, uint32_t thread_count = 0)
        : tasks(std::move(tasks))
    {
        uint32_t n = thread_count;
        if (n <= 0)
            n = std::max(std::thread::hardware_concurrency(), 1u);
        for (uint32_t i = 0; i < n; ++i)
            execution_threads.emplace_back(&ParallelExecturor::ExecutionThread, this);
    }

    void WaitForJobDone()
    {
        for (auto& thred : execution_threads)
            thred.join();
    }

private:
    void ExecutionThread()
    {
        while (!tasks.empty())
        {
            Task task;
            {
                std::unique_lock<std::mutex> lock(tasks_mutex);
                task.swap(tasks.front());
                tasks.pop_front();
            }
            task();
        }
    }

    std::deque<Task>         tasks;
    std::mutex               tasks_mutex;
    std::vector<std::thread> execution_threads;
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

void AddHeuristicSolution(const BnBHelper& helper, const Graph& graph, Solution& initial_solution)
{
    std::set<uint32_t> curr_clique;
    Heuristic(graph, curr_clique);
    for (auto vert : curr_clique)
        initial_solution.vars[vert] = 1.0;

    initial_solution.max_non_int_index = helper.SelectBranch(initial_solution.vars);
    initial_solution.upper_bound = static_cast<double>(curr_clique.size());
    initial_solution.int_count = curr_clique.size();
}

std::vector<uint32_t> FindMaxCliqueBnB(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Float);
    auto n = graph.get().num_vertices();

    Printer printer(graph);
    BranchingCallback first_branching(7);

    BnBHelper bnbh(model, printer, 1, &first_branching);
    auto initial_solution = bnbh.Solve();
    std::cout << "Initial solution: " << initial_solution.upper_bound << std::endl << std::endl;

    AddHeuristicSolution(bnbh, graph, initial_solution);

    bnbh.BnB(initial_solution);

    ParallelExecturor::Tasks tasks;
    std::atomic<size_t> finished_index = 0;
    std::atomic<size_t> tasks_total = 0;
    for (const auto& first_branches : first_branching.GetBranches())
    {
        const auto& path = first_branches.first;
        const auto& solution = first_branches.second;
        tasks.push_back([path, solution, &model, &printer, &finished_index, &tasks_total]() {
            ModelData copy_model(model);

            std::vector<std::unique_ptr<ConstrainsGuard>> path_constr;
            for (const auto& branch : path)
                path_constr.push_back(copy_model.AddScopedConstrainsPtr(branch.first, branch.second, branch.second));

            auto helper = std::make_shared<BnBHelper>(copy_model, printer, printer.GetBest());
            printer.AddCallback(helper);
            helper->BnB(solution);

            std::cout << "Task finished: " << ++finished_index << " / " << tasks_total << std::endl << std::endl;
        });
    }
    tasks_total = tasks.size();

    ParallelExecturor executor(std::move(tasks));
    executor.WaitForJobDone();

    std::cout << "Total branches: " << bnbh.GetBranchCount() << std::endl << std::endl;

    return{};

}
