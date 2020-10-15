#include "solver.h"
#include "graph.h"

#define NOMINMAX
#include <Windows.h>
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

constexpr size_t g_invalid_index = std::numeric_limits<size_t>::max();

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

std::set<uint32_t> ExtractClique(const std::vector<double>& variables)
{
    std::set<uint32_t> clique;
    for (uint64_t i = 0; i < variables.size(); ++i)
        if (EpsValue(variables[i]) == 1.0)
            clique.emplace(i);
    return clique;
}

struct Solution
{
    double              upper_bound = 0.0;
    std::vector<double> variables;
    size_t              branching_index = g_invalid_index;
    uint64_t            integer_count = 0;
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
        , m_size(graph.GetSize())
        , m_type(type)
        , m_model(m_env)
        , m_variables(m_env, m_size)
    {
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::default);
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::maxdegree);
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::mindegree);
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::random);
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::maxdegree_random);
        m_sum_vert_less_one += m_graph.GetHeuristicConstr(Graph::ColorizationType::mindegree_random);
        AddNonEdgePairs();
        AddIndependetConst(Graph::ColorizationType::default);
        AddIndependetConst(Graph::ColorizationType::maxdegree);
        AddIndependetConst(Graph::ColorizationType::mindegree);
        AddIndependetConst(Graph::ColorizationType::random);
        AddIndependetConst(Graph::ColorizationType::maxdegree_random);
        AddIndependetConst(Graph::ColorizationType::mindegree_random);
        InitModel(m_sum_vert_less_one);
    }

    ModelData(const ModelData& r)
        : m_graph(r.m_graph)
        , m_size(r.m_size)
        , m_type(r.m_type)
        , m_model(m_env)
        , m_variables(m_env, m_size)
    {
        InitModel(r.m_sum_vert_less_one);
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

    double ExtractValue(const IloCplex& cplex, size_t variable_index) const
    {
        return cplex.getValue(m_variables[variable_index]);
    }

    size_t GetSize() const
    { 
        return m_size; 
    }

    bool AddHeuristicConstrains(uint32_t max_depth)
    {
        // If v1 .. vk all has not edge -> v1 + v2 + ... + vk <= 1
        uint32_t break_count = uint32_t(double(m_size) * std::exp(3.14159));
        // searching for optimal "k" < "m_break_count" to no add too much constrains
        std::atomic_bool stop = false;
        NonEdgeKHelper helper{ max_depth, break_count, m_size, m_graph, stop };

        auto int_start = std::chrono::system_clock::now();
        std::thread timer_thread([&stop]() {
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(500ms);
            stop = true;
        });

        helper.AddNonEdgeRec();

        if (stop)
        {
            std::cout << "AddHeuristicConstrains interrupted" << std::endl;
            timer_thread.join();
            return false;
        }

        std::set<std::set<uint32_t>> add_constr;
        for (int i = 3; i < m_size; ++i)
        {
            if (helper.res[i].size() >= break_count)
                continue;

            for (const auto& constr : helper.res[i])
                add_constr.emplace(std::set<uint32_t>{constr.begin(), constr.end()});
        }
        if (!add_constr.empty())
            AddConstrains(add_constr);

        auto int_end = std::chrono::system_clock::now();
        auto int_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(int_end - int_start).count();
        std::cout << "AddHeuristicConstrains " << int_elapsed << " ms" << std::endl;

        TerminateThread(timer_thread.native_handle(), 0);
        timer_thread.join();
        return true;
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
            std::string name = "y[" + std::to_string(i + 1) + "]";
            m_variables[i] = IloNumVar(m_env, 0, 1, m_type, name.c_str());
        }

        IloExpr obj_expr(m_env);
        for (uint32_t i = 0; i < m_size; ++i)
            obj_expr += m_variables[i];
        IloObjective obj(m_env, obj_expr, IloObjective::Maximize);

        AddConstrains(m_sum_vert_less_one);
        m_model.add(obj);
    }

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
                if (m_graph.HasEdge(i, j))
                    continue;

                m_sum_vert_less_one.emplace(std::set<uint32_t>{ i, j });
            }
        }
    }

    struct NonEdgeKHelper
    {
        uint32_t          max_depth = 0;
        uint32_t          m_break_count = 1000;
        size_t            m_size = 0;
        const Graph&      m_graph;
        std::atomic_bool& stop;

        std::vector<std::vector<std::vector<uint32_t>>> res;

        void AddNonEdgeRec()
        {
            check.clear();
            check.reserve(m_size);
            res.resize(m_size);
            for (auto& x : res)
                x.reserve(m_break_count);
            AddNonEdgeRec(0);
        }

        std::vector<uint32_t> check;

        bool Skip(uint32_t start)
        {
            for (auto x : check)
            {
                if (m_graph.HasEdge(x, start))
                    return true;
            }
            return false;
        }

        void AddNonEdgeRec(uint32_t start)
        {
            if (stop || check.size() > max_depth)
                return;

            if (res[check.size()].size() < m_break_count)
                res[check.size()].emplace_back(check);

            for (uint32_t i = start; i < m_size; ++i)
            {
                if (Skip(i))
                    continue;

                check.push_back(i);
                AddNonEdgeRec(i + 1);
                check.pop_back();
            }
        }
    };


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

struct BranchingStateTracker
{
    BranchingStateTracker(size_t ones_cnt, size_t md)
        : max_depth(md)
        , ones_cnt(ones_cnt)
    {
    }

    DepthGuard OnBranch(const Solution& solution, size_t vertex, double constr)
    {
        curr_path.emplace_back(vertex, constr);
        if (solution.integer_count >= ones_cnt || curr_path.size() >= max_depth)
        {
            branches.emplace_back(curr_path, solution);
            return DepthGuard(curr_path, true);
        }
        return DepthGuard(curr_path, false);
    }

    const std::deque<std::pair<Path, Solution>>& GetBranches()
    {
        return branches;
    }

private:
    std::deque<std::pair<Path, Solution>> branches;
    size_t                   max_depth = 0;
    size_t                   ones_cnt = 0;
    Path                     curr_path;
};

class BnBHelper;

struct ThreadsAccumulator
{
    ThreadsAccumulator(const Graph& g)
        : graph(g)
    {
    }

    void OnBestSolution(const Solution& solution)
    {
        auto clique = ExtractClique(solution.variables);
        {
            std::lock_guard<std::mutex> lg(solutions_mutex);
            solutions.emplace_back(solution);
        }

        std::stringstream ss;
        ss << "Found " << solution.integer_count << " " << std::endl;

        for (auto v : clique)
            ss << v << " ";

        ss << std::endl;
        ss << "Valid " << (graph.IsClique(clique) ? "True" : "False") << std::endl;

        ss << std::endl;
        ss << std::endl;

        Print(ss.str());
    }

    void ProgressBegin(size_t total)
    {
        tasks_total = total;
    }

    void OnTaskComplete(const Path& path)
    {
        std::stringstream ss;
        ss << "Task finished: " << ++finished_index << " / " << tasks_total << "Path";
        for (const auto& branch : path)
            ss << " " << branch.second;
        ss << std::endl;
        Print(ss.str());
    }

    void Print(const std::string& str)
    {
        std::lock_guard<std::mutex> lg(print_mutex);
        std::cout << str << std::endl;
    }

    Solution GetBestSolution() const
    {
        Solution best{};
        for (const auto& solution : solutions)
        {
            if (solution.integer_count <= best.integer_count)
                continue;

            best = solution;
        }
        return best;
    }

private:
    std::atomic<size_t> finished_index = 0;
    std::atomic<size_t> tasks_total = 0;

    std::vector<Solution> solutions;
    std::mutex            solutions_mutex;

    std::mutex print_mutex;
    const Graph& graph;
};

class BnBHelper
{
    std::atomic_bool&      stop;
    std::atomic<uint64_t>& best_obj_value;
    uint64_t               branches_count = 0;
    ModelData&             model;
    ThreadsAccumulator&    accumulator;
    BranchingStateTracker* branching_state_tracker = nullptr;

    struct BranchingData
    {
        using CRef = std::reference_wrapper<const BranchingData>;
        double way = 0.0;
        Solution res;
    };
    using Branches = std::pair<BranchingData, BranchingData>;

    void Branching(const BranchingData& way, size_t vertex)
    {
        auto guard = model.AddScopedConstrains(vertex, way.way, way.way);
        BnB(way.res);
    }

    void Branching(BranchingStateTracker& callback, const BranchingData& way, size_t vertex)
    {
        auto guard = callback.OnBranch(way.res, vertex, way.way);
        if (guard.stop)
            return;
        Branching(way, vertex);
    }

    void Branching(BranchingStateTracker* callback, const BranchingData& way, size_t vertex)
    {
        if (callback)
            return Branching(*callback, way, vertex);
        Branching(way, vertex);
    }

    Branches GetWays(size_t branch_index)
    {
        Branches ways = { { 1.0 }, { 0.0 } };
        {
            auto guard = model.AddScopedConstrains(branch_index, ways.first.way, ways.first.way);
            ways.first.res = Solve();
        }
        {
            auto guard = model.AddScopedConstrains(branch_index, ways.second.way, ways.second.way);
            ways.second.res = Solve();
        }

        bool go_right = ways.first.res.integer_count > ways.second.res.integer_count;
        if (ways.first.res.integer_count == ways.second.res.integer_count)
            go_right = ways.first.res.upper_bound > ways.second.res.upper_bound;

        if (!go_right)
            std::swap(ways.first, ways.second);

        return ways;
    }

    size_t SelectBranch(const std::vector<double>& vars) const
    {
        double min = std::numeric_limits<double>::max();
        size_t res = g_invalid_index;

        for (size_t i = 0; i < vars.size(); ++i)
        {
            auto val = vars[i];
            if (IsInteger(val))
                continue;

            if (val >= min)
                continue;

            min = val;
            res = i;
        }

        return res;
    }

    Solution Solve()
    {
        Solution res{};

        IloCplex cplex = model.CreateSolver();
        if (!branching_state_tracker)
            cplex.setParam(IloCplex::Param::Threads, 1);
        if (!cplex.solve())
            return res;

        res.variables.resize(model.GetSize(), 0.0);
        for (size_t i = 0; i < model.GetSize(); ++i)
        {
            res.variables[i] = model.ExtractValue(cplex, i);
            res.integer_count += uint64_t(EpsValue(res.variables[i]) == 1.0);
        }

        res.branching_index = SelectBranch(res.variables);
        res.upper_bound = cplex.getObjValue();
        cplex.end();
        return res;
    }

public:
    BnBHelper(ModelData& m, ThreadsAccumulator& p, std::atomic<uint64_t>& best, std::atomic_bool& stop, BranchingStateTracker* callback = nullptr)
        : model(m)
        , accumulator(p)
        , branching_state_tracker(callback)
        , best_obj_value(best)
        , stop(stop)
    {
    };

    uint64_t GetBranchCount() const
    {
        return branches_count;
    }

    void BnB(const Solution& solution)
    {
        if (stop)
            return;

        ++branches_count;

        if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
            return;

        if (IsInteger(solution.upper_bound) && solution.integer_count >= best_obj_value)
        {
            best_obj_value = static_cast<uint64_t>(solution.integer_count);
            accumulator.OnBestSolution(solution);
        }

        size_t branch_index = solution.branching_index;
        if (branch_index == g_invalid_index)
            return;

        auto ways = GetWays(branch_index);
        Branching(branching_state_tracker, ways.first, branch_index);
        Branching(branching_state_tracker, ways.second, branch_index);
    }

    void BnB(const std::set<uint32_t>& heuristic_clique)
    {
        auto int_start = std::chrono::system_clock::now();
        auto initial_solution = Solve();
        auto int_end = std::chrono::system_clock::now();
        auto int_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(int_end - int_start).count();
        std::cout << "Initial solution: " << initial_solution.upper_bound << std::endl;
        std::cout << "Heuristic solution: " << heuristic_clique.size() << std::endl;
        std::cout << "Time: " << int_elapsed << " ms" << std::endl << std::endl;

        if (static_cast<uint32_t>(initial_solution.upper_bound) + 1 >
            static_cast<uint32_t>(1.2 * double(heuristic_clique.size())))
        {
            std::cout << "Too bad, try to add more constrains..." << std::endl;
            if (model.AddHeuristicConstrains(static_cast<uint32_t>(initial_solution.upper_bound) + 1))
            {
                initial_solution = Solve();
                std::cout << "Seconditional solution: " << initial_solution.upper_bound << std::endl << std::endl;
            }
        }

        if (initial_solution.integer_count < heuristic_clique.size())
        {
            initial_solution.branching_index = SelectBranch(initial_solution.variables);

            for (auto vert : heuristic_clique)
                initial_solution.variables[vert] = 1.0;
            initial_solution.upper_bound = static_cast<double>(heuristic_clique.size());
            initial_solution.integer_count = static_cast<uint64_t>(heuristic_clique.size());
        }

        BnB(initial_solution);
    }
};

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

bool Heuristic(const Graph& graph, std::set<uint32_t>& curr_clique, Graph::ColorizationType strategy)
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
    if (!Heuristic(graph.GetSubGraph(min_degree_vert), curr_clique, strategy))
    {
        for (uint32_t i = 0; i < graph.GetSize(); ++i)
        {
            bool is_clique_i = true;
            for (auto v : curr_clique)
            {
                if (!graph.HasEdge(i, min_degree_vert))
                {
                    is_clique_i = false;
                    break;
                }
            }

            if (!is_clique_i)
                continue;

            curr_clique.emplace(i);
            return true;
        }
    }
    return true;
}

std::set<uint32_t> Heuristic(const Graph& graph)
{
    std::set<uint32_t> best_clique;
    auto strategies = {
        Graph::ColorizationType::default,
        Graph::ColorizationType::random,
        Graph::ColorizationType::maxdegree_random,
        Graph::ColorizationType::mindegree_random,
    };
    for (auto str : strategies)
    {
        std::set<uint32_t> curr_clique;
        Heuristic(graph, curr_clique, str);
        if (curr_clique.size() > best_clique.size())
            best_clique = std::move(curr_clique);
    }
    return best_clique;
}

std::vector<uint32_t> CliqueFinder::FindMaxCliqueInteger(const Graph& graph)
{
    ModelData model(graph, IloNumVar::Int);

    IloCplex cplex = model.CreateSolver();
    cplex.exportModel("model.lp");

    bool solved = cplex.solve();
    std::vector<uint32_t> res;
    std::vector<double> vars;

    auto n = graph.GetSize();
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

std::vector<uint32_t> CliqueFinder::FindMaxCliqueBnB(const Graph& graph)
{
    ModelData             model(graph, IloNumVar::Float);
    ThreadsAccumulator    accumulator(graph);
    BranchingStateTracker branching_state_tracker(1, graph.GetSize());

    std::atomic<uint64_t> best_obj_value = 1;
    BnBHelper bnb_helper(model, accumulator, best_obj_value, stop, &branching_state_tracker);
    bnb_helper.BnB(Heuristic(graph));

    ParallelExecturor::Tasks tasks;
    std::atomic<uint64_t>    total_branches = bnb_helper.GetBranchCount();
    for (const auto& first_branches : branching_state_tracker.GetBranches())
    {
        const auto& path = first_branches.first;
        const auto& solution = first_branches.second;
        tasks.push_back([this, path, solution, &model, &accumulator, &total_branches, &best_obj_value]() {
            ModelData copy_model(model);

            std::vector<std::unique_ptr<ConstrainsGuard>> path_constr;
            for (const auto& branch : path)
                path_constr.push_back(copy_model.AddScopedConstrainsPtr(branch.first, branch.second, branch.second));

            auto helper = std::make_shared<BnBHelper>(copy_model, accumulator, best_obj_value, stop);
            helper->BnB(solution);

            accumulator.OnTaskComplete(path);
            total_branches += helper->GetBranchCount();
        });
    }
    accumulator.ProgressBegin(tasks.size());

    ParallelExecturor executor(std::move(tasks));
    executor.WaitForJobDone();

    std::cout << "Total branches: " << total_branches << std::endl << std::endl;
    branch_count = total_branches;

    auto max_clique = ExtractClique(accumulator.GetBestSolution().variables);
    return { max_clique.begin(), max_clique.end() };
}
