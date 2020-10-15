#include "utils.hpp"
#include "heuristic.hpp"
#include "model.hpp"

class BnBHelper
{
    std::atomic_bool&      stop;
    std::atomic<uint64_t>& best_obj_value;
    ModelData&             model;
    ThreadsAccumulator&    accumulator;
    BranchingStateTracker* branching_state_tracker = nullptr;
    uint64_t               branches_count = 0;

    struct BranchingData
    {
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
        Timer timer;
        auto initial_solution = Solve();
        std::cout << "Initial solution:   " << initial_solution.upper_bound << std::endl;
        std::cout << "Heuristic solution: " << heuristic_clique.size() << std::endl;
        std::cout << "Time:               " << timer.Stop() << " ms" << std::endl << std::endl;

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

std::vector<uint32_t> CliqueFinder::FindMaxCliqueBnB(const Graph& graph)
{
    ModelData             model(graph, IloNumVar::Float);
    ThreadsAccumulator    accumulator(graph);
    BranchingStateTracker branching_state_tracker(graph.GetSize());
    std::atomic<uint64_t> best_obj_value = 1;
    BnBHelper             bnb_helper(model, accumulator, best_obj_value, stop, &branching_state_tracker);

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
