#include "utils.hpp"
#include "heuristic.hpp"
#include "model.hpp"

class BnCHelper
{
    std::atomic<uint64_t>& best_obj_value;
    ModelData&             model;
    ThreadsAccumulator&    accumulator;
    uint64_t               branches_count = 0;
    const Graph&           graph;
    uint32_t               depth = 0;
    AvgTimer               average_heuristic_timer;
    AvgTimer               average_loop_timer;
    AvgTimer               average_solve_timer;

    std::vector<uint32_t> vars_stats;

    struct IndependetConstrain
    {
        std::vector<uint32_t> nodes;
        std::unique_ptr<ConstrainsGuard> constrain;
    };

    //std::vector<IndependetConstrain> additional_constrains;

    struct BranchingData
    {
        double way = 0.0;
        Solution res;
    };
    using Branches = std::pair<BranchingData, BranchingData>;

    void Branching(BranchingData& way, size_t vertex, bool is_first = false)
    {
        ++depth;
        auto guard = model.AddScopedConstrains(vertex, way.way, way.way);
        BnC(is_first ? &way.res : nullptr);
        --depth;
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
        double max = 0;
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

    Solution Solve()
    {
        TimerGuard tg(average_solve_timer);
        Solution res{};

        IloCplex cplex = model.CreateSolver();
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

    void OnBestSolution(const Solution& solution)
    {
        best_obj_value = static_cast<uint64_t>(solution.integer_count);
        accumulator.OnBestSolution(solution);
    }

    void MoveAppend(std::vector<IndependetConstrain>& src, std::vector<IndependetConstrain>& dst)
    {
        if (dst.empty())
        {
            dst = std::move(src);
        }
        else
        {
            dst.reserve(dst.size() + src.size());
            std::move(std::begin(src), std::end(src), std::back_inserter(dst));
            src.clear();
        }
    }

    void CleanUp(const Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        std::erase_if(additional_constrains, [&solution](const auto& constr) {
            return EpsValue(solution.GetNodesWeight(constr.nodes)) < 1.0;
        });
    }

    void Separate(const Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        TimerGuard tg(average_heuristic_timer);
        graph.GetWeightHeuristicConstrFor(solution.branching_index, solution.variables, [this, &additional_constrains](auto&& constr) {
            additional_constrains.emplace_back();
            additional_constrains.back().constrain = model.AddScopedConstrain(constr);
            additional_constrains.back().nodes = std::move(constr);
        });
    }

    void Cutting(Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        TimerGuard tg(average_loop_timer);
        std::deque<double> sols;
        while (true)
        {
            if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
                return;

            if (solution.branching_index == g_invalid_index)
                return;

            size_t base_index = additional_constrains.size();
            Separate(solution, additional_constrains);

            if (base_index == additional_constrains.size())
                return;

            sols.emplace_back(solution.upper_bound);
            if (sols.size() > 10)
            {
                sols.pop_front();
                if (sols.front() - sols.back() < 0.1)
                    return;
            }

            solution = Solve();

            if (sols.back() - solution.upper_bound < 0.001)
                return;
        }
    }

    void BnC(const Solution* initial = nullptr)
    {
        Solution solution = initial ? *initial : Solve();
        ++branches_count;

        std::vector<IndependetConstrain> additional_constrains;
        Cutting(solution, additional_constrains);

        if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
            return;

        size_t branch_index = solution.branching_index;
        if (branch_index == g_invalid_index)
        {
            auto non_edge_pairs = graph.CheckSolution(ExtractClique(solution.variables));
            if (non_edge_pairs.empty())
            {
                if (solution.integer_count >= best_obj_value)
                    OnBestSolution(solution);
            }
            else
            {
                for (auto& constr : non_edge_pairs)
                {
                    additional_constrains.emplace_back();
                    additional_constrains.back().constrain = model.AddScopedConstrain(constr);
                    additional_constrains.back().nodes = constr;
                }

                BnC();
            }
            return;
        }

        CleanUp(solution, additional_constrains);

        std::cout << depth << std::fixed << std::setprecision(2)
            << " " << additional_constrains.size()
            << " " << solution.upper_bound
            << " " << average_heuristic_timer.GetValue()
            << " " << average_loop_timer.GetValue()
            << " " << average_solve_timer.GetValue()
            << std::endl;

        auto ways = GetWays(branch_index);
        Branching(ways.first,  branch_index);
        Branching(ways.second, branch_index);
    }

public:
    BnCHelper(ModelData& m, ThreadsAccumulator& p, std::atomic<uint64_t>& best, const Graph& graph)
        : model(m)
        , accumulator(p)
        , best_obj_value(best)
        , graph(graph)
    {
        vars_stats.resize(graph.GetSize(), 0);
    };

    uint64_t GetBranchCount() const
    {
        return branches_count;
    }

    void BnC(const std::set<uint32_t>& heuristic_clique)
    {
        Timer timer;
        auto initial_solution = Solve();
        std::cout << "Initial solution:   " << initial_solution.upper_bound << std::endl;
        std::cout << "Heuristic solution: " << heuristic_clique.size() << std::endl;
        std::cout << "Time:               " << timer.Stop() << " ms" << std::endl << std::endl;

        if (initial_solution.integer_count < heuristic_clique.size())
        {
            initial_solution.branching_index = SelectBranch(initial_solution.variables);
            for (auto& x : initial_solution.variables)
                x = 0.0;
            for (auto vert : heuristic_clique)
                initial_solution.variables[vert] = 1.0;
            initial_solution.upper_bound = static_cast<double>(heuristic_clique.size());
            initial_solution.integer_count = static_cast<uint64_t>(heuristic_clique.size());
            OnBestSolution(initial_solution);
        }

        BnC();
    }
};

std::vector<uint32_t> CliqueFinder::FindMaxCliqueBnB(const Graph& graph)
{
    ModelData             model(graph, IloNumVar::Float);
    ThreadsAccumulator    accumulator(graph);
    BranchingStateTracker branching_state_tracker(graph.GetSize());
    std::atomic<uint64_t> best_obj_value = 1;
    BnCHelper             bnc_helper(model, accumulator, best_obj_value, graph);

    bnc_helper.BnC(Heuristic(graph));

    std::cout << "Total branches: " << bnc_helper.GetBranchCount() << std::endl << std::endl;

    auto max_clique = accumulator.GetBestSolution();
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
