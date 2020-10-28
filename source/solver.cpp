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
    double                 avg = 0;
    uint32_t               cnt = 0;

    std::vector<uint32_t> vars_stats;

    struct IndependetConstrain
    {
        std::set<uint32_t> nodes;
        std::unique_ptr<ConstrainsGuard> constrain;
    };

    std::vector<IndependetConstrain> additional_constrains;

    struct BranchingData
    {
        double way = 0.0;
        Solution res;
    };
    using Branches = std::pair<BranchingData, BranchingData>;

    void Branching(BranchingData& way, size_t vertex)
    {
        ++depth;
        auto guard = model.AddScopedConstrains(vertex, way.way, way.way);
        BnC();
        --depth;
    }

    Branches GetWays(size_t branch_index)
    {
        Branches ways = { { 0.0 }, { 1.0 } };
        //{
        //    auto guard = model.AddScopedConstrains(branch_index, ways.first.way, ways.first.way);
        //    ways.first.res = Solve();
        //}
        //{
        //    auto guard = model.AddScopedConstrains(branch_index, ways.second.way, ways.second.way);
        //    ways.second.res = Solve();
        //}
        //
        //bool go_right = ways.first.res.integer_count > ways.second.res.integer_count;
        //if (ways.first.res.integer_count == ways.second.res.integer_count)
        //    go_right = ways.first.res.upper_bound > ways.second.res.upper_bound;
        //
        //if (!go_right)
        //    std::swap(ways.first, ways.second);
        //
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

    bool Separation(const Solution& solution)
    {
        auto heur_constrs = graph.GetHMIS(solution.variables);
        //auto heur_constrs = graph.GetWeightHeuristicConstrFor(solution.branching_index, solution.variables);
        for (auto& constr : heur_constrs)
        {
            for (auto v : constr.nodes)
                ++vars_stats[v];

            additional_constrains.emplace_back();
            additional_constrains.back().constrain = model.AddScopedConstrain(constr.nodes);
            additional_constrains.back().nodes = constr.nodes;
        }

        return !heur_constrs.empty();
    }

    void CleanUp(const Solution& solution)
    {
        //std::erase_if(additional_constrains, [&solution](const IndependetConstrain& constr) {
        //    double sum = 0;
        //    for (auto x : constr.nodes)
        //        sum += solution.variables[x];
        //
        //    return EpsValue(sum) < 0.5;
        //});
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

    void BnC()
    {
        Solution solution = Solve();
        ++branches_count;

        if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
            return;

        std::deque<double> sols;
        Timer timer_total;
        while (true)
        {
            Timer timer;
            auto heur_constrs = graph.GetWeightHeuristicConstr(solution.variables);
            //auto heur_constrs = graph.GetHMIS(solution.variables);
            auto elapsed = timer.Stop();
            if (cnt++ == 0)
            {
                avg = double(elapsed);
            }
            else
            {
                avg += (1.0 / (double(cnt))) * (double(elapsed) - avg);
            }

            if (heur_constrs.empty())
                break;

            sols.emplace_back(solution.upper_bound);
            if (sols.size() > 10)
            {
                sols.pop_front();
                if (sols.front() - sols.back() < 0.1)
                    break;
            }

            std::vector<IndependetConstrain> additional_constrains_t;
            for (auto& constr : heur_constrs)
            {
                additional_constrains_t.emplace_back();
                additional_constrains_t.back().constrain = model.AddScopedConstrain(constr.nodes);
                additional_constrains_t.back().nodes = constr.nodes;
            }

            solution = Solve();

            std::erase_if(additional_constrains_t, [&solution](const IndependetConstrain& constr) {
                double sum = 0;
                for (auto x : constr.nodes)
                    sum += solution.variables[x];

                return EpsValue(sum) < 1.0;
            });

            MoveAppend(additional_constrains_t, additional_constrains);

            if (sols.back() - solution.upper_bound < 0.001)
                break;

            if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
                return;
        }
        CleanUp(solution);

        std::cout << depth << " " << EpsValue(solution.upper_bound) << " " << EpsValue(avg) << " " << timer_total.Stop() << std::endl;

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

        auto ways = GetWays(branch_index);
        Branching(ways.first, branch_index);
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
