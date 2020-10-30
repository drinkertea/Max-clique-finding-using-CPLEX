#include "utils.hpp"
#include "heuristic.hpp"
#include "model.hpp"
#include <future>

static std::atomic<uint32_t> threads_count = 0;

class BnCHelper
{
    std::atomic<uint64_t>& best_obj_value;
    ModelData              model;
    ThreadsAccumulator&    accumulator;
    uint64_t               branches_count = 0;
    const Graph&           graph;
    uint32_t               depth = 0;
    uint32_t               ones_count = 0;

public:
    BnCHelper(const BnCHelper& r, bool is_one)
        : best_obj_value(r.best_obj_value)
        , model(r.model)
        , accumulator(r.accumulator)
        , graph(r.graph)
        , depth(r.depth)
        , ones_count(r.ones_count + uint32_t(is_one))
    {
    }

private:
    AvgTimer               average_heuristic_timer;
    AvgTimer               average_loop_timer;
    AvgTimer               average_solve_timer;

    struct IndependetConstrain
    {
        std::vector<uint32_t> nodes;
        std::unique_ptr<ConstrainsGuard> constrain;
    };

    struct BranchingData
    {
        double way = 0.0;
        Solution res;
    };
    using Branches = std::pair<BranchingData, BranchingData>;

    void Branching(BranchingData& way, size_t vertex, bool is_first = false)
    {
        if (way.way == 1.0)
            ++ones_count;
        auto guard = model.AddScopedConstrains(vertex, way.way, way.way);
        BnC(is_first ? &way.res : nullptr, way.way == 1.0);
        if (way.way == 1.0)
            --ones_count;
    }

    auto BranchingAsync(const BranchingData& way, size_t vertex)
    {
        auto sub_helper = std::make_shared<BnCHelper>(*this, way.way == 1.0);
        return std::async(std::launch::async, [sub_helper, way_data = way, vertex] {
            auto& helper = *sub_helper;
            auto guard = helper.model.AddScopedConstrains(vertex, way_data.way, way_data.way);
            helper.BnC(way_data.res.variables.empty() ? nullptr : &way_data.res, way_data.way == 1.0);
            return helper.GetBranchCount();
        });
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
        TimerGuard tg(average_solve_timer);
        Solution res{};

        IloCplex cplex = model.CreateSolver();
        //if (depth > 1)
        //    cplex.setParam(IloCplex::Param::Threads, 1);
        if (!cplex.solve())
            return res;

        res.variables.resize(model.GetSize(), 0.0);
        for (size_t i = 0; i < model.GetSize(); ++i)
        {
            res.variables[i] = model.ExtractValue(cplex, i);
            if (EpsValue(res.variables[i]) == 1.0)
                ++res.integer_count;
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

    void CleanUp(const Solution& solution, size_t base, std::vector<IndependetConstrain>& additional_constrains)
    {
        additional_constrains.erase(
            std::remove_if(additional_constrains.begin() + base, additional_constrains.end(),
                [&solution](const IndependetConstrain& constr) {
                    double sum = 0;
                    for (auto x : constr.nodes)
                        sum += solution.variables[x];

                    return EpsValue(sum) < 1.0;
                }
            ),
            additional_constrains.end()
        );
    }

    void Separate(const Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        TimerGuard tg(average_heuristic_timer);
        size_t max_size = 0;
        graph.GetWeightHeuristicConstr(solution.variables, [this, &additional_constrains, &max_size](auto&& constr) {
            //max_size = std::min(max_size, constr.size());
            //if (constr.size() + 3 < max_size)
            //    return;
            additional_constrains.emplace_back();
            additional_constrains.back().constrain = model.AddScopedConstrain(constr);
            additional_constrains.back().nodes = std::move(constr);
        });
    }

    void Cutting(Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        TimerGuard tg(average_loop_timer);
        std::deque<double> sols;
        size_t base_index = additional_constrains.size();
        auto start = std::chrono::system_clock::now();
        while (true)
        {
            auto curr = std::chrono::system_clock::now();
            auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(curr - start).count();
            //if (diff > 1000)
            //    return;

            if (EpsValue(solution.upper_bound) < static_cast<double>(best_obj_value + 1))
                return;

            if (solution.branching_index == g_invalid_index)
                return;

            CleanUp(solution, base_index, additional_constrains);

            base_index = additional_constrains.size();
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

    void BnC(const Solution* initial = nullptr, bool is_one_way = false)
    {
        Solution solution = initial ? *initial : Solve();
        ++branches_count;
        auto force_cut = std::max(uint32_t(sqrt(double(graph.GetSize()))), 1u);

        std::vector<IndependetConstrain> additional_constrains;
        //if (depth < force_cut || ones_count > 0)
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

                BnC(nullptr, is_one_way);
            }
            return;
        }

        CleanUp(solution, 0, additional_constrains);

        std::stringstream ss;
        ss << depth << std::fixed << std::setprecision(2)
            << " " << threads_count
            << " " << ones_count
            << " " << additional_constrains.size()
            << " " << solution.upper_bound
            << " " << average_heuristic_timer.GetValue()
            << " " << average_loop_timer.GetValue()
            << " " << average_solve_timer.GetValue()
            << std::endl;
        accumulator.Print(ss.str());

        bool is_zero_branch = ones_count == 0;
        bool is_hard_branch = ones_count == 1 && is_one_way;
        bool is_close_to_root = depth < 4;
        depth++;
        if (is_zero_branch)
        {
            accumulator.Print("\n ZERO SPLIT \n");
            threads_count += 1;

            Branches ways = { { 1.0 }, { 0.0 } };
            auto left = BranchingAsync(ways.first, branch_index);
            Branching(ways.second, branch_index);
            left.wait();
            branches_count += left.get();

            threads_count -= 1;
        }
        else
        {
            if (is_hard_branch || is_close_to_root)
            {
                accumulator.Print("\n THREAD SPLIT \n");
                threads_count += 2;
                auto ways = GetWays(branch_index);
                auto left = BranchingAsync(ways.first, branch_index);
                auto right = BranchingAsync(ways.second, branch_index);
                left.wait();
                right.wait();
                branches_count += left.get();
                branches_count += right.get();
                threads_count -= 2;
            }
            else
            {
                auto ways = GetWays(branch_index);
                Branching(ways.first, branch_index);
                Branching(ways.second, branch_index);
            }
        }
        depth--;
    }

public:
    BnCHelper(const ModelData& m, ThreadsAccumulator& p, std::atomic<uint64_t>& best, const Graph& graph)
        : model(m)
        , accumulator(p)
        , best_obj_value(best)
        , graph(graph)
    {
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

        bool initial_better = initial_solution.branching_index == g_invalid_index &&
                              initial_solution.integer_count < heuristic_clique.size() &&
                              graph.CheckSolution(ExtractClique(initial_solution.variables)).empty();

        if (!initial_better)
        {
            for (auto& x : initial_solution.variables)
                x = 0.0;
            for (auto vert : heuristic_clique)
                initial_solution.variables[vert] = 1.0;
            initial_solution.upper_bound = static_cast<double>(heuristic_clique.size());
            initial_solution.integer_count = static_cast<uint64_t>(heuristic_clique.size());
            OnBestSolution(initial_solution);
        }

        BnC(nullptr, true);
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
