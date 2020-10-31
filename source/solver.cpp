#include "utils.hpp"
#include "heuristic.hpp"
#include "model.hpp"
#include <future>

struct Statistics
{
    uint64_t branches_count = 0;
    AvgTimer average_heuristic_timer;
    AvgTimer average_loop_timer;
    AvgTimer average_solve_timer;

    Statistics& operator+=(const Statistics& rhs)
    {
        branches_count += rhs.branches_count;
        average_heuristic_timer += rhs.average_heuristic_timer;
        average_loop_timer += rhs.average_loop_timer;
        average_solve_timer += rhs.average_solve_timer;
        return *this;
    }
};

class BnCHelper
{
    ModelData               model;
    const Graph&            graph;
    ThreadingData&          threading;
    const std::atomic_bool& stop;
    Statistics              statistics{};
    uint32_t                depth = 0;
    uint32_t                ones_count = 0;

public:
    BnCHelper(const BnCHelper& r)
        : model(r.model)
        , graph(r.graph)
        , depth(r.depth)
        , threading(r.threading)
        , stop(r.stop)
    {
    }

private:
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

    void Branching(const BranchingData& way, size_t vertex)
    {
        ++depth;
        ones_count += static_cast<uint32_t>(way.way);
        {
            auto guard = model.AddScopedConstrains(vertex, way.way, way.way);
            BnC(way.res.variables.empty() ? nullptr : &way.res);
        }
        ones_count -= static_cast<uint32_t>(way.way);
        --depth;
    }

    void BranchingZero(const BranchingData& way, size_t vertex)
    {
        if (depth > 0)
            return Branching(way, vertex);

        auto guard = threading.GetAsyncGuard(depth);
        Branching(way, vertex);
    }

    auto BranchingAsync(const BranchingData& way, size_t vertex)
    {
        auto sub_helper = std::make_shared<BnCHelper>(*this);
        return std::async(std::launch::async, [sub_helper, way_data = way, vertex] {
            auto guard = sub_helper->threading.GetAsyncGuard(sub_helper->depth);
            sub_helper->Branching(way_data, vertex);
            return sub_helper->GetStats();
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
        TimerGuard tg(statistics.average_solve_timer);
        Solution res{};

        IloCplex cplex = model.CreateSolver();
        if (threading.max_workers_count <= threading.active_workers_count)
            cplex.setParam(IloCplex::Param::Threads, 1);

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
        OnBestSolution(ExtractClique(solution.variables));
    }

    void OnBestSolution(const std::set<uint32_t>& clique)
    {
        if (clique.size() <= threading.best_obj_value)
            return;

        threading.best_obj_value = clique.size();
        threading.accumulator.OnBestSolution(clique);
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
        TimerGuard tg(statistics.average_heuristic_timer);
        graph.GetWeightHeuristicConstr(solution.variables, [this, &additional_constrains](auto&& constr) {
            additional_constrains.emplace_back();
            additional_constrains.back().constrain = model.AddScopedConstrain(constr);
            additional_constrains.back().nodes = std::move(constr);
        });
    }

    void Cutting(Solution& solution, std::vector<IndependetConstrain>& additional_constrains)
    {
        TimerGuard tg(statistics.average_loop_timer);
        std::deque<double> sols;
        size_t base_index = additional_constrains.size();
        while (true)
        {
            if (EpsValue(solution.upper_bound) < static_cast<double>(threading.best_obj_value + 1))
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

    void PrintStats(const Solution& solution, size_t constr_size)
    {
        if (depth > 0 && ones_count != 0)
            return;

        std::stringstream ss;
        ss << depth << std::fixed << std::setprecision(2)
            << " " << ones_count
            << " " << constr_size
            << " " << solution.upper_bound
            << " " << statistics.average_heuristic_timer.GetValue()
            << " " << statistics.average_loop_timer.GetValue()
            << " " << statistics.average_solve_timer.GetValue()
            << std::endl;
        threading.accumulator.Print(ss.str());
    }

    void BnC(const Solution* initial = nullptr)
    {
        if (stop)
            return;

        Solution solution = initial ? *initial : Solve();
        ++statistics.branches_count;

        if (EpsValue(solution.upper_bound) < static_cast<double>(threading.best_obj_value + 1))
            return;

        std::vector<IndependetConstrain> additional_constrains;
        additional_constrains.reserve(graph.GetSize() * 10);
        Cutting(solution, additional_constrains);

        if (EpsValue(solution.upper_bound) < static_cast<double>(threading.best_obj_value + 1))
            return;

        size_t branch_index = solution.branching_index;
        if (branch_index == g_invalid_index)
        {
            auto possible_clique = ExtractClique(solution.variables);
            auto non_edge_pairs  = graph.CheckSolution(possible_clique);
            if (non_edge_pairs.empty())
            {
                OnBestSolution(possible_clique);
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

        CleanUp(solution, 0, additional_constrains);
        PrintStats(solution, additional_constrains.size());

        if (ones_count == 0)
        {
            Branches ways = { { 1.0 }, { 0.0 } };
            auto left = BranchingAsync(ways.first, branch_index);
            BranchingZero(ways.second, branch_index);

            statistics += left.get();
        }
        else
        {
            auto ways = GetWays(branch_index);
            Branching(ways.first, branch_index);
            Branching(ways.second, branch_index);
        }
    }

public:
    BnCHelper(const ModelData& m, const Graph& graph, const std::atomic_bool& stop, ThreadingData& data)
        : model(m)
        , graph(graph)
        , threading(data)
        , stop(stop)
    {
    };

    const Statistics& GetStats() const
    {
        return statistics;
    }

    void BnC(const std::set<uint32_t>& heuristic_clique)
    {
        auto initial_solution = Solve();
        std::cout << "Initial solution:   " << initial_solution.upper_bound << std::endl;
        std::cout << "Heuristic solution: " << heuristic_clique.size() << std::endl;

        bool initial_better = initial_solution.branching_index == g_invalid_index &&
                              initial_solution.integer_count < heuristic_clique.size() &&
                              graph.CheckSolution(ExtractClique(initial_solution.variables)).empty();

        if (!initial_better)
            OnBestSolution(heuristic_clique);

        BnC();
    }
};

std::vector<uint32_t> CliqueFinder::FindMaxCliqueBnB(const Graph& graph)
{
    ModelData             model(graph, IloNumVar::Float);
    ThreadingData         threading_data{ graph };
    BnCHelper             bnc_helper(model, graph, stop, threading_data);

    bnc_helper.BnC(Heuristic(graph));

    const auto& stats      = bnc_helper.GetStats();
    branch_count           = stats.branches_count;
    average_heuristic_time = stats.average_heuristic_timer.GetValue();
    average_loop_time      = stats.average_loop_timer.GetValue();
    average_solve_time     = stats.average_solve_timer.GetValue();

    auto max_clique = threading_data.accumulator.GetBestSolution();
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
