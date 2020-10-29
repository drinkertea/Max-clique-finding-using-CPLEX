#pragma once
#include "common.hpp"

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

struct BranchingStateTracker
{
    BranchingStateTracker(size_t md)
        : max_depth(md)
    {
    }

    DepthGuard OnBranch(const Solution& solution, size_t vertex, double constr)
    {
        curr_path.emplace_back(vertex, constr);
        if (solution.integer_count >= 1 || curr_path.size() >= max_depth)
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
    Path                     curr_path;
};

struct ThreadsAccumulator
{
    ThreadsAccumulator(const Graph& g)
        : graph(g)
    {
    }

    void OnBestSolution(const Solution& solution)
    {
        OnBestSolution(ExtractClique(solution.variables));
    }

    void OnBestSolution(const std::set<uint32_t>& clique)
    {
        {
            std::lock_guard<std::mutex> lg(solutions_mutex);
            solutions.emplace_back(clique);
        }

        std::stringstream ss;
        ss << "Found " << clique.size() << " " << std::endl;

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

    std::set<uint32_t> GetBestSolution() const
    {
        std::set<uint32_t> best;
        for (const auto& solution : solutions)
        {
            if (solution.size() <= best.size())
                continue;

            best = solution;
        }
        return best;
    }

private:
    std::atomic<size_t> finished_index = 0;
    std::atomic<size_t> tasks_total = 0;

    std::vector<std::set<uint32_t>> solutions;
    std::mutex            solutions_mutex;

    std::mutex print_mutex;
    const Graph& graph;
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

struct TimeoutThread
{
    std::atomic_bool stopped = false;

    TimeoutThread(std::chrono::milliseconds wait_for)
        : m_thread([this, wait_for]() {
            std::this_thread::sleep_for(wait_for);
            stopped = true;
        })
    {
    }

    ~TimeoutThread()
    {
        if (!stopped)
            TerminateThread(m_thread.native_handle(), 0);
        m_thread.join();
    }

private:
    std::thread m_thread;
};

struct Timer
{
    Timer()
    {
        Reset();
    }

    void Reset()
    {
        start = std::chrono::system_clock::now();
    }

    uint64_t Stop()
    {
        auto int_end = std::chrono::system_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(int_end - start).count();
    }

private:
    std::chrono::system_clock::time_point start;
};

struct AvgTimer
{
    AvgTimer()
    {
    }

    void Reset()
    {
        timer.Reset();
    }

    void Stop()
    {
        auto elapsed = timer.Stop();
        if (cnt++ == 0)
        {
            avg = double(elapsed);
        }
        else
        {
            avg += (1.0 / (double(cnt))) * (double(elapsed) - avg);
        }
    }

    double GetValue() const
    {
        return avg;
    }

private:
    Timer    timer{};
    double   avg = 0;
    uint32_t cnt = 0;
};

struct TimerGuard
{
    TimerGuard(AvgTimer& t)
        : timer(t)
    {
        timer.Reset();
    }

    ~TimerGuard()
    {
        timer.Stop();
    }

private:
    AvgTimer& timer;
};

struct NonEdgeKHelper
{
    uint32_t                max_depth = 0;
    uint32_t                m_break_count = 1000;
    size_t                  m_size = 0;
    const Graph& m_graph;
    const std::atomic_bool& stop;

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
