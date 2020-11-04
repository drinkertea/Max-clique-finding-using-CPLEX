#pragma once
#include "common.hpp"

struct ThreadsAccumulator
{
    ThreadsAccumulator(const Graph& g)
        : graph(g)
    {
    }

    void OnBestSolution(const std::set<uint32_t>& clique)
    {
        {
            std::lock_guard<std::mutex> lg(solution_mutex);
            if (clique.size() > best_solution.size())
                best_solution = clique;
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

    void Print(const std::string& str)
    {
        std::lock_guard<std::mutex> lg(print_mutex);
        std::cout << str << std::endl;
    }

    const std::set<uint32_t>& GetBestSolution() const
    {
        return best_solution;
    }

private:
    std::atomic<size_t> finished_index = 0;
    std::atomic<size_t> tasks_total = 0;
    std::set<uint32_t>  best_solution;
    std::mutex          solution_mutex;

    std::mutex print_mutex;
    const Graph& graph;
};

struct ThreadingData
{
    ThreadsAccumulator      accumulator;
    std::atomic<uint64_t>   best_obj_value = 1;

    std::condition_variable cv;
    std::mutex              cv_m;
    std::atomic<uint32_t>   active_workers_count = 0;
    const uint32_t          max_workers_count = std::thread::hardware_concurrency();

    struct ScopeGuard
    {
        ThreadingData& threading;
        uint32_t depth = 0;
        ScopeGuard(ThreadingData& data, uint32_t depth)
            : threading(data)
            , depth(depth)
        {
            std::unique_lock<std::mutex> lk(threading.cv_m);
            threading.cv.wait(lk, [this] {
                return threading.active_workers_count < threading.max_workers_count;
            });
            ++threading.active_workers_count;

            std::stringstream ss;
            ss << std::endl << depth << ": started " << threading.active_workers_count << " / " << threading.max_workers_count;
            threading.accumulator.Print(ss.str());
        }

        ~ScopeGuard()
        {
            std::stringstream ss;
            ss << std::endl << depth << ": finished " << threading.active_workers_count << " / " << threading.max_workers_count;
            threading.accumulator.Print(ss.str());

            --threading.active_workers_count;
            threading.cv.notify_one();
        }
    };

    ScopeGuard GetAsyncGuard(uint32_t depth)
    {
        return { *this, depth };
    }
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

    AvgTimer& operator+=(const AvgTimer& rhs)
    {
        uint32_t sum_cnt = cnt + rhs.cnt;
        avg = avg * (double(cnt) / double(sum_cnt)) + rhs.avg * (double(rhs.cnt) / double(sum_cnt));
        return *this;
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
