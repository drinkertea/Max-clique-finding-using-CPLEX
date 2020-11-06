#include "graph.h"
#include "solver.h"
#include "utils.hpp"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <Windows.h>

void Test(const std::string& file_path, const std::string& path)
{
    std::cout << "Test graph: " << path << std::endl;

    uint64_t int_elapsed = 0;
    uint64_t bnb_elapsed = 0;
    std::string test_result = "CRASHED";
    std::vector<uint32_t> max_clique;

    CliqueFinder cf{ false };
    try
    {
        Graph g(path);

        Timer timer;
        TimeoutThread break_timer(std::chrono::seconds(7200));
        std::vector<uint32_t> int_result = cf.FindMaxCliqueInteger(g);
        int_elapsed = timer.Stop();

        timer.Reset();
        auto bnb_result = cf.FindMaxCliqueBnB(g);
        bnb_elapsed = timer.Stop();

        test_result = int_result.size() == bnb_result.size() ? "PASSED" : "FAILED";
        if (cf.stop)
            test_result += " INTERRUPTED";

        max_clique = bnb_result;
    }
    catch (std::runtime_error&)
    {
        std::cerr << "Something went wrong, please debug me :)" << std::endl;
    }

    {
        std::stringstream ss_far;
        ss_far << "Test graph:                      " << path << std::endl;
        ss_far << "Result:                          " << test_result << std::endl;
        ss_far << "CPLEX integer algorithm time:    " << int_elapsed << " ms" << std::endl;
        ss_far << "Branch and Bound algorithm time: " << bnb_elapsed << " ms" << std::endl;
        ss_far << "Branch and Bound branch count:   " << cf.branch_count << std::endl;
        ss_far << "Average heurisrtic time:         " << cf.average_heuristic_time << std::endl;
        ss_far << "Average solve time:              " << cf.average_solve_time << std::endl;
        ss_far << "Average loop time:               " << cf.average_loop_time << std::endl;
        ss_far << "Max clique size:                 " << max_clique.size() << std::endl;
        ss_far << "Max clique:                      ";
        for (auto v : max_clique)
            ss_far << v << " ";
        ss_far << std::endl;
        ss_far << std::endl;

        std::ofstream outfile(file_path, std::ios_base::app);
        outfile << ss_far.str();
    }
    {
        std::stringstream ss_short;
        ss_short << path << ";" << bnb_elapsed;
        ss_short << std::endl;

        std::ofstream outfile_short("short_" + file_path, std::ios_base::app);
        outfile_short << ss_short.str();
    }
}

std::vector<std::string> ReadTasks(const std::string& path)
{
    std::ifstream myfile(path);
    std::string line;
    std::vector<std::string> tasks;
    while (std::getline(myfile, line))
    {
        tasks.push_back(line);
    }
    return tasks;
}

void RunTests(const std::string& prefix)
{
    for (const auto& test : ReadTasks("tasks/" + prefix + ".txt"))
        Test(prefix + "_test_results.txt", "graphs/" + test);
}

int main()
{
    RunTests("custom");
    RunTests("simple");
    RunTests("medium");
    RunTests("hardcr");
    return 0;
}
