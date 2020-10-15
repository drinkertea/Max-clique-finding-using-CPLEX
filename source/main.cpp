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

    std::stringstream ss;
    ss << "Test graph:                      " << path << std::endl;
    ss << "Result:                          " << test_result << std::endl;
    ss << "CPLEX integer algorithm time:    " << int_elapsed << " ms" << std::endl;
    ss << "Branch and Bound algorithm time: " << bnb_elapsed << " ms" << std::endl;
    ss << "Branch and Bound branch count:   " << cf.branch_count << std::endl;
    ss << "Max clique size:                 " << max_clique.size() << std::endl;
    ss << "Max clique:                      ";
    for (auto v : max_clique)
        ss << v << " ";
    ss << std::endl;
    ss << std::endl;

    std::ofstream outfile(file_path, std::ios_base::app);;
    outfile << ss.str();
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
