#include "graph.h"
#include "solver.h"

#include <iostream>
#include <algorithm>
#include <chrono>

std::vector<std::string> g_simple_graphs = {
    "c-fat200-1.clq",
    "c-fat200-2.clq",
    "c-fat200-5.clq",
    "c-fat500-1.clq",
    "c-fat500-10.clq",
    "c-fat500-2.clq",
    "c-fat500-5.clq",
    "MANN_a9.clq",
    "hamming6-2.clq",
    "hamming6-4.clq",
    "C125.9.clq",
    "keller4.clq",
    "brock200_1.clq",
    "brock200_2.clq",
    "brock200_3.clq",
    "brock200_4.clq",
    "gen200_p0.9_44.clq",
    "gen200_p0.9_55.clq",
    "san200_0.7_1.clq",
    "san200_0.7_2.clq",
    "san200_0.9_1.clq",
    "san200_0.9_2.clq",
    "san200_0.9_3.clq",
    "sanr200_0.7.clq",
    "sanr200_0.9.clq",
    "p_hat300-1.clq",
    "p_hat300-2.clq",
    "p_hat300-3.clq",
};

std::vector<std::string> g_medium_graphs = {
    "hamming8-2.clq",
    "hamming8-4.clq",
    "san400_0.5_1.clq",
    "san400_0.7_1.clq",
    "san400_0.7_2.clq",
    "san400_0.7_3.clq",
    "san400_0.9_1.clq",
    "sanr400_0.5.clq",
    "sanr400_0.7.clq",
    "brock400_1.clq",
    "brock400_2.clq",
    "brock400_3.clq",
    "brock400_4.clq",
    "C250.9.clq",
    "C500.9.clq",
    "p_hat500-1.clq",
    "p_hat500-2.clq",
    "p_hat500-3.clq",
    "gen400_p0.9_55.clq",
    "gen400_p0.9_65.clq",
    "gen400_p0.9_75.clq",
    "DSJC500_5.clq",
    "johnson8-2-4.clq",
    "johnson8-4-4.clq",
    "MANN_a27.clq",
};

std::vector<std::string> g_hard_graphs = {
    "p_hat1000-1.clq",
    "p_hat1000-2.clq",
    "p_hat1000-3.clq",
    "p_hat1500-1.clq",
    "p_hat1500-2.clq",
    "p_hat1500-3.clq",
    "p_hat700-1.clq",
    "p_hat700-2.clq",
    "p_hat700-3.clq",
    "san1000.clq",
    "brock800_1.clq",
    "brock800_2.clq",
    "brock800_3.clq",
    "brock800_4.clq",
    "C1000.9.clq",
    "C2000.5.clq",
    "C2000.9.clq",
    "C4000.5.clq",
    "DSJC1000_5.clq",
    "hamming10-2.clq",
    "hamming10-4.clq",
    "johnson16-2-4.clq",
    "johnson32-2-4.clq",
    "keller5.clq",
    "keller6.clq",
    "MANN_a45.clq",
    "MANN_a81.clq",
};

void Test(const std::string& file_path, const std::string& path)
{
    uint64_t int_elapsed = 0;
    uint64_t bnb_elapsed = 0;
    std::string test_result = "CRASHED";

    try
    {
        Graph g(path);

        auto int_start = std::chrono::system_clock::now();
        auto int_result = FindMaxCliqueInteger(g);
        auto int_end = std::chrono::system_clock::now();
        int_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(int_end - int_start).count();

        auto bnb_start = std::chrono::system_clock::now();
        auto bnb_result = FindMaxCliqueBnB(g);
        auto bnb_end = std::chrono::system_clock::now();
        bnb_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(bnb_end - bnb_start).count();

        test_result = int_result.size() == bnb_result.size() ? "PASSED" : "FAILED";
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
    ss << std::endl;

    std::ofstream outfile(file_path, std::ios_base::app);;
    outfile << ss.str();
}

void RunTests(const std::string& prefix, const std::vector<std::string>& tests)
{
    for (const auto& test : tests)
        Test(prefix + "_test_results.txt", "graphs/" + test);
}

int main()
{
    RunTests("simple", g_simple_graphs);
    RunTests("medium", g_medium_graphs);
    RunTests("hardcr", g_hard_graphs);
    return 0;
}
