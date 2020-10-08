#include "graph.h"
#include "solver.h"

#include <iostream>
#include <algorithm>

void HardcodedTest(const std::string& path, const std::vector<uint32_t>& vericies)
{

}

int main()
{
    try
    {
        Graph g("graphs/C125.9.clq");
        //Graph g("graphs/brock200_2.clq");
        std::vector<uint32_t> brock200_2_test = { 26, 120 ,119 ,157, 69 ,182 ,47 ,148 ,104 ,134 , 54, 144 };
        std::sort(brock200_2_test.begin(), brock200_2_test.end());
        //Graph g("graphs/keller4.clq");
        //Graph g("graphs/p_hat300-1.clq");
        //Graph g("graphs/simple.clq");
        auto result = FindMaxClique(g);

        std::cout << "Max clique verticies:" << std::endl;
        for (auto y : result)
            std::cout << y << " ";
        std::cout << std::endl;

        std::cout << "Clique equal to answer: " << int(result == brock200_2_test) << std::endl;
    }
    catch (std::runtime_error&)
    {
        std::cerr << "Something went wrong, please debug me :)" << std::endl;
    }
    return 0;
}
