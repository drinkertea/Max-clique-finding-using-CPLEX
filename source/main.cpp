#include "graph.h"
#include "solver.h"

#include <iostream>

int main()
{
    try
    {
        //Graph g("graphs/C125.9.clq");
        //Graph g("graphs/brock200_2.clq");
        //Graph g("graphs/keller4.clq");
        //Graph g("graphs/p_hat300-1.clq");
        Graph g("graphs/simple.clq");
        auto result = FindMaxClique(g);

        std::cout << "Max clique verticies:" << std::endl;
        for (auto y : result)
            std::cout << y << " ";
        std::cout << std::endl;
    }
    catch (std::runtime_error&)
    {
        std::cerr << "Something went wrong, please debug me :)" << std::endl;
    }
    return 0;
}
