#include "tridiagonal_solver.hpp"
#include <fstream>

std::vector<double> tridiag_solve(std::string name){
        std::ifstream file;
        file.open(name);
        int n;
        file >> n;
        std::vector<double> system(n * 6 + 2, 0);
        std::vector<double> x(n + 1, 0);

        for(int i = 1; i < 3 * n - 1; i++){   //a b c
                file >> system[i];
        }
        for(int i = 3 * n; i < 4 * n; i++){   //d
                file >> system[i];
        }


        for(int i = 0; i < n; i++){
                system[4 * n + i + 1] = -(system[2 * n + i] / (system[i] * system[4 * n + i] + system[n + i]));                                         //p
                system[5 * n + i + 2] = (system[3 * n + i] - system[i] * system[5 * n + i + 1]) / (system[i] * system[4 * n + i] + system[n + i]);      //q
        }

        for(int i = 0; i < n; i++){
                x[n - i - 1] = x[n - i] * system[5 * n - i] + system[6 * n + 1 - i];    //x
        }

        x.pop_back();

        file.close();

        return x;
}
