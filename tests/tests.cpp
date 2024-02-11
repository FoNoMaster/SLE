#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include "../src/tridiagonal_solver.hpp"

TEST(Test, Test1){
	std::vector<double> x = tridiag_solve("../data/tridiag_1.txt");		//solution
	x.push_back(0);	
	x.emplace(x.begin(), 0);
	
	std::ifstream file;
        file.open("../data/tridiag_1.txt");
        int n;
        file >> n;
        std::vector<double> system(n * 4, 0);

        for(int i = 1; i < 3 * n - 1; i++){   //a b c
                file >> system[i];
        }
        for(int i = 3 * n; i < 4 * n; i++){   //d
                file >> system[i];
        }
	
	std::vector<double> d(n, 0); //new d = Ax
	for(int i = 0; i < n; i++){
		d[i] = system[i] * x[i] + system[n + i] * x[i + 1] + system[2 * n + i] * x[i + 2];
	}

	for(int i = 0; i < n; i++){
		ASSERT_NEAR(d[i], system[3 * n + i], pow(10, -14));	//new d == d
	}
}
