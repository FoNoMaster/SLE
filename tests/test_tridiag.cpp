#include <gtest/gtest.h>
#include <code/tridiag.hpp>
#include <cmath>

TEST(TriDiagSolve, Test){
	std::vector<double> a{1, 2, 3};
	std::vector<double> b{4, 5, 60, 7};
	std::vector<double> c{8, 9, 10};
	std::vector<double> d{20, 38, 224, 37};
	std::vector<double> ans{1, 2, 3, 4};

	TriDiagMatrix matr(a, b, c);
	std::vector<double> res = TriDiagSolve(matr, d);

	for(std::size_t i = 0; i < ans.size(); i++){
                ASSERT_NEAR(res[i], ans[i], pow(10, -15));
        }
}
