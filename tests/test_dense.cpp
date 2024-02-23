#include <gtest/gtest.h>
#include <code/dense.hpp>
#include <cmath>


TEST(Dense, right_product){
        std::vector<std::vector<double>> a{{1, 2}, {3, 4}};
	DenseMatrix A(a);
        A = A * 2;
        std::vector<double> b{2, 4, 6, 8};
        for(std::size_t i = 0; i < b.size(); i++){
                ASSERT_NEAR(A.get_values()[i], b[i], pow(10, -15));
        }
}


TEST(Dense, left_product){
        std::vector<std::vector<double>> a{{1, 2}, {3, 4}};
        DenseMatrix A(a);
        A = 2 * A;
        std::vector<double> b{2, 4, 6, 8};
        for(std::size_t i = 0; i < b.size(); i++){
                ASSERT_NEAR(A.get_values()[i], b[i], pow(10, -15));
        }
}


TEST(Dense, right_multiply_by_vector){
	std::vector<std::vector<double>> a{{1, 2, 0, 3, 5}, {4, 0, 6, 0, 13}, {7, 0, 0, 8, 6}, {0, 3, 0, 0, 9}};
	std::vector<double> b{5, 7, 0, 3, 0};
	DenseMatrix A(a);
	std::vector<double> res = A * b;
	std::vector<double> real_res{28, 20, 59, 21};
	for(std::size_t i = 0; i < res.size(); i++){
		ASSERT_NEAR(res[i], real_res[i], pow(10, -15));
	}
}

TEST(Dense, left_multiply_by_vector){
	std::vector<std::vector<double>> a{{1, 2, 0, 3, 5}, {4, 0, 6, 0, 13}, {7, 0, 0, 8, 6}, {0, 3, 0, 0, 9}};
        std::vector<double> b{5, 7, 0, 3};
        DenseMatrix A(a);
        std::vector<double> res = b * A;
        std::vector<double> real_res{33, 19, 42, 15, 143};
        for(std::size_t i = 0; i < res.size(); i++){
                ASSERT_NEAR(res[i], real_res[i], pow(10, -15));
        }
}

TEST(Dense, add){
	std::vector<double> a{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        std::vector<double> b{13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
	std::vector<double> c{14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36};
	DenseMatrix A(a, 3, 4);
	DenseMatrix B(b, 3, 4);
	DenseMatrix C = A + B;
	for(std::size_t i = 0; i < c.size(); i++){
                ASSERT_NEAR(C.get_values()[i], c[i], pow(10, -15));
        }
}

TEST(Dense, product){
	std::vector<double> a{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        std::vector<double> b{13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
	std::vector<double> c{190, 200, 210, 470, 496, 522, 750, 792, 834};
	DenseMatrix A(a, 3, 4);
        DenseMatrix B(b, 4, 3);
	DenseMatrix C = A * B;
	for(std::size_t i = 0; i < c.size(); i++){
                ASSERT_NEAR(C.get_values()[i], c[i], pow(10, -13));
        }
}
