#include <gtest/gtest.h>
#include <code/solvers.hpp>
#include <cmath>


TEST(CSR, right_product){
        std::vector<std::vector<double>> a{{1, 2}, {3, 4}};
        CSR_Matrix A(a);
        A = A * 2.0;
        std::vector<double> b{2, 4, 6, 8};
        for(std::size_t i = 0; i < b.size(); i++){
                ASSERT_NEAR(A.get_vals()[i], b[i], pow(10, -15));
        }
}


TEST(CSR, left_product){
        std::vector<std::vector<double>> a{{1, 2}, {3, 4}};
        CSR_Matrix A(a);
        A = 2.0 * A;
        std::vector<double> b{2, 4, 6, 8};
        for(std::size_t i = 0; i < b.size(); i++){
                ASSERT_NEAR(A.get_vals()[i], b[i], pow(10, -15));
        }
}


TEST(CSR, right_multiply_by_vector){
        std::vector<std::vector<double>> a{{1, 2, 0, 3, 5}, {4, 0, 6, 0, 13}, {7, 0, 0, 8, 6}, {0, 3, 0, 0, 9}};
        std::vector<double> b{5, 7, 0, 3, 0};
        CSR_Matrix A(a);
        std::vector<double> res = A * b;
        std::vector<double> real_res{28, 20, 59, 21};
        for(std::size_t i = 0; i < res.size(); i++){
                ASSERT_NEAR(res[i], real_res[i], pow(10, -15));
        }
}

TEST(CSR, Simple_Iteration_Method){
	std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    CSR_Matrix<double> A(vals, 3, 3);

	std::vector<double> x0 = {4, 4, 4};
	std::vector<double> b = {5, 17, 32};
    double tol = 1e-20;
	std::vector<double> x = Simple_Iteration_Method(A, b, x0, tol, 10000);
	std::vector<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        ASSERT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR, Jacobi_Method){
	std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
        CSR_Matrix<double> A(vals, 3, 3);


        std::vector<double> x0 = {4, 4, 4};
        std::vector<double> b = {5, 17, 32};
        double tol = 1e-20;
        std::vector<double> x = Jacobi_Method(A, b, x0, tol, 10000);
        std::vector<double> expected = {1, 2, 3};
        for (std::size_t j = 0; j < 3; ++j)
                ASSERT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR, Gauss_Sejdel_Method){
        std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
        CSR_Matrix<double> A(vals, 3, 3);

        std::vector<double> x0 = {4, 4, 4};
        std::vector<double> b = {5, 17, 32};
        double tol = 1e-20;
        std::vector<double> x = Gauss_Sejdel_Method(A, b, x0, tol, 10000);
        std::vector<double> expected = {1, 2, 3};
        for (std::size_t j = 0; j < 3; ++j)
                ASSERT_NEAR(expected[j], x[j], 0.01);
}

TEST(Chebyshev, Test1) {
    std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    CSR_Matrix<double> A(vals, 3, 3);

    std::vector<double> x0 = {4, 4, 4};
    std::vector<double> b = {5, 17, 32};
    double tol = 1e-20;

    double lambda_min = 0.287;
    double lambda_max = 10.261;
    std::vector<double> x = chebyshev(A, b, x0, lambda_max, lambda_min, 7, tol);
    std::vector<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        ASSERT_NEAR(expected[j], x[j], 0.01);
}


TEST(FGD, Test1)
{
    std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    CSR_Matrix<double> A(vals, 3, 3);

    std::vector<double> x0 = {4, 4, 4};
    std::vector<double> b = {5, 17, 32};
    double tol = 1e-10;
    std::vector<double> x = FGD(A, b, x0, tol);
    std::vector<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        ASSERT_NEAR(expected[j], x[j], 0.01);
}

TEST(Sym_Gauss_Zeidel_Method, Test1)
{
    std::vector<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    CSR_Matrix<double> A(vals, 3, 3);

    std::vector<double> x0 = {4, 4, 4};
    std::vector<double> b = {5, 17, 32};
    double tol = 1e-20;
    std::vector<double> x = Sym_Gauss_Sejdel_Method(A, b, x0, tol, 10000);
    std::vector<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j) {
        ASSERT_NEAR(expected[j], x[j], 0.01);
    }
}
