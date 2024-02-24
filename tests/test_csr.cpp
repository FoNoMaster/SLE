#include <gtest/gtest.h>
#include <code/CSR.hpp>
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
