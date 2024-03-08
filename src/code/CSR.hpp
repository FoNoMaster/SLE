#pragma once

#include "vector.hpp"


template<typename T>
class CSR_Matrix{
private:
	std::vector<T> vals_;
	std::vector<std::size_t> cols_;
	std::vector<std::size_t> rows_;
public:
	CSR_Matrix(const std::vector<T>& vals, const std::vector<std::size_t>& cols, const std::vector<std::size_t>& rows);
	CSR_Matrix(const std::vector<std::vector<T>>& data);
	CSR_Matrix (const std::vector<T>& data, std::size_t M, std::size_t N);

	const std::vector<T>& get_vals() const{return vals_;}
	const std::vector<std::size_t>& get_cols() const{return cols_;}
	const std::vector<std::size_t>& get_rows() const{return rows_;}
	const T& v(const std::size_t& i) const{return vals_[i];}
	const std::size_t& c(const std::size_t& i) const{return cols_[i];}
	const std::size_t& r(const std::size_t& i) const{return rows_[i];}

	T operator()(std::size_t i, std::size_t j) const;

	std::vector<T> operator*(const std::vector<T>& b) const;

};

template<typename T>
CSR_Matrix<T>::CSR_Matrix(const std::vector<T>& vals, const std::vector<std::size_t>& cols, const std::vector<std::size_t>& rows): vals_(vals), cols_(cols), rows_(rows){}

template<typename T>
CSR_Matrix<T>::CSR_Matrix(const std::vector<std::vector<T>>& data){
        std::size_t counter = 0;
        rows_.push_back(counter);
        for(std::size_t i = 0; i < data.size(); i++){
                for(std::size_t j = 0; j < data[0].size(); j++){
                        if (data[i][j] != 0) {
                                vals_.push_back(data[i][j]);
                                cols_.push_back(j);
                                counter++;
                        }
                }
                rows_.push_back(counter);
        }
}

template<typename T>
CSR_Matrix<T>::CSR_Matrix (const std::vector<T>& data, std::size_t M, std::size_t N){
        std::size_t counter = 0;
        rows_.push_back(counter);
        for(std::size_t i = 0; i < M; i++){
                for(std::size_t j = 0; j < N; j++){
                        if (data[i * M + j] != 0) {
                                vals_.push_back(data[i * M + j]);
                                cols_.push_back(j);
                                counter++;
                        }
                }
                rows_.push_back(counter);
        }
}


template<typename T>
T CSR_Matrix<T>::operator()(std::size_t i, std::size_t j) const{
        for (std::size_t k = rows_[i]; k < rows_[i + 1]; k++) {
                if (cols_[k] == j) {
                        return vals_[k];
                }
        }
        return 0;
}

template<typename T>
std::vector<T> CSR_Matrix<T>::operator*(const std::vector<T>& b) const{
        std::vector<T> res(rows_.size() - 1, 0);
        for(std::size_t i = 0; i < b.size(); i++){
                for(std::size_t k = rows_[i]; k < rows_[i + 1]; k++){
                        res[i] += vals_[k] * b[cols_[k]];
                }
        }
        return res;
}

template<typename T>
CSR_Matrix<T> operator*(const CSR_Matrix<T>& A, const T& b){
	std::vector<T> res1 = A.get_vals() * b;
        std::vector<std::size_t> res2 = A.get_cols();
        std::vector<std::size_t> res3 = A.get_rows();
        CSR_Matrix<T> RES(res1, res2, res3);
        return RES;
}


template<typename T>
CSR_Matrix<T> operator*(const T& b, const CSR_Matrix<T>& A){
        std::vector<T> res1 = A.get_vals() * b;
        std::vector<std::size_t> res2 = A.get_cols();
        std::vector<std::size_t> res3 = A.get_rows();
        CSR_Matrix<T> RES(res1, res2, res3);
        return RES;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const CSR_Matrix<T>& A){
        std::cout << "vals: ";
        std::cout << A.get_vals();
        std::cout << "cols: ";
        std::cout << A.get_cols();
        std::cout << "rows: ";
        std::cout << A.get_rows();
        return os;
}

template<typename T>
T max_lamda(const CSR_Matrix<T>& A){
        std::vector<T> v(Max_element(A.get_cols()), 1);
	std::vector<T> mu(2);
	mu[0] = 5;
    	mu[1] = 100;
    	while (mu[1] - mu[0] > 1e-5){
        	v = A * v / norm(A * v);
        	mu[0] = mu[1];
        	mu[1] = dot(v, A * v) / dot(v, v);
    	}
    	return mu[1];
}

template<typename T>
std::vector<T> Simple_Iteration_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& tol, const std::size_t& Niter){
    double tau = 1 / max_lamda(A);
    std::vector<T> x = x0;
    std::vector<T> r(x0.size());
    for (std::size_t i = 0; i < Niter; ++i){
        r = A * x - b;
        if (norm(r) < tol){
            break;
        }
        x = x - tau * r;
    }
    return x;
}

template<typename T>
std::vector<T> Jacobi_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& tol, const std::size_t& Niter){
    std::vector<T> x = x0;
    std::vector<T> res(x.size());
    T tmp;
    for (std::size_t j = 0; j < Niter; j++) {
        res = b;
        for (std::size_t i = 0; i < x.size(); i++) {
            for (std::size_t k = A.get_rows()[i]; k < A.get_rows()[i + 1]; k++) {
                if(i != A.c(k))
                    res[i] -= A.get_vals()[k] * x[A.get_cols()[k]];
                else
                    tmp = A.v(k);
            }
            res[i] /= tmp;
        }

        x = res;
        if(Max_element(abs(A * x - b)) < tol)
            break;
    }

    return x;
}

template<typename T>
std::vector<T> Gauss_Sejdel_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& tol, const std::size_t& Niter){
    std::vector<T> x = x0;
    T tmp;
    
    for (std::size_t j = 0; j < Niter; j++) {
        for(std::size_t i = 0; i < x.size(); i++){
            x[i] = b[i];
            for(std::size_t k = A.r(i); k < A.r(i + 1); k++){
                if(i != A.c(k))
                    x[i] -= A.get_vals()[k] * x[A.get_cols()[k]];
                else
                    tmp = A.v(k);
            }
            x[i] /= tmp;
        }

        if(Max_element(abs(A * x - b)) < tol)
            break;
    }
    return x;
}
