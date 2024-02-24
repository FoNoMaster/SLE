#pragma once

#include "vector.hpp"

template<typename T>
class DenseMatrix{
private:
	std::vector<T> data_;
	std::size_t M, N;
public:
	DenseMatrix(const std::vector<T>& data, const std::size_t m, const std::size_t n):
		data_(data), M(m), N(n){}
	DenseMatrix(const std::vector<std::vector<T>>& data);

	std::size_t get_size_M() const{return M;}
	std::size_t get_size_N() const{return N;}

	const std::vector<T>& get_values() const{return data_;}

	T operator()(std::size_t i, std::size_t j) const {return data_[i * N + j];}
};


template<typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<std::vector<T>>& data){
        M = data.size();
        N = data[0].size();
        for(std::size_t i = 0; i < M; i++){
                for(std::size_t j = 0; j < N; j++){
                        data_.push_back(data[i][j]);
		}
        }
}


template<typename T>
std::vector<T> operator*(const DenseMatrix<T>& A, const std::vector<T>& b){
        std::vector<T> res(A.get_size_M(), 0);
        for(std::size_t i = 0; i < A.get_size_M(); i++){
                for(std::size_t j = 0; j < b.size(); j++){
                        res[i] += b[j] * A(i, j);
                }
        }
        return res;
}


template<typename T>
std::vector<T> operator*(const std::vector<T>& b, const DenseMatrix<T>& A){
        std::vector<T> res(A.get_size_N(), 0);
        for(std::size_t i = 0; i < A.get_size_N(); i++){
                for(std::size_t j = 0; j < A.get_size_M(); j++){
                        res[i] += b[j] * A(j, i);
                }
        }
        return res;
}


template<typename T>
DenseMatrix<T> operator*(const DenseMatrix<T>& A, const T& b){
        std::vector<T> res;
        res = A.get_values() * b;
        DenseMatrix<T> RES(res, A.get_size_M(), A.get_size_N());
        return RES;
}


template<typename T>
DenseMatrix<T> operator*(const T& a, const DenseMatrix<T>& A){
        std::vector<T> res;
        res = A.get_values() * a;
        DenseMatrix<T> RES(res, A.get_size_M(), A.get_size_N());
        return RES;
}


template<typename T>
DenseMatrix<T> operator+(const DenseMatrix<T>& A, const DenseMatrix<T>& B){
	std::vector<T> res(A.get_size_M() * A.get_size_N());
        for(std::size_t i = 0; i < A.get_size_M(); i++){
                for(std::size_t j = 0; j < A.get_size_N(); j++){
                        res[i * A.get_size_N() + j] = (A(i, j) + B(i, j));
                }
        }
        DenseMatrix<T> RES(res, A.get_size_M(), A.get_size_N());
        return RES;
}


template<typename T>
DenseMatrix<T> operator-(const DenseMatrix<T>& A, const DenseMatrix<T>& B){
	DenseMatrix<T> RES = A;
        for(std::size_t i = 0; i < A.get_size_M(); i++){
                for(std::size_t j = 0; j < A.get_size_N(); j++){
                        RES(i, j) =(A(i, j) - B(i, j));
                }
        }
        return RES;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<T>& A){
        std::cout << "[" << std::endl;
        for(std::size_t i = 0; i < A.get_size_M(); i++){
                std::cout << "[";
                for(std::size_t j = 0; j < A.get_size_N() - 1; j++){
                        std::cout << A(i, j) << ", ";
                }
                std::cout << A(i, A.get_size_N() - 1);
                std::cout << "]" << std::endl;
        }
        std::cout << "]" << std::endl;
	return os;
}


template<typename T>
DenseMatrix<T> operator*(const DenseMatrix<T>& A, const DenseMatrix<T>& B){
        std::vector<T> res(A.get_size_M() * B.get_size_N(), 0);
        for(std::size_t i = 0; i < A.get_size_M(); i++){
		for(std::size_t k = 0; k < B.get_size_N(); k++){
                	for(std::size_t j = 0; j < A.get_size_N(); j++){
				res[i * B.get_size_N() + k] += A(i, j) * B(j, k);
			}
		}
        }
        DenseMatrix<T> RES(res, A.get_size_M(), B.get_size_N());
        return RES;
}

