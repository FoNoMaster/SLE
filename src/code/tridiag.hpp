#pragma once

#include "vector.hpp"

template<typename T>
class TriDiagMatrix{
private:
	std::vector<T> a_;
	std::vector<T> b_;
	std::vector<T> c_;
public:
	TriDiagMatrix(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c): a_(a), b_(b), c_(c){}

	const std::vector<T>& get_a() const{return a_;}
	const std::vector<T>& get_b() const{return b_;}
	const std::vector<T>& get_c() const{return c_;}

	const T& a(const std::size_t& i) const{return a_[i];}
	const T& b(const std::size_t& i) const{return b_[i];}
	const T& c(const std::size_t& i) const{return c_[i];}
};


template<typename T>
std::vector<T> TriDiagSolve(const TriDiagMatrix<T>& A, const std::vector<T>& d){
	std::vector<T> p(d.size());
	std::vector<T> q(d.size());
	std::vector<T> x(d.size());

	p[0] = -A.c(0) / A.b(0);
	q[0] = d[0] / A.b(0);

	for(std::size_t i = 1; i < d.size() - 1; i++){
		p[i] = -(A.c(i) / (A.a(i - 1) * p[i - 1] + A.b(i)));
		q[i] = (d[i] - A.a(i - 1) * q[i - 1]) / (A.a(i - 1) * p[i - 1] + A.b(i));
	}

	q[d.size() - 1] = (d[d.size() - 1] - A.a(d.size() - 2) * q[d.size() - 2]) / (A.a(d.size() - 2) * p[d.size() - 2] + A.b(d.size() - 1));

	x[d.size() - 1] = q[d.size() - 1];

	for(std::size_t i = 1; i < d.size(); i++){
		x[d.size() - 1 - i] = p[d.size() - 1 - i] * x[d.size() - i] + q[d.size() - 1 - i];
	}
	
	return x;
}

