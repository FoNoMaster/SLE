#pragma once

#include "CSR.hpp"

template<typename T>
T max_lamda(const CSR_Matrix<T>& A){
    std::vector<T> v(Max_element(A.get_cols()), 1);
    std::vector<T> mu(2);
    mu[0] = 5;
    mu[1] = 100;
    while (mu[1] - mu[0] > 1e-15){
        v = A * v / norm(A * v);
        mu[0] = mu[1];
        mu[1] = dot(v, A * v) / dot(v, v);
    }
    return mu[1];
}


template<typename T>
std::vector<T> Simple_Iteration_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter, const T& lambda_max){
    double tau = 2 * 0.9 / (lambda_max);
    std::vector<T> x = x0;
    std::vector<T> r(x0.size());
    for (std::size_t i = 0; i < Niter; ++i){
        r = A * x - b;
        if (norm(r) < percision){
            break;
        }
        x = x - tau * r;
    }
    return x;
}


template<typename T>
std::vector<T> Simple_Iteration_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter){
    double tau = 1 / max_lamda(A);
    std::vector<T> x = x0;
    std::vector<T> r(x0.size());
    for (std::size_t i = 0; i < Niter; ++i){
        r = A * x - b;
        if (norm(r) < percision){
            break;
        }
        x = x - tau * r;
    }
    return x;
}

template<typename T>
std::vector<T> Jacobi_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter){
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
        if(norm(A * x - b) < percision)
            break;
    }

    return x;
}

template<typename T>
std::vector<T> Gauss_Sejdel_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter){
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

        if(norm(A * x - b) < percision)
            break;
    }
    return x;
}


std::vector<std::size_t> redecorate(std::size_t r){
    std::size_t n = static_cast<std::size_t>(std::pow(2, r));
    std::vector<std::size_t> index(n, 0);
    index[static_cast<std::size_t>(std::pow(2, r - 1))] = 1;
    std::size_t step = 0;

    for (std::size_t j = 2; j <= r; j++) {
        step = static_cast<std::size_t>(std::pow(2, r - j));
        for (std::size_t k = 0; k < n; k += 2 * step) {
            index[k + step] = (static_cast<std::size_t>(std::pow(2, j))) - index[k] - 1;
        }
    }

    return index;
}

std::vector<double> find_tau(const std::size_t r, const double lambda_max, const double lambda_min){
    std::size_t n = static_cast<size_t>(std::pow(2, r));
    std::vector<double> roots(n, 0);
    const double cos_n = std::cos(3.141592653589793238462643383279 / static_cast<int>(n));
    const double sin_n = std::sin(3.141592653589793238462643383279 / static_cast<int>(n));
    const double cos_2n = std::cos(3.141592653589793238462643383279 / (2 * static_cast<int>(n)));
    double sin_i = std::sin(3.141592653589793238462643383279 / (2 * static_cast<int>(n)));
    roots[0] = cos_2n;
    for(std::size_t i = 1; i < n / 2 + 1; i++){
        roots[i] = roots[i - 1] * cos_n - sin_i * sin_n;
        sin_i = sin_i * cos_n + roots[i - 1] * sin_n;
        roots[n - i] = -roots[i - 1];
        roots[i - 1] = (lambda_min + lambda_max) / 2 + ((lambda_max - lambda_min) / 2) * roots[i - 1];
        roots[n - i] = (lambda_min + lambda_max) / 2 + ((lambda_max - lambda_min) / 2) * roots[n - i];
    }
    for (auto& it : roots){
        it = 1 / it;
    }
    return roots;
}


template <typename T>
std::vector<double> chebyshev(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const double lambda_max, const double lambda_min, size_t r_, const double percision){
    std::vector<double> x = x0;
    std::vector<double> r(b.size());
    r =  A * x - b;
    double r_norm = norm(r);

    std::vector<size_t> index = redecorate(r_);
    std::vector<double> tau = find_tau(r_, lambda_max, lambda_min);

    while (r_norm > percision){
        for (auto& it : index) {
            r = A * x - b;
            x = x - tau[it] * r;
            r_norm = norm(r);
        }
    }

    return x;
}

template <typename T>
std::vector<double> chebyshevv(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const double lambda_max, const double lambda_min, size_t r_, const std::size_t& Niter){
    std::vector<double> x = x0;
    std::vector<double> r(b.size());

    std::vector<size_t> index = redecorate(r_);
    std::vector<double> tau = find_tau(r_, lambda_max, lambda_min);

    for(std::size_t i = 0; i < Niter; i++){
        for (auto& it : index) {
            r = A * x - b;
            x = x - tau[it] * r;
        }
    }

    return x;
}

template<typename T>
std::vector<double> FGD(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const double percision, const std::size_t &Niter){
    std::vector<double> x = x0;
    std::size_t n = 0;
    std::vector<double> r(b.size());

    r = A * x - b;
    double tau;

    while (norm(r) > percision){
        r = A * x - b;
        tau = dot(r, r) / dot(r, (A * r));
        x = x - tau * r;
        n++;
        if(n > Niter)
            break;
    }

    return x;
}

template<typename T>
std::vector<double> FGD(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const double percision){
    std::vector<double> x = x0;
    std::vector<double> r(b.size());

    r = A * x - b;
    double tau;

    while (norm(r) > percision){
        r = A * x - b;
        tau = dot(r, r) / dot(r, (A * r));
        x = x - tau * r;
    }

    return x;
}


template<typename T>
std::vector<T> Sym_Gauss_Sejdel_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter){
    std::vector<T> x = x0;
    T tmp;

    for (std::size_t j = 0; j < Niter; j++) {
        for(std::size_t i = 0; i < x.size(); i++){
            x[i] = b[i];
            for(std::size_t k = A.r(i); k < A.r(i + 1); k++){
                if(i != A.c(k))
                    x[i] -= A.v(k) * x[A.c(k)];
                else
                    tmp = A.v(k);
            }
            x[i] /= tmp;
        }

        for (std::size_t i = x.size() - 1; i > 0; i--) {
            x[i] = b[i];
            for (std::size_t k = A.r(i); k < A.r(i + 1); k++) {
                if (i != A.c(k))
                    x[i] -= A.v(k) * x[A.c(k)];
                else
                    tmp = A.v(k);
            }
            x[i] /= tmp;
        }

        if(norm(A * x - b) < percision)
            break;
    }
    return x;
}


template<typename T>
std::vector<T> chebyshevSGS(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter) {
    T rho = 0.9;
    std::vector<T> mu = {1, 1 / rho, 0};
    std::vector<std::vector<T>> x(3, x0);

    T tmp;

    for (std::size_t j = 0; j < Niter; j++) {

        mu[2] = 2 * mu[1] / rho - mu[0];

        for (std::size_t i = 0; i < x[1].size(); i++) {
            x[1][i] = b[i];
            for (std::size_t k = A.r(i); k < A.r(i + 1); k++) {
                if (i != A.c(k))
                    x[1][i] -= A.v(k) * x[1][A.c(k)];
                else
                    tmp = A.v(k);
            }
            x[1][i] /= tmp;
        }

        for (std::size_t i = x[1].size() - 1; i > 0; i--) {
            x[1][i] = b[i];
            for (std::size_t k = A.r(i); k < A.r(i + 1); k++) {
                if (i != A.c(k))
                    x[1][i] -= A.v(k) * x[1][A.c(k)];
                else
                    tmp = A.v(k);
            }
            x[1][i] /= tmp;
        }

        x[2] = (2 * mu[1])/(rho * mu[2]) * x[1] - (mu[0] / mu[2]) * x[0];
        mu[0] = mu[1];
        mu[1] = mu[2];
        x[0] = x[1];
        x[1] = x[2];

        if (norm(A * x[1] - b) < percision)
            break;
    }
    return x[2];
}


template<typename T>
std::vector<T> Conjugate_Gradient_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0) {
    std::vector<std::vector<T>> x(2, x0);
    std::vector<std::vector<T>> r(2, A * x[0] - b);
    std::vector<T> d = r[0];

    while(norm(r[1]) != 0){
        x[1] = x[0] - (dot(r[0], r[0]))/(dot((A * d), d)) * d;
        r[1] = A * x[1] - b;
        d = r[1] + dot(r[1], r[1]) / dot(r[0], r[0]) * d;

        x[0] = x[1];
        r[0] = r[1];
    }

    return x[1];
}


template<typename T>
std::vector<T> Heavy_Ball_Method(const CSR_Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T& percision, const std::size_t& Niter, const T& betta){
    double tau = 1 / max_lamda(A);
    std::vector<T> x = x0;
    std::vector<T> x1 = x0;
    std::vector<T> x2 = x0;
    std::vector<T> r(x0.size());
    for (std::size_t i = 0; i < Niter; ++i){
        r = A * x1 - b;
        if (norm(r) < percision){
            break;
        }
        x = x1 - tau * r + betta * (x1 - x2);
        x2 = x1;
        x1 = x;
    }
    return x;
}

