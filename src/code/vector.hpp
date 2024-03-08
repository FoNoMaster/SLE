#pragma once

#include <vector>
#include <iostream>
#include <cmath>

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> res(a.size());
        for(std::size_t i = 0; i < a.size(); i++)
                res[i] = (a[i] + b[i]);
        return res;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> res(a.size());
        for(std::size_t i = 0; i < a.size(); i++)
                res[i] = (a[i] - b[i]);
        return res;

}

template<typename T>
std::vector<T> operator*(const T& a, const std::vector<T>& b){
        std::vector<T> res(b.size());
        for(std::size_t i = 0; i < b.size(); i++)
                res[i] = (a * b[i]);
        return res;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& a, const T& b){
        std::vector<T> res(a.size());
        for(std::size_t i = 0; i < a.size(); i++)
                res[i] = (a[i] * b);
        return res;
}

template<typename T>
std::vector<T> operator/(const std::vector<T>& a, const T& b){
    std::vector<T> res(a.size());
    for(std::size_t i = 0; i < a.size(); i++)
        res[i] = (a[i] / b);
    return res;
}

template<typename T>
T dot(const std::vector<T>& a, const std::vector<T>& b){
        T res = 0;
        for(std::size_t i = 0; i < a.size(); i++)
                res += a[i] * b[i];
        return res;
}

template<typename T>
T norm(const std::vector<T>& v){
    return sqrt(dot(v, v));
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& b){
        std::cout << "[";
        for(std::size_t i = 0; i < b.size() - 1; i++){
                std::cout << b[i] << ", ";
        }
        std::cout << b[b.size() - 1] << "]" << std::endl;
        return os;
}

template<typename T>
T Max_element(const std::vector<T>& v){
    T max = v[0];
    for(std::size_t i = 0; i < v.size(); i++){
        if (v[i] > max)
            max = v[i];
    }
    return max;
}

template<typename T>
std::vector<T> abs(const std::vector<T>& a){
    std::vector<T> res(a.size());
    for(std::size_t i = 0; i < a.size(); i++){
        if(a[i] < 0)
            res[i] = -a[i];
        else
            res[i] = a[i];
    }
    return res;
}
