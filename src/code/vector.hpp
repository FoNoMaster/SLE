#pragma once

#include <vector>
#include <iostream>

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> res;
        for(std::size_t i = 0; i < a.size(); i++)
                res.push_back(a[i] + b[i]);
        return res;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> res;
        for(std::size_t i = 0; i < a.size(); i++)
                res.push_back(a[i] - b[i]);
        return res;

}

template<typename T>
std::vector<T> operator*(const T& a, const std::vector<T>& b){
        std::vector<T> res;
        for(std::size_t i = 0; i < b.size(); i++)
                res.push_back(a * b[i]);
        return res;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& a, const T& b){
        std::vector<T> res;
        for(std::size_t i = 0; i < a.size(); i++)
                res.push_back(a[i] * b);
        return res;
}

template<typename T>
T dot(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> res;
        for(std::size_t i = 0; i < a.size(); i++)
                res.push_back(a[i] * b[i]);
        return res;
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
