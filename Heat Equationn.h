#pragma once

#define _USE_MATH_DEFINES

#include "Function.h"

#include <iostream>
#include <chrono>
#include <future>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <utility>

using _answer = std::pair<std::vector<std::vector<double>>, std::chrono::milliseconds>;
using _parameteres = std::pair<size_t, size_t>;

// Solving Au=f for tridiagonal matrix A
std::vector<double> SolveTridiag(std::vector<std::vector<double>>& A, std::vector<double>& f);


// Solving tridiagonal matrix system using two threads
std::vector<double> SolveTridiagDoubleThread(std::vector<std::vector<double>>& A, std::vector<double>& f);

std::pair<std::vector<double>, std::vector<double>> CalcCoeffsRight(std::vector<std::vector<double>>& A, size_t p, std::vector<double>& f);

std::pair <std::vector<double>, std::vector<double>> CalcCoeffsLeft(std::vector<std::vector<double>>& A, size_t p, std::vector<double>& f);

// Real solution function
inline double u(double x, double t);

// Comparing the real solution with approximation, return the max difference
double Compare(std::vector<std::vector<double>> myanswer,
	size_t N_x, size_t N_t, double thau, double h);

class Solver {
	Function f;
	Function d2f;

	// boundary conditions: mu = u(x, 0), mu1 = u(0,t), mu2 = u(1,t)
	Function mu, mu1, mu2;
	const double a;

	std::vector<_answer> answers;
	std::vector<_parameteres> params;
	std::vector<double> errors;
public:
	Solver(const char* func, const char* d2func, double a,
		const char* mu, const char* mu1, const char* mu2);

	void addParams(size_t N_x, size_t N_t);
	void solveAll();
	void printResult(std::ostream& out, size_t i) const;
	size_t getAnswersSize() const;
private:
	// f + h^2/12 * f''(x, t + thau/2)
	double half_f(double x, double t, double thau, double h) const;

	void Solve(size_t N_x, size_t N_t);

	// Initializing matrix of linear operator
	void initA(std::vector<std::vector<double>>& A, size_t N_x, double thau, double h, double sigma) const;
};
