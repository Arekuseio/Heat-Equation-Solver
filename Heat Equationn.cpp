﻿#include "Heat Equationn.h"

inline double u(double x, double t) {
	return -0.5 * pow(x, 4.f) + pow(x, 2.f) - x + t * x + 2.f * pow(t, 5.f) - t * pow(M_E, x);
}

double Solver::half_f(double x, double t, double thau, double h) const {
	return f.calc(x, t) + h * h / 12.f * d2f.calc(x, t + thau / 2.f);
}

Solver::Solver(const char* func, const char* d2func, double a,
			   const char* mu, const char* mu1, const char* mu2) 
			: f(func), d2f(d2func), a(a), mu(mu), mu1(mu1), mu2(mu2) { }

void Solver::addParams(int N_x, int N_t) {
	params.push_back({ N_x, N_t });
}

void Solver::solveAll() {
	for (auto& p : params) {
		Solve(p.first, p.second);
	}
}

void Solver::printResult(std::ofstream& out, int i) const {
	auto& A = answers[i].first;

	int N_x = params[i].first;
	int N_t = params[i].second;
	const double thau = 1.f / N_t;
	const double h = 1.f / N_x;


	out << std::setfill('_') << std::setw(10) << 0;
	for (size_t i = 1; i < N_t + 1; ++i) {
		out << std::setw(11) << i * thau;
	}
	out << '\n';
	out << std::setprecision(4);
	out << std::setfill(' ');
	for (size_t i = 0; i < A.size(); ++i) {
		out << std::setw(4) << i * h << "|";
		for (size_t j = 0; j < A[i].size(); ++j) {
			out << std::setw(10) << A[i][j] << " ";
		}
		out << '\n';
	}
	out << "Time: " << answers[i].second << '\n';
	out << "Error: " << errors[i] << '\n';
}

void Solver::Solve(int N_x, int N_t) {
	std::cout << "Solving for N_x = " << N_x << " , N_t = " << N_t << '\n';
	auto start = std::chrono::steady_clock::now();

	const double thau = 1.f / N_t;
	const double h = 1.f / N_x;
	const double sigma = 0.5 - h * h / (12.f * a * thau);

	// Initializing matrix for linear operator
	std::vector<std::vector<double>> A2;
	initA(A2, N_x, thau, h, sigma);

	// Initializing answer matrix with boundary conditions
	std::vector<std::vector<double>> answer(N_x + 1);
	for (auto& ans : answer) {
		ans.resize(N_t + 1);
	}

	for (size_t j = 0; j < answer.size(); ++j) {
		answer[j][0] = mu.calc(j * h, 0);
	}

	for (size_t j = 0; j < answer[0].size(); ++j) {
		answer[0][j] = mu1.calc(j * thau, 0);
		answer[N_x][j] = mu2.calc(j * thau, 0);
	}
	//----------------------------------

	std::vector<double> u_0(N_x - 1);
	std::vector<double> v_j(N_x - 1, 0.f);
	for (size_t j = 1; j < N_t + 1; ++j) {

		for (size_t i = 0; i < v_j.size(); ++i) {
			v_j[i] = ((answer[i + 2][j - 1] - 2.f * answer[i + 1][j - 1] + answer[i][j - 1]) / (h * h)) * (1.f - sigma) * a;
			v_j[i] = v_j[i] + answer[i + 1][j - 1] / thau + half_f(h * (i + 1), j * thau, thau, h);
		}

		// Adding boundary conditions
		v_j[0] -= -a * sigma / (h * h) * answer[0][j];
		v_j[v_j.size() - 1] -= -a * sigma / (h * h) * answer[N_x][j];
		//-----------------------------------------------------------

		u_0 = SolveTridiag(A2, v_j);

		for (size_t k = 0; k < u_0.size(); ++k) {
			answer[k + 1][j] = u_0[k];
		}
	}

	auto end = std::chrono::steady_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Done! Time: " << time.count() << '\n';

	answers.push_back({ answer, time });

	double error = Compare(answer, N_x, N_t, thau, h);
	errors.push_back(error);
}


std::vector<double> SolveTridiag(std::vector<std::vector<double>>& A, std::vector<double>& f) {

	if (A.size() != f.size()) throw std::runtime_error("Mismatch of dimensions");

	std::vector<double> u(A.size());

	std::vector<double> a(A.size() - 1);
	std::vector<double> b(A.size());

	a[0] = A[0][1] / A[0][0];
	b[0] = f[0] / A[0][0];

	for (size_t i = 1; i < A.size(); ++i) {
		if (i != A.size() - 1) {
			a[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * a[i - 1]);
		}
		b[i] = (f[i] - A[i][i - 1] * b[i - 1]) / (A[i][i] - A[i][i - 1] * a[i - 1]);
	}

	u[A.size() - 1] = b[A.size() - 1];
	for (int32_t i = A.size() - 2; i >= 0; --i) {
		u[i] = b[i] - a[i] * u[i + 1];
	}

	return u;
}

void Solver::initA(std::vector<std::vector<double>>& A, int N_x, 
					double thau, double h, double sigma) const {
	A.resize(N_x - 1);

	for (int i = 0; i < N_x - 1; ++i) {
		A[i].resize(N_x - 1, 0.f);

		for (int j = 0; j < N_x - 1; ++j) {
			if (j == i) {
				A[i][j] = (1.f / thau + 2.f * sigma * a / (h * h));
			}
			if (abs(i - j) == 1) {
				A[i][j] = -a * sigma / (h * h);
			}
		}
	}
}

double Compare(std::vector<std::vector<double>> myanswer, 
			int N_x, int N_t, double thau, double h) {

	double max_diff = 0.f;
	double diff;
	std::vector<std::vector<double>> answer(N_x + 1);
	for (size_t i = 0; i < answer.size(); ++i) {
		answer[i].resize(N_t + 1);
		for (size_t k = 0; k < N_t + 1; ++k) {
			answer[i][k] = u(h * i, thau * k);
			diff = fabs(answer[i][k] - myanswer[i][k]);
			if (diff > max_diff) {
				max_diff = diff;
			}
		}
	}
	return max_diff;
}


int main() {
	return 0;
}