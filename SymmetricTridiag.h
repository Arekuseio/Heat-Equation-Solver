#pragma once

#include <algorithm>

class SymmetricTridiag {
	size_t n = 0;
	double diag_elem = 0.f;
	double side_elem = 0.f;
	
	struct Line {
		size_t line;
		double diag_elem = 0.f;
		double side_elem = 0.f;
		
		double operator[](size_t column) const;
	};

public:
	SymmetricTridiag(double diag, double side, size_t size);
	
	size_t size() const;

	Line operator[](size_t line) const;
};