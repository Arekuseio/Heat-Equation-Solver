#include "SymmetricTridiag.h"
#include <algorithm>


size_t SymmetricTridiag::size() const {
	return n;
}

SymmetricTridiag::SymmetricTridiag(double diag, double side, size_t size) {
	n = size;
	diag_elem = diag;
	side_elem = side;
}

SymmetricTridiag::Line SymmetricTridiag::operator[](size_t line) const {

	return Line{ line, diag_elem, side_elem };
}

double SymmetricTridiag::Line::operator[](size_t column) const {
	if (line == column) {
		return diag_elem;
	}
	else {
		if (std::max(line, column) - std::min(line, column) == 1) {
			return side_elem;
		}
		else {
			return 0.f;
		}
	}
}