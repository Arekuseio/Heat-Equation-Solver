
#include "Function.h"

Function::Function(const char* a) {
	Parser p(a);
	func = p.parse();
}

double Function::calc(double x, double t) const {
	return eval(func, x, t);
}