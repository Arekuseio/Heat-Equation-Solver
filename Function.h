#pragma once

#include "Parser.h"

class Function {
	Expression func;
public:
	Function(const char* a);
	double calc(double x, double t) const;
};