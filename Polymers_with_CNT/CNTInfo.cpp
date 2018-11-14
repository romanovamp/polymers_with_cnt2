#include "CNTInfo.h"
#include "CNT.h"
#include <cmath>
#define M_PI 3.14159265358979323846
CNTInfo::CNTInfo()
{
	x1 = 0;
	y1 = 0;
	x2 = 0;
	y2 = 0;
	x3 = 0;
	y3 = 0;
	x4 = 0;
	y4 = 0;

}

CNTInfo::CNTInfo(double _x1, double _y1, double _x2, double _y2, double _x3, double _y3, double _x4, double _y4)
{
	x1 = _x1;
	y1 = _y1;
	x2 = _x2;
	y2 = _y2;
	x3 = _x3;
	y3 = _y3;
	x4 = _x4;
	y4 = _y4;
}
double coordX(double x, double k, double a)
{
	return x + k * cos((a*M_PI) / 180.0);
}
double coordY(double y, double k, double a)
{
	return y + k * sin((a*M_PI) / 180.0);
}
CNTInfo::CNTInfo(double _x, double _y, double _k, int _a, double r)
{
	x1 = coordX(_x, r, _a + 90);
	y1 = coordY(_y, r, _a + 90);

	x2 = coordX(_x, r, _a - 90);
	y2 = coordY(_y, r, _a - 90);

	x3 = coordX(x2, _k, _a);
	y3 = coordY(y2, _k, _a);

	x4 = coordX(x1, _k, _a);
	y4 = coordY(y1, _k, _a);

	a = _a;
	k = _k;
	x = _x;
	y = _y;

}
void CNTInfo::top(double &a, double &b, double &c)
{
	a = y3 - y4;
	b = x4 - x3;
	c = x3 * y4 - x4 * y3;
}
void CNTInfo::bottom(double &a, double &b, double &c)
{
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}
void CNTInfo::left(double &a, double &b, double &c)
{
	a = y1 - y4;
	b = x4 - x1;
	c = x1 * y4 - x4 * y1;
}
void CNTInfo::right(double &a, double &b, double &c)
{
	a = y3 - y2;
	b = x2 - x3;
	c = x3 * y2 - x2 * y3;
}

double CNTInfo::tA()
{
	return y3 - y4;
}
double CNTInfo::tB()
{
	return x4 - x3;
}
double CNTInfo::tC()
{
	return x3 * y4 - x4 * y3;
}

double CNTInfo::bA()
{
	return y1 - y2;
}
double CNTInfo::bB()
{
	return x2 - x1;
}
double CNTInfo::bC()
{
	return x1 * y2 - x2 * y1;
}

double CNTInfo::lA()
{
	return y1 - y4;
}
double CNTInfo::lB()
{
	return x4 - x1;
}
double CNTInfo::lC()
{
	return x1 * y4 - x4 * y1;
}

double CNTInfo::rA()
{
	return y3 - y2;
}
double CNTInfo::rB()
{
	return x2 - x3;
}
double CNTInfo::rC()
{
	return x3 * y2 - x2 * y3;
}
