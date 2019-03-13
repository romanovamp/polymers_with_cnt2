#include "cnt.h"
#include <cmath>
#define M_PI 3.14159265358979323846

//parent = 0 - не дочерний, не родитель
//parent = 1 - родитель
//parent = 2 - дочерний

cnt::cnt()
{
	x1 = 0;
	y1 = 0;
	x2 = 0;
	y2 = 0;
	x3 = 0;
	y3 = 0;
	x4 = 0;
	y4 = 0;
	idClus = 0;
	idParent = -1;
	parent = 0;
}

double coordX(double x, double k, double a)
{
	return x + k * cos((a*M_PI) / 180.0);
}
double coordY(double y, double k, double a)
{
	return y + k * sin((a*M_PI) / 180.0);
}
cnt::cnt(double _x, double _y, double _k, int _a, double r, int _idParent, int _parent)
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
	idClus = 0;
	idParent = _idParent;
	parent = _parent;
}
cnt::cnt(double _x, double _y, double _k, int _a, double r)
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
	idClus = 0;
	idParent = -1;
	parent = 0;
}
void cnt::top(double &a, double &b, double &c)
{
	a = y3 - y4;
	b = x4 - x3;
	c = x3 * y4 - x4 * y3;
}
void cnt::bottom(double &a, double &b, double &c)
{
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}
void cnt::left(double &a, double &b, double &c)
{
	a = y1 - y4;
	b = x4 - x1;
	c = x1 * y4 - x4 * y1;
}
void cnt::right(double &a, double &b, double &c)
{
	a = y3 - y2;
	b = x2 - x3;
	c = x3 * y2 - x2 * y3;
}

double cnt::tA()
{
	return y3 - y4;
}
double cnt::tB()
{
	return x4 - x3;
}
double cnt::tC()
{
	return x3 * y4 - x4 * y3;
}

double cnt::bA()
{
	return y1 - y2;
}
double cnt::bB()
{
	return x2 - x1;
}
double cnt::bC()
{
	return x1 * y2 - x2 * y1;
}

double cnt::lA()
{
	return y1 - y4;
}
double cnt::lB()
{
	return x4 - x1;
}
double cnt::lC()
{
	return x1 * y4 - x4 * y1;
}

double cnt::rA()
{
	return y3 - y2;
}
double cnt::rB()
{
	return x2 - x3;
}
double cnt::rC()
{
	return x3 * y2 - x2 * y3;
}
