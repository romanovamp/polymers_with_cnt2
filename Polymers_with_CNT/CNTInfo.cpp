#include "CNTInfo.h"
#include "CNT.h"

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