#pragma once
#include "CNT.h"
class CNTInfo
{
public:
	double x1, y1, x2, y2, x3, y3, x4, y4;
	CNTInfo();
	CNTInfo(double _x1, double _y1, //����� ������ ����
			double _x2, double _y2, //������ ������ ����
			double _x3, double _y3, //������ ������� ����
			double _x4, double _y4); //����� ������� ����
	void top(double &a, double &b, double &c);
	void bottom(double &a, double &b, double &c);
	void left(double &a, double &b, double &c);
	void right(double &a, double &b, double &c);
};

