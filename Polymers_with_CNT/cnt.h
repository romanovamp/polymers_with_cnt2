#pragma once
#include <vector>
using std::vector;

class cnt
{
public:
	double x1, y1, x2, y2, x3, y3, x4, y4, x, y, k;
	int a, idClus, idParent;
	int parent; 
	vector <int> m;
	//parent = 0 - не дочерний, не родитель
	//parent = 1 - родитель
	//parent = 2 - дочерний
	cnt();
	cnt(double x, double y, double k, int a, double r, int idParent, int parent);
	cnt(double x, double y, double k, int a, double r);
	void top(double &a, double &b, double &c);
	void bottom(double &a, double &b, double &c);
	void left(double &a, double &b, double &c);
	void right(double &a, double &b, double &c);
	double tA();
	double tB();
	double tC();
	double bA();
	double bB();
	double bC();
	double lA();
	double lB();
	double lC();
	double rA();
	double rB();
	double rC();
};

