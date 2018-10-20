#include <fstream>
#include <iostream>
#include <cmath>
#include "mersennetwister.h"
#include <Windows.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include "Header.h"
#include "CNT.h"
#include "CNTInfo.h"
#include "Interval.h"

#include <ctime>

using namespace std;

//************************************************************************
int L;
double mean, devi, radius;
double p;
vector <CNT> cnt(0);
vector <CNT> cntTrans(0);
vector <CNTInfo> cntInfo(0);
vector <CNTInfo> cntTransInfo(0);
ofstream file, raspr, dd;
bool flag;
MtRng64 mt;
bool ready = false;
double second = 0.0;
Interval *intervals;
#define M_PI 3.14159265358979323846
int changeX = 300, changeY=50;


//************************************************************************
int GraphInConsole()
{
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);
	MoveToEx(hDC, changeX, changeY, NULL); // -
	LineTo(hDC, changeX + L, changeY);

	MoveToEx(hDC, changeX, changeY + L, NULL); // -
	LineTo(hDC, changeX + L, changeY + L);

	MoveToEx(hDC, changeX + L, changeY, NULL); // |
	LineTo(hDC, changeX + L, changeY + L);

	MoveToEx(hDC, changeX, changeY, NULL); // |
	LineTo(hDC, changeX, changeY + L);

	return 0;
}

bool parall(double a1, double a2, double b1, double b2)
{
	return (a1*b2 == a2*b1); //true - прямые параллельны
}
void intersect(double a1, double a2, double b1, double b2, double c1, double c2, double& x, double& y) //точка пересечения прямых
{
	x = (b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
	y = (a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1);
}
double coord_x(double x, double k, double a)
{
	return x + k * cos((a*M_PI)/180.0);
}
double coord_y(double y, double k, double a)
{
	return y + k * sin((a*M_PI) / 180.0);
}

void drawCNT(vector <CNTInfo> loc, int id)
{
	
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);

	MoveToEx(hDC, loc[id].x1 + changeX, loc[id].y1 + changeY, NULL);
	LineTo(hDC, loc[id].x4 + changeX, loc[id].y4 + changeY);

	MoveToEx(hDC, loc[id].x2 + changeX, loc[id].y2 + changeY, NULL);
	LineTo(hDC, loc[id].x3 + changeX, loc[id].y3 + changeY);

	MoveToEx(hDC, loc[id].x1 + changeX, loc[id].y1 + changeY, NULL);
	LineTo(hDC, loc[id].x2 + changeX, loc[id].y2 + changeY);

	MoveToEx(hDC, loc[id].x4 + changeX, loc[id].y4 + changeY, NULL);
	LineTo(hDC, loc[id].x3 + changeX, loc[id].y3 + changeY);

}
double d(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
bool minDistance(double x1, double y1, double x2, double y2, double x, double y) //true - допустимое расстояние
{
	double a = d(x2, y2, x, y);
	if (a <= radius * 0.154) return false;
	double b = d(x1, y1, x2, y2);
	double c = d(x1, y1, x, y);
	if (c <= radius * 0.154) return false;
	double A = (180.0 / M_PI) * acos(((b*b+c*c-a*a) / (2.0*b*c))*M_PI/180.0);
	if (A == 0) return false;
	double C = (180.0 / M_PI) * acos(((b*b+a*a-c*c) / (2.0*b*a))*M_PI / 180.0);
	if (C == 0) return false;
	if (A > 90 || C > 90) return true; // можно не проверять
	double p = (d(x1, y1, x2, y2) + d(x1, y1, x, y) + d(x2, y2, x, y)) / 2.0;
	if ((2.0 / d(x1, y1, x2, y2))*sqrt(p*(p - d(x1, y1, x2, y2))*(p - d(x1, y1, x, y))*(p - d(x2, y2, x, y))) <= radius * 0.154) return false;
	return true;
}


void equation(double x1, double y1, double x2, double y2, double &a, double &b, double &c) // уравнение прямой 
{ 
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}

bool belong(double x, double y) //true - точка лежит внутри основного квадрата
{
	if (x >= 0 && x <= L && y >= 0 && y <= L) return true;
	else return false;
}
bool belong(double x, double y, double _x, double _y, double k, double a) //true - точка (x, y) лежит внутри трубки или недопустимо близко
{
	if (!belong(x, y)) return false;
	double x1 = coord_x(_x, radius, a + 90);
	double y1 = coord_y(_y, radius, a + 90);

	double x2 = coord_x(_x, radius, a - 90);
	double y2 = coord_y(_y, radius, a - 90);

	double x3 = coord_x(x2, k, a);
	double y3 = coord_y(y2, k, a);

	double x4 = coord_x(x1, k, a);
	double y4 = coord_y(y1, k, a);
	
	double k1 = (y2 - y1) / (x2 - x1);
	double k2 = (y4 - y1) / (x4 - x1);

	if (((k1*(x - x1) - (y - y1))*(k1*(x - x3) - (y - y3)) < 0 && (k2*(x - x1) - (y - y1))*(k2*(x - x2) - (y - y2)) < 0) ||
		!minDistance(x1, y1, x2, y2, x, y) ||
		!minDistance(x2, y2, x3, y3, x, y) ||
		!minDistance(x3, y3, x4, y4, x, y) ||
		!minDistance(x4, y4, x1, y1, x, y)) return true;


	return false;
}

bool check(double x1, double y1, double x2, double y2, double k1, double al1, double k2, double al2) //true - пересекаются
{
	double a1, a2, b1, b2, c1, c2;
	double x, y; //(x,y) - точка пересечения прямых
	/*
	double x3 = coord_x(x1, k1, al1);
	double y3 = coord_y(y1, k1, al1);

	double x4 = coord_x(x2, k2, al2);
	double y4 = coord_y(y2, k2, al2);

	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 0, 120));
	SelectObject(hDC, Pen);


	MoveToEx(hDC, x1 + 50 + 1, y1 + 130 + 1, NULL);
	LineTo(hDC, coord_x(x1, k1, al1) + 50 + 1, coord_y(y1, k1, al1) + 130 + 1);

	MoveToEx(hDC, x2 + 50 + 1, y2 + 130 + 1, NULL);
	LineTo(hDC, coord_x(x2, k2, al2) + 50 + 1, coord_y(y2, k2, al2) + 130 + 1);
	*/

	
	equation(x1, y1, coord_x(x1, k1, al1), coord_y(y1, k1, al1), a1, b1, c1);
	equation(x2, y2, coord_x(x2, k2, al2), coord_y(y2, k2, al2), a2, b2, c2);

	if (!parall(a1, a2, b1, b2))
	{
		intersect(a1, a2, b1, b2, c1, c2, x, y); //(x,y) - точка пересечения прямых

		double k_1 = (y2 - y1) / (x2 - x1);
		double k_2 = (coord_y(y2, k2, al2) - y1) / (coord_x(x2, k2, al2) - x1);

		if ((k_1*(x - x1) - (y - y1))*(k_1*(x - coord_x(x1, k1, al1)) - (y - coord_y(y1, k1, al1))) < 0 && 
			(k_2*(x - x1) - (y - y1))*(k_2*(x - x2) - (y - y2)) < 0)
			return true;
	}
	return false;
}

bool allTest(double x, double y, double k, double a, vector<CNT>loc) //true - удачное расположение
{
	for (int i = 0; i < loc.size(); i++)
	{
		if ((sqrt(pow(k, 2) + pow(radius, 2)) + sqrt(pow(loc[i].k, 2) + pow(radius, 2))) < d(x, y, loc[i].x, loc[i].y)) continue;


		/*left -> left, bottom, right, top*/
		/*
		double x1_l = coord_x(x, radius, a + 90); 
		double y1_l = coord_y(y, radius, a + 90); 

		double x2_l = coord_x(loc[i].x, radius, loc[i].a + 90); 
		double y2_l = coord_y(loc[i].y, radius, loc[i].a + 90); 

		double x1_r = coord_x(x, radius, a - 90); 
		double y1_r = coord_y(y, radius, a - 90); 

		double x2_r = coord_x(loc[i].x, radius, loc[i].a - 90);
		double y2_r = coord_y(loc[i].y, radius, loc[i].a - 90);
		*/
		
		if (check(coord_x(x, radius, a + 90), coord_y(y, radius, a + 90), coord_x(loc[i].x, radius, loc[i].a + 90), coord_y(loc[i].y, radius, loc[i].a + 90), k, a, loc[i].k, loc[i].a)) return false;/*left-left*/
		if (check(coord_x(x, radius, a + 90), coord_y(y, radius, a + 90), coord_x(loc[i].x, radius, loc[i].a - 90), coord_y(loc[i].y, radius, loc[i].a - 90), k, a, loc[i].k, loc[i].a)) return false; /*left-right*/
		if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), coord_x(loc[i].x, radius, loc[i].a + 90), coord_y(loc[i].y, radius, loc[i].a + 90), k, a, loc[i].k, loc[i].a)) return false; /*right-left*/
		if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), coord_x(loc[i].x, radius, loc[i].a - 90), coord_y(loc[i].y, radius, loc[i].a - 90), k, a, loc[i].k, loc[i].a)) return false; /*right-right*/
		
		
		if (belong(coord_x(x, radius, a + 90), coord_y(y, radius, a + 90), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(coord_x(coord_x(x, radius, a + 90), k, a), coord_y(coord_y(y, radius, a + 90), k, a), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(coord_x(coord_x(x, radius, a - 90), k, a), coord_y(coord_y(y, radius, a - 90), k, a), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;


		if (belong(coord_x(loc[i].x, radius, loc[i].a + 90), coord_y(loc[i].y, radius, loc[i].a + 90), x, y, k, a)) return false;
		if (belong(coord_x(loc[i].x, radius, loc[i].a - 90), coord_y(loc[i].y, radius, loc[i].a - 90), x, y, k, a)) return false;
		if (belong(coord_x(coord_x(loc[i].x, radius, loc[i].a + 90), loc[i].k, loc[i].a), coord_y(coord_y(loc[i].y, radius, loc[i].a + 90), loc[i].k, loc[i].a), x, y, k, a)) return false;
		if (belong(coord_x(coord_x(loc[i].x, radius, loc[i].a - 90), loc[i].k, loc[i].a), coord_y(coord_y(loc[i].y, radius, loc[i].a - 90), loc[i].k, loc[i].a), x, y, k, a)) return false;

	}
	return true;
}
bool test(double x, double y, double k, double a) //true - удачное расположение
{
	if (!allTest(x, y, k, a, cnt) || !allTest(x, y, k, a, cntTrans)) return false;

	if (coord_x(x, k, a) < 0) if (!allTest(x + L, y, k, a, cnt) || !allTest(x + L, y, k, a, cntTrans)) return false;
	if (coord_x(x, k, a) > L) if (!allTest(x - L, y, k, a, cnt) || !allTest(x - L, y, k, a, cntTrans)) return false;
	if (coord_y(y, k, a) < 0) if (!allTest(x, y + L, k, a, cnt) || !allTest(x, y + L, k, a, cntTrans)) return false;
	if (coord_y(y, k, a) > L) if (!allTest(x, y - L, k, a, cnt) || !allTest(x, y - L, k, a, cntTrans)) return false;

	flag = true; //удачное расположение + рисуем
	return true;

}
double bm()
{
	if (ready)
	{
		ready = false;
		return (second * devi + mean);
	}
	double s = 0, u = 0, v = 0;
	do
	{
		u = 2.0 * mt.getReal1() - 1.0;
		v = 2.0 * mt.getReal1() - 1.0;
		s = u * u + v * v;
	} while (s > 1.0 || s == 0.0);

	double r = sqrt(-2.0 * log(s) / s);
	second = r * u;
	ready = true;
	return (r * v * devi + mean);
}
string toStr(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

bool coincides(double x, double y) //true - трубка с такими координатами уже добавлена
{
	for (int i = 1; i < 4 && cntTrans.size() >= i; i++)
		if (cntTrans[cntTrans.size() - i].x == x && cntTrans[cntTrans.size() - i].y == y) return true;
	return false;
}
void addInfoCNT(double x, double y, int a, double k)
{
	cntTransInfo.push_back(CNTInfo( coord_x(x, radius, a + 90), coord_y(y, radius, a + 90), /*(x1, y1)*/
							   coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), /*(x2, y2)*/
							   coord_x(coord_x(x, radius, a - 90), k, a), coord_y(coord_y(y, radius, a - 90), k, a), /*(x3, y3)*/
							   coord_x(coord_x(x, radius, a + 90), k, a), coord_y(coord_y(y, radius, a + 90), k, a)));/*(x4, y4)*/
}
void trans(double x, double y, int id) //добавление трубки в cntTrans
{
	
	if (x < 0 && !coincides(cnt[id].x + L, cnt[id].y)) //лево
	{
		cntTrans.push_back(CNT(cnt[id].x + L, cnt[id].y, cnt[id].a, cnt[id].k));
		addInfoCNT(cnt[id].x + L, cnt[id].y, cnt[id].a, cnt[id].k);
		//drawCNT(cntTransInfo, cntTransInfo.size()-1);

		if (y < 0 && !coincides(cnt[id].x + L, cnt[id].y + L)) //левый нижний угол
		{
			cntTrans.push_back(CNT(cnt[id].x + L, cnt[id].y + L, cnt[id].a, cnt[id].k));
			addInfoCNT(cnt[id].x + L, cnt[id].y + L, cnt[id].a, cnt[id].k);
			//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
		}
		if (y > L && !coincides(cnt[id].x + L, cnt[id].y - L)) //правый нижний угол
		{
			cntTrans.push_back(CNT(cnt[id].x + L, cnt[id].y - L, cnt[id].a, cnt[id].k));
			addInfoCNT(cnt[id].x + L, cnt[id].y - L, cnt[id].a, cnt[id].k);
			//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
		}
	}

	if (x > L && !coincides(cnt[id].x - L, cnt[id].y)) //право
	{
		cntTrans.push_back(CNT(cnt[id].x - L, cnt[id].y, cnt[id].a, cnt[id].k));
		addInfoCNT(cnt[id].x - L, cnt[id].y, cnt[id].a, cnt[id].k);
		//drawCNT(cntTransInfo, cntTransInfo.size() - 1);

		if (y > L && !coincides(cnt[id].x - L, cnt[id].y - L)) //правый верхний угол
		{
			cntTrans.push_back(CNT(cnt[id].x - L, cnt[id].y - L, cnt[id].a, cnt[id].k));
			addInfoCNT(cnt[id].x - L, cnt[id].y - L, cnt[id].a, cnt[id].k);
			//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
		}
		if (y < 0 && !coincides(cnt[id].x - L, cnt[id].y + L)) //левый верхний угол
		{
			cntTrans.push_back(CNT(cnt[id].x - L, cnt[id].y + L, cnt[id].a, cnt[id].k));
			addInfoCNT(cnt[id].x - L, cnt[id].y + L, cnt[id].a, cnt[id].k);
			//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
		}
	}

	if (y < 0 && !coincides(cnt[id].x, cnt[id].y + L)) //низ
	{
		cntTrans.push_back(CNT(cnt[id].x, cnt[id].y + L, cnt[id].a, cnt[id].k));
		addInfoCNT(cnt[id].x, cnt[id].y + L, cnt[id].a, cnt[id].k);
		//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
	}
	if (y > L && !coincides(cnt[id].x, cnt[id].y - L)) //верх
	{
		cntTrans.push_back(CNT(cnt[id].x, cnt[id].y - L, cnt[id].a, cnt[id].k));
		addInfoCNT(cnt[id].x, cnt[id].y - L, cnt[id].a, cnt[id].k);
		//drawCNT(cntTransInfo, cntTransInfo.size() - 1);
	}
}
void packaging()
{
	double x, y, k;
	int a, kol = 0;
	double S = 0;
	while (S < L * L * p)
	{
		flag = false;
		kol = 0;
		k = bm();
		while (k < mean - devi || k > mean + devi)
			k = bm();
		a = (int)(mt.getReal1() * 360);
		do
		{
			if (kol <= L*L)
			{
				x = mt.getReal1()*L;
				y = mt.getReal1()*L;
			}
			else
			{
				
				cout << "Реальная плотность: "<< S / (L*L) << endl;
				S = L * L * p;
				break;
			}
			kol++;
		} while (!test(x, y, k, a));

		if (flag)
		{
			cnt.push_back(CNT(x, y, a, k));
			int id = cnt.size() - 1;
			cntInfo.push_back(CNTInfo(  coord_x(x, radius, a + 90), coord_y(y, radius, a + 90), /*(x1, y1)*/
										coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), /*(x2, y2)*/
										coord_x(coord_x(x, radius, a - 90), k, a), coord_y(coord_y(y, radius, a - 90), k, a), /*(x3, y3)*/
										coord_x(coord_x(x, radius, a + 90), k, a), coord_y(coord_y(y, radius, a + 90), k, a)));/*(x4, y4)*/
			drawCNT(cntInfo, id);
			if (!belong(cntInfo[id].x1, cntInfo[id].y1)) trans(cntInfo[id].x1, cntInfo[id].y1, id);
			if (!belong(cntInfo[id].x2, cntInfo[id].y2)) trans(cntInfo[id].x2, cntInfo[id].y2, id);
			if (!belong(cntInfo[id].x3, cntInfo[id].y3)) trans(cntInfo[id].x3, cntInfo[id].y3, id);
			if (!belong(cntInfo[id].x4, cntInfo[id].y4)) trans(cntInfo[id].x4, cntInfo[id].y4, id);

			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
			S += k * radius * 2.0;
		}
	}
	cout << "Реальная плотность: " << S / (L*L) << endl;
}
int numIntervals()
{
	return pow(2.0 * cnt.size() / 3.0, 1.0 / 3.0);
}

void sort(double** A, int n, int m)
{
	double* t = NULL;
	for (int i=0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			if (A[i][0] >= A[j][0])
			{
				t = A[i];
				A[i] = A[j];
				A[j] = t;
			}
}

void cut(int j, double &x1, double &y1, double &x2, double &y2, double A1, double B1, double C1)
{
	double A2, B2, C2;
	bool dot1 = false;
	double x, y;

	cntInfo[j].left(A2, B2, C2);
	intersect(A1, A2, B1, B2, C1, C2, x, y);
	if (cntInfo[j].y1 <= y && y <= cntInfo[j].y4 || cntInfo[j].y4 <= y && y <= cntInfo[j].y1) //точка подходит
	{
		x1 = x;
		y1 = y;
		dot1 = true; //первая точка есть
	}

	cntInfo[j].bottom(A2, B2, C2);
	intersect(A1, A2, B1, B2, C1, C2, x, y);
	if (cntInfo[j].y1 <= y && y <= cntInfo[j].y2 || cntInfo[j].y2 <= y && y <= cntInfo[j].y1) //точка подходит
		if (!dot1)
		{
			x1 = x;
			y1 = y;
			dot1 = true; //первая точка есть
		}
		else 
		{
			x2 = x;
			y2 = y;
			return;
		}

	cntInfo[j].right(A2, B2, C2);
	intersect(A1, A2, B1, B2, C1, C2, x, y);
	if (cntInfo[j].y2 <= y && y <= cntInfo[j].y3 || cntInfo[j].y3 <= y && y <= cntInfo[j].y2) //точка подходит
		if (!dot1)
		{
			x1 = x;
			y1 = y;
			dot1 = true; //первая точка есть
		}
		else
		{
			x2 = x;
			y2 = y;
			return;
		}

	cntInfo[j].top(A2, B2, C2);
	intersect(A1, A2, B1, B2, C1, C2, x, y);
	if (cntInfo[j].y3 <= y && y <= cntInfo[j].y4 || cntInfo[j].y4 <= y && y <= cntInfo[j].y3) //точка подходит
		if (!dot1)
		{
			x1 = x;
			y1 = y;
			dot1 = true; //первая точка есть
		}
		else
		{
			x2 = x;
			y2 = y;
			return;
		}
}

void equability()
{
	int num = 2 + numIntervals();
	intervals = new Interval[num];
	for (int i = 0; i < num; i++)
	{
		intervals[i].x1 = (i-1) *L / numIntervals();
		intervals[i].x2 = i*L / numIntervals();
		intervals[i].summ = 0;

		HDC hDC = GetDC(GetConsoleWindow());
		HPEN Pen = CreatePen(PS_SOLID, 2, RGB(0, 255, 255));
		SelectObject(hDC, Pen);

		MoveToEx(hDC, changeX +intervals[i].x1, changeY, NULL); // |
		LineTo(hDC, changeX + intervals[i].x1, changeY + L);

		Pen = CreatePen(PS_SOLID, 2, RGB(0, 0, 210));
		SelectObject(hDC, Pen);

		MoveToEx(hDC, changeX + intervals[i].x2, changeY, NULL); // |
		LineTo(hDC, changeX + intervals[i].x2, changeY + L);
	}
	int pn = 4, pm = 3;
	double **point = new double*[pn]; //в каком интервале лежит точка
	for (int i = 0; i < pn; i++)
		point[i] = new double[pm];

	for (int j = 0; j < cnt.size(); j++)
	{

		for (int i = 0; i < num; i++)
		{
			if ( intervals[i].x1 <= cntInfo[j].x1  && cntInfo[j].x1 <= intervals[i].x2) { point[0][0] = i; point[0][1] = cntInfo[j].x1; point[0][2] = cntInfo[j].y1; }
			if ( intervals[i].x1 <= cntInfo[j].x2  && cntInfo[j].x2 <= intervals[i].x2) { point[1][0] = i; point[1][1] = cntInfo[j].x2; point[1][2] = cntInfo[j].y2; }
			if ( intervals[i].x1 <= cntInfo[j].x3  && cntInfo[j].x3 <= intervals[i].x2) { point[2][0] = i; point[2][1] = cntInfo[j].x3; point[2][2] = cntInfo[j].y3; }
			if ( intervals[i].x1 <= cntInfo[j].x4  && cntInfo[j].x4 <= intervals[i].x2) { point[3][0] = i; point[3][1] = cntInfo[j].x4; point[3][2] = cntInfo[j].y4; }
		}
		//1. унт внутри интервала
		if (point[0][0] == point[1][0] && point[1][0] == point[2][0] && point[2][0] == point[3][0])
		{
			intervals[(int)point[0][0]].summ = intervals[(int)point[0][0]].summ + cnt[j].k*radius*2.0;
			continue;
		}
		//2. если в разных интервалах 

		//сотировка по первой строке 
		sort(point, pn, pm);
		

		double A = 1, B = 0, C = -intervals[(int)point[0][0]].x2;
		double x1, y1, x2, y2;
		cut(j, x1, y1, x2, y2, A, B, C); // координаты точек пересечения с границей интервала
		double s = 0;

		//2.1 случай 1222 треугольник-многоугольник
		if (point[0][0] != point[1][0] && point[1][0] == point[2][0] && point[2][0] == point[3][0])
		{
			double h = abs(A*point[0][1] + B*point[0][2] + C) / sqrt(A*A+B*B); //высота треугольника
			s = d(x1, y1, x2, y2)*h/2.0; //площадь треугольника
			//cout << "ТЛ   " << s << "     " << cnt[j].k * 2.0 * radius - s << "     " << cnt[j].k * 2.0 * radius << endl;
		}
		//2.2 случай 1122 трапеция-трапеция
		if (point[0][0] == point[1][0] && point[1][0] != point[2][0] && point[2][0] == point[3][0])
		{
			s = d((point[0][1] + point[1][1]) / 2.0, (point[0][2] + point[1][2]) / 2.0, (x1 + x2) / 2.0, (y1 + y2) / 2.0) * radius * 2.0; //площадь треугольника
			//cout << "ПП   " << s << "     " << cnt[j].k * 2.0 * radius - s << "     " << cnt[j].k * 2.0 * radius << endl;
		}
		//2.3 случай 1112 трапеция-многоугольник
		if (point[0][0] == point[1][0] && point[1][0] == point[2][0] && point[2][0] != point[3][0])
		{
			double h = abs(A*point[3][1] + B * point[3][2] + C) / sqrt(A * A + B * B); //высота треугольника справа
			s = (cnt[j].k * 2 * radius) - d(x1, y1, x2, y2)*h / 2.0; //площадь треугольника	
			//cout << "ТП   " << s << "     " << cnt[j].k * 2.0 * radius - s << "     " << cnt[j].k * 2.0 * radius << endl;
		}
		intervals[(int)point[0][0]].summ = intervals[(int)point[0][0]].summ + s; //площадь треугольника в левый интервал
		intervals[((int)point[0][0] + 1)].summ = intervals[((int)point[0][0] + 1)].summ + (cnt[j].k * 2 * radius - s);//площадь оставшейся фигуры в правый инервал

		//3 унт лежит более чем в двух интервалах
		if (point[0][0] != point[3][0])
		{

		}

	}
	delete[]point;
}

void main()
{
	CreateDirectoryW(L"files", NULL);
	file.open("./files/coordinates.txt");
	raspr.open("./files/rasp.txt");
	dd.open("./files/coos.txt");

	//ofstream aa("a.txt"), kk("k.txt");
	setlocale(LC_ALL, "rus");
	cout << "Размер квадрата: ";
	cin >> L;
	cout << "Средняя длина трубки: ";
	cin >> mean;
	cout << "Радиус трубки: ";
	cin >> radius;
	cout << "Плотность: ";
	cin >> p;
	
	file << "Размер квадрата: " << L << endl;
	file << "Средняя длина трубки: " << mean << endl;
	file << "Радиус: " << radius << endl;
	file << "Плотность: " << p << endl;

	GraphInConsole();
	devi = mean * 0.1;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	double start = clock();
	packaging();

	double finish = clock();
	cout << "meow! " << endl;
	/*for (int i = 0; i < cntTrans.size(); i++)
	{
		dd << cntTrans[i].x << "  " << cntTrans[i].y << endl;
	}*/
	
	file.close();
	raspr.close();
	dd.close();
	cout << "Упаковано " << cnt.size() << " за " << (finish-start)/CLOCKS_PER_SEC << "c." << endl;
	equability();
	double summ1 = 0;
	for (int j = 0; j < 2 + numIntervals(); j++)
	{
		summ1 = summ1 + intervals[j].summ;
		cout << intervals[j].x1 << "//" << intervals[j].x2 << "//" << intervals[j].summ << endl;
	}

	double summ2 = 0;
	for (int i = 0; i < cnt.size(); i++)
		summ2 = summ2 + cnt[i].k*2.0*radius;
	cout << summ1 << " = " << summ2 << endl;
	delete[]intervals;
	cin >> p; 


}


