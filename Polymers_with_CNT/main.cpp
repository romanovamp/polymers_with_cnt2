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
int L, N, n;
double mean, devi, radius;
double p;

vector <CNT> cnt(0);
vector <CNT> cntTrans(0);
vector <CNTInfo> cntInfo(0);
vector <CNTInfo> cntTransInfo(0);
ofstream file, raspr, dd, aa;
bool flag;
bool *transFlag;
MtRng64 mt;
bool ready = false, firstTest;
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
double coordX(double x, double k, int a)
{
	return x + k * cos((a*M_PI)/180.0);
}
double coordY(double y, double k, int a)
{
	return y + k * sin((a*M_PI) / 180.0);
}

void drawCNT(CNTInfo loc, bool red)
{
	
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen;
	if(red)
		Pen = CreatePen(PS_SOLID, 2, RGB(255, 0, 0));
	else 
		Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);

	MoveToEx(hDC, loc.x1 + changeX, loc.y1 + changeY, NULL);
	LineTo(hDC, loc.x4 + changeX, loc.y4 + changeY);

	MoveToEx(hDC, loc.x2 + changeX, loc.y2 + changeY, NULL);
	LineTo(hDC, loc.x3 + changeX, loc.y3 + changeY);

	MoveToEx(hDC, loc.x1 + changeX, loc.y1 + changeY, NULL);
	LineTo(hDC, loc.x2 + changeX, loc.y2 + changeY);

	MoveToEx(hDC, loc.x4 + changeX, loc.y4 + changeY, NULL);
	LineTo(hDC, loc.x3 + changeX, loc.y3 + changeY);

}
double d(double x1, double y1, double x2, double y2) //расстояние от точки до точки
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
double dTP(double x, double y, double A, double B, double C) //расстояние от точки до прямой
{
	return abs(A * x + B * y + C) / sqrt(A * A + B * B);
}

double D(double x1, double y1, double x2, double y2, double x, double y)
{
	return (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1);
}
bool minDistance(CNTInfo cI, double x, double y) //true - допустимое расстояние
{
	//double c_c = 0.154;
	double c_c = 2;

	if ((cI.lA()*x + cI.lB()*y + cI.lC())*(cI.rA()*x + cI.rB()*y + cI.rC()) >= 0)
	{
		if (dTP(x, y, cI.bA(), cI.bB(), cI.bC()) < radius * c_c) return false;
		if (dTP(x, y, cI.tA(), cI.tB(), cI.tC()) < radius * c_c) return false;
	}
	if ((cI.bA()*x + cI.bB()*y + cI.bC())*(cI.tA()*x + cI.tB()*y + cI.tC()) >= 0) 
	{
		if (dTP(x, y, cI.rA(), cI.rB(), cI.rC()) < radius * c_c) return false;
		if (dTP(x, y, cI.lA(), cI.lB(), cI.lC()) < radius * c_c) return false;

		
	}
	if (d(x, y, cI.x1, cI.y1) < radius * c_c) return false;
	if (d(x, y, cI.x2, cI.y2) < radius * c_c) return false;
	if (d(x, y, cI.x3, cI.y3) < radius * c_c) return false;
	if (d(x, y, cI.x4, cI.y4) < radius * c_c) return false;

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
bool belong(double x, double y, CNTInfo cI) //true - точка (x, y) лежит внутри трубки или недопустимо близко
{
	if ((cI.lA()*x + cI.lB()*y + cI.lC())*(cI.rA()*x + cI.rB()*y + cI.rC()) > 0 &&
		(cI.bA()*x + cI.bB()*y + cI.bC())*(cI.tA()*x + cI.tB()*y + cI.tC()) > 0) return true;
	return false;
}

bool check(double x1, double y1, double x2, double y2, double k1, double al1, double k2, double al2) //true - пересекаются
{
	double a1, a2, b1, b2, c1, c2;
	double x, y; //(x,y) - точка пересечения прямых

	equation(x1, y1, coordX(x1, k1, al1), coordY(y1, k1, al1), a1, b1, c1);
	equation(x2, y2, coordX(x2, k2, al2), coordY(y2, k2, al2), a2, b2, c2);

	if (!parall(a1, a2, b1, b2))
	{
		intersect(a1, a2, b1, b2, c1, c2, x, y); //(x,y) - точка пересечения прямых

		double k_1 = (y2 - y1) / (x2 - x1);
		double k_2 = (coordY(y2, k2, al2) - y1) / (coordX(x2, k2, al2) - x1);

		if ((k_1*(x - x1) - (y - y1))*(k_1*(x - coordX(x1, k1, al1)) - (y - coordY(y1, k1, al1))) <= 0 && 
			(k_2*(x - x1) - (y - y1))*(k_2*(x - x2) - (y - y2)) <= 0)
			return true;
	}
	return false;
}

bool allTest(CNTInfo cntNew, vector<CNT>loc, vector <CNTInfo> locInfo) //true - удачное расположение
{
	for (int i = 0; i < loc.size(); i++)
	{
		if ((sqrt(pow(cntNew.k, 2) + pow(radius, 2)) + sqrt(pow(loc[i].k, 2) + pow(radius, 2))) < d(cntNew.x, cntNew.y, loc[i].x, loc[i].y)) continue;

		double c_c = 0.154*radius;
		
		//if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), locInfo[i].x1, locInfo[i].y1, k, a, loc[i].k, loc[i].a)) return false; /*right-left*/
		//if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), locInfo[i].x2, locInfo[i].y2, k, a, loc[i].k, loc[i].a)) return false; /*right-right*/
		
		/*от точки до прямой*/
		
		double x=0, y=0;

		x = coordX(locInfo[i].x, c_c, locInfo[i].a - 180);
		y = coordY(locInfo[i].y, c_c, locInfo[i].a - 180);

		CNTInfo cI = CNTInfo(x,y, locInfo[i].k + 2 * c_c, locInfo[i].a, radius + c_c);

		if (belong(cntNew.x1, cntNew.y1, cI)) return false;
		if (belong(cntNew.x2, cntNew.y2, cI)) return false;
		if (belong(cntNew.x3, cntNew.y3, cI)) return false;
		if (belong(cntNew.x4, cntNew.y4, cI)) return false;
		
		x = coordX(cntNew.x, c_c, cntNew.a - 180);
		y = coordY(cntNew.y, c_c, cntNew.a - 180);

		cI = CNTInfo(x, y, cntNew.k + 2 * c_c, cntNew.a, radius + c_c);

		if (belong(locInfo[i].x1, locInfo[i].y1, cI)) return false;
		if (belong(locInfo[i].x2, locInfo[i].y2, cI)) return false;
		if (belong(locInfo[i].x3, locInfo[i].y3, cI)) return false;
		if (belong(locInfo[i].x4, locInfo[i].y4, cI)) return false;
		
		//if (check(cntNew.x1, cntNew.y1, locInfo[i].x1, locInfo[i].y1, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false;/*left-left*/
		if (check(cntNew.x1, cntNew.y1, locInfo[i].x2, locInfo[i].y2, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false; /*left-right*/

	}
	return true;
}
bool test(double x, double y, double k, int a) //true - удачное расположение
{
	CNTInfo cIM = CNTInfo(x, y, k, a, radius);
	CNTInfo cI;
	if (!allTest(cIM, cnt, cntInfo) || !allTest(cIM, cntTrans, cntTransInfo)) return false;

	if (cIM.x1 < 0 || cIM.x2 < 0 || cIM.x3 < 0 || cIM.x4 < 0)
	{
		cI = CNTInfo(x + L, y, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[0] = true;
	}
	if (cIM.x1 > L || cIM.x2 > L || cIM.x3 > L || cIM.x4 > L)
	{
		cI = CNTInfo(x - L, y, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[1] = true;
	}
	if (cIM.y1 < 0 || cIM.y2 < 0 || cIM.y3 < 0 || cIM.y4 < 0)
	{
		cI = CNTInfo(x, y + L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[2] = true;
	}
	if (cIM.y1 > L || cIM.y2 > L || cIM.y3 > L || cIM.y4 > L)
	{
		cI = CNTInfo(x, y - L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[3] = true;
	}



	if (cIM.x1 < 0 && cIM.y1 < 0 || cIM.x2 < 0 && cIM.y2 < 0 || cIM.x3 < 0 && cIM.y3 < 0 || cIM.x4 < 0 && cIM.y4 < 0)
	{
		cI = CNTInfo(x + L, y + L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[4] = true;
	}
	if (cIM.x1 > L && cIM.y1 < 0 || cIM.x2 > L && cIM.y2 < 0 || cIM.x3 > L && cIM.y3 < 0 || cIM.x4 > L && cIM.y4 < 0)
	{
		cI = CNTInfo(x - L, y + L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[5] = true;
	}

	if (cIM.x1 > L && cIM.y1 > L || cIM.x2 > L && cIM.y2 > L || cIM.x3 > L && cIM.y3 > L || cIM.x4 > L && cIM.y4 > L)
	{
		cI = CNTInfo(x + L, y + L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[6] = true;
	}
	if (cIM.x1 < 0 && cIM.y1 > L || cIM.x2 < 0 && cIM.y2 > L || cIM.x3 < 0 && cIM.y3 > L || cIM.x4 < 0 && cIM.y4 > L)
	{
		cI = CNTInfo(x + L, y - L, k, a, radius);
		if (!allTest(cI, cnt, cntInfo) || !allTest(cI, cntTrans, cntTransInfo)) return false;
		transFlag[7] = true;
	}



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

void addTransCNT(double x, double y, int a, double k)
{
	cntTrans.push_back(CNT(x, y, a, k));
	cntTransInfo.push_back(CNTInfo(x, y, k, a, radius));
	if (firstTest)drawCNT(cntTransInfo[(cntTransInfo.size() - 1)], true);
}
void packaging()
{
	double x, y, k;
	int a, kol = 0;
	double S = 0;

	for(int i=0; i<n; i++)
	{
		flag = false;
		kol = 0;
		k = bm();
		while (k < mean - devi || k > mean + devi)
			k = bm();
		a = (int)(mt.getReal1() * 360);
		do
		{
			for (int j = 0; j<8; j++)
				transFlag[j] = false;
			if (kol <= L*L)
			{
				
				x = mt.getReal1()*L;
				y = mt.getReal1()*L;
			}
			else
			{
				//cout << "Реальная плотность: "<< S / (L*L) << endl;
				//?????
				n = n - 1;
				break;
			}
			kol++;
		} while (!test(x, y, k, a));

		if (flag)
		{
			cnt.push_back(CNT(x, y, a, k));
			cntInfo.push_back(CNTInfo(x, y, k, a, radius));

			if(firstTest)drawCNT(cntInfo[(cntInfo.size()-1)], false);

			if (transFlag[0]) addTransCNT(x + L, y, a, k);
			if (transFlag[1]) addTransCNT(x - L, y, a, k);
			if (transFlag[2]) addTransCNT(x, y + L, a, k);
			if (transFlag[3]) addTransCNT(x, y - L, a, k);
			if (transFlag[4]) addTransCNT(x + L , y + L, a, k);
			if (transFlag[5]) addTransCNT(x - L, y + L, a, k);
			if (transFlag[6]) addTransCNT(x - L, y - L, a, k);
			if (transFlag[7]) addTransCNT(x + L, y - L, a, k);


			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
			S += k * radius * 2.0;
		}
	}
	file << "Реальная плотность: " << S / (L*L) << endl;

	cout << "реальная p: " << S / (L*L) << " ";
}
int numIntervals()
{
	return pow(2.0 * n / 3.0, 1.0 / 3.0);
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
		if (firstTest)
		{
			HDC hDC = GetDC(GetConsoleWindow());
			HPEN Pen = CreatePen(PS_SOLID, 2, RGB(0, 0, 210));
			SelectObject(hDC, Pen);

			MoveToEx(hDC, changeX + intervals[i].x1, changeY, NULL); // |
			LineTo(hDC, changeX + intervals[i].x1, changeY + L);

			Pen = CreatePen(PS_SOLID, 2, RGB(0, 0, 210));
			SelectObject(hDC, Pen);

			MoveToEx(hDC, changeX + intervals[i].x2, changeY, NULL); // |
			LineTo(hDC, changeX + intervals[i].x2, changeY + L);
		}
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
			double h =dTP(point[0][1], point[0][2], A, B, C); //высота треугольника (расстояние от точки до прямой)
			s = d(x1, y1, x2, y2)*h/2.0; //площадь треугольника
		}
		//2.2 случай 1122 трапеция-трапеция
		if (point[0][0] == point[1][0] && point[1][0] != point[2][0] && point[2][0] == point[3][0])
			s = d((point[0][1] + point[1][1]) / 2.0, (point[0][2] + point[1][2]) / 2.0, (x1 + x2) / 2.0, (y1 + y2) / 2.0) * radius * 2.0; //площадь треугольника

		//2.3 случай 1112 трапеция-многоугольник
		if (point[0][0] == point[1][0] && point[1][0] == point[2][0] && point[2][0] != point[3][0])
		{
			double h = abs(A*point[3][1] + B * point[3][2] + C) / sqrt(A * A + B * B); //высота треугольника справа
			s = (cnt[j].k * 2 * radius) - d(x1, y1, x2, y2)*h / 2.0; //площадь треугольника	
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

void intervalsAlpha(int *summAlpha, int numInterAlpha)
{
	int*a = new int[numInterAlpha];
	for (int i = 0; i < numInterAlpha; i++)
		a[i] = 0;
	for (int j = 0; j < n; j++)
	{

		if (cnt[j].a >= 0 && cnt[j].a <= 40) a[0] += 1;
		if (cnt[j].a >= 41 && cnt[j].a <= 80) a[1] += 1;
		if (cnt[j].a >= 81 && cnt[j].a <= 120) a[2] += 1;
		if (cnt[j].a >= 121 && cnt[j].a <= 160) a[3] += 1;
		if (cnt[j].a >= 161 && cnt[j].a <= 200) a[4] += 1;
		if (cnt[j].a >= 201 && cnt[j].a <= 240) a[5] += 1;
		if (cnt[j].a >= 241 && cnt[j].a <= 280) a[6] += 1;
		if (cnt[j].a >= 281 && cnt[j].a <= 320) a[7] += 1;
		if (cnt[j].a >= 321 && cnt[j].a <= 360) a[8] += 1;
	}
	for (int j = 0; j < numInterAlpha; j++)
		summAlpha[j] = summAlpha[j] + a[j];
	delete[]a;
}
void main()
{
	CreateDirectoryW(L"files", NULL);
	file.open("./files/coordinates.txt");
	raspr.open("./files/rasp_s.txt");
	dd.open("./files/raspr_a.txt");
	aa.open("./files/test.txt");


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
	cout << "Количество испытаний: ";
	cin >> N;

	file << "Размер квадрата: " << L << endl;
	file << "Средняя длина трубки: " << mean << endl;
	file << "Радиус: " << radius << endl;
	file << "Плотность: " << p << endl;
	file << "Количество испытаний: " << N << endl;

	n = p * L * L / (mean * 2 * radius);
	cout << "Всего: " << n << endl;
	file << "Всего: " << n << endl;
	GraphInConsole();
	devi = mean * 0.1;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	transFlag = new bool[8];
	int numInterRaspr = numIntervals(), numInterAlpha = 9;
	double* summInterRaspr = new double[numInterRaspr];
	int* summAlpha = new int[numInterAlpha];
	for (int j = 0; j < numInterAlpha; j++)
		summAlpha[j] = 0;
	for (int j = 0; j < numInterRaspr; j++)
		summInterRaspr[j] = 0;
	for (int i = 0; i < N; i++)
	{
		n = p * L * L / (mean * 2 * radius);
		if (i == 0) firstTest = true;
		else firstTest = false;
		file << "********************************************************************" << endl;
		file << "Испытание" + (i+1) << endl;
		file << "********************************************************************" << endl;
		cout << (i + 1) << "Исп: ";
		double start = clock();
		packaging();
		double finish = clock();
		cout << "уп. за " << (finish - start) / CLOCKS_PER_SEC << "c." << endl;
		equability();
		
		if (firstTest)
			cout << "Примерное время: " << N * (((finish - start) / CLOCKS_PER_SEC) / 60) / 60 << "ч" << endl;
		

		summInterRaspr[0] = summInterRaspr[0] + intervals[(numInterRaspr + 1)].summ;
		summInterRaspr[(numInterRaspr - 1)] = summInterRaspr[(numInterRaspr - 1)] + intervals[ 0].summ;
		
		for (int j = 0; j < numInterRaspr; j++)
			summInterRaspr[j] = summInterRaspr[j] + intervals[j+1].summ;
		
		intervalsAlpha(summAlpha, numInterAlpha);
		
		delete[]intervals;
		cnt.clear();
		cntTrans.clear();
		cntInfo.clear();
		cntTransInfo.clear();
	}

	for (int j = 0; j < numInterRaspr; j++)
		raspr << summInterRaspr[j]/N << endl;

	for (int j = 0; j < numInterAlpha; j++)
		dd << (double)summAlpha[j] / N << endl;

	delete[]summAlpha;
	delete[]summInterRaspr;
	file.close();
	raspr.close();
	dd.close();
	cout << "meow! " << endl;
	cin >> p; 


}


