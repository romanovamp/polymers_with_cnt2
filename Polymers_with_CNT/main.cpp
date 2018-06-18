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
#include <ctime>

using namespace std;

//************************************************************************
int L;
double mean, devi, radius;
int n;
vector <CNT> cnt(0);
vector <CNT> cnt_trans(0);

ofstream file("coordinates.txt"), raspr("rasp.txt");
ofstream dd("coos.txt");
//ofstream aa("a.txt"), kk("k.txt");
bool flag;
bool *f = new bool[4];
MtRng64 mt;
bool ready = false;
double second = 0.0;
//************************************************************************
int GraphInConsole()
{
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);
	MoveToEx(hDC, 10, 100, NULL); // -
	LineTo(hDC, 10 + L, 100);

	MoveToEx(hDC, 10, 100 + L, NULL); // -
	LineTo(hDC, 10 + L, 100 + L);

	MoveToEx(hDC, 10 + L, 100, NULL); // |
	LineTo(hDC, 10 + L, 100 + L);

	MoveToEx(hDC, 10, 100, NULL); // |
	LineTo(hDC, 10, 100 + L);

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
double coord_x(double x, double k, int a)
{
	return x + k * cos(a);
}
double coord_y(double y, double k, int a)
{
	return y + k * sin(a);
}

void draw_CNT(double x, double y, double k, int a)
{
	double x1_l = coord_x(x, radius, a + 90);
	double y1_l = coord_y(y, radius, a + 90);

	double x2_l = coord_x(x1_l , k, a);
	double y2_l = coord_y(y1_l , k, a);

	double x1_r = coord_x(x, radius, a - 90);
	double y1_r = coord_y(y, radius, a - 90);

	double x2_r = coord_x(x1_r, k, a);
	double y2_r = coord_y(y1_r, k, a);

	double x2 = coord_x(x, k, a);
	double y2 = coord_y(y, k, a);

	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);

	MoveToEx(hDC, x1_l + 10, y1_l + 100, NULL);
	LineTo(hDC, x2_l + 10, y2_l + 100);

	MoveToEx(hDC, x1_r + 10, y1_r + 100, NULL);
	LineTo(hDC, x2_r + 10, y2_r + 100);


	MoveToEx(hDC, x1_l + 10, y1_l + 100, NULL);
	LineTo(hDC, x1_r + 10, y1_r + 100);

	MoveToEx(hDC, x2_l + 10, y2_l + 100, NULL);
	LineTo(hDC, x2_r + 10, y2_r + 100);


}
double d(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
bool min_distance(double x1, double y1, double x2, double y2, double x, double y) //true - допустимое расстояние
{
	double p = (d(x1, y1, x2, y2) + d(x1, y1, x, y) + d(x2, y2, x, y)) / 2.0;
	if ((2.0 / d(x1, y1, x2, y2))*sqrt(p*(p - d(x1, y1, x2, y2))*(p - d(x1, y1, x, y))*(p - d(x2, y2, x, y))) < 0.5) return false;
	return true;
}	


void equation(double x1, double y1, double x2, double y2, double &a, double &b, double &c) // уравнение прямой 
{
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}

bool belong(double x, double y)
{
	if (x >= 0 && x <= L && y >= 0 && y <= L) return true;
	else return false;
}
bool belong(double x, double y, double _x, double _y, double k, double a) //true - точка лежит внутри трубки или недопустимо близко
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
		!min_distance(x1, y1, x2, y2, x, y) ||
		!min_distance(x2, y2, x3, y3, x, y) ||
		!min_distance(x3, y3, x4, y4, x, y) ||
		!min_distance(x4, y4, x1, y1, x, y)) return true;


	return false;
}

bool check(double x1, double y1, double x2, double y2, double k1, int al1, double k2, int al2) //true - пересекаются
{
	double a1, a2, b1, b2, c1, c2;
	double x, y; //(x,y) - точка пересечения прямых

	double x3 = coord_x(x1, k1, al1);
	double y3 = coord_y(y1, k1, al1);

	double x4 = coord_x(x2, k2, al2);
	double y4 = coord_y(y2, k2, al2);

	equation(x1, y1, x3, y3, a1, b1, c1);
	equation(x2, y2, x4, y4, a2, b2, c2);

	if (!parall(a1, a2, b1, b2))
	{
		intersect(a1, a2, b1, b2, c1, c2, x, y); //(x,y) - точка пересечения прямых

		double k1 = (y2 - y1) / (x2 - x1);
		double k2 = (y4 - y1) / (x4 - x1);

		if ((k1*(x - x1) - (y - y1))*(k1*(x - x3) - (y - y3)) < 0 && (k2*(x - x1) - (y - y1))*(k2*(x - x2) - (y - y2)) < 0)
			return true;
	}
	return false;
}

bool all_test(double x, double y, double k, int a, vector<CNT>loc) //true - удачное расположение
{
	for (int i = 0; i < loc.size(); i++)
	{
		if ((sqrt(pow(k, 2) + pow(radius, 2)) + sqrt(pow(loc[i].k, 2) + pow(radius, 2))) < d(x, y, loc[i].x, loc[i].y)) continue;


		/*left -> left, bottom, right, top*/
		double x1_l = coord_x(x, radius, a + 90);
		double y1_l = coord_y(y, radius, a + 90);

		double x2_l = coord_x(loc[i].x, radius, loc[i].a + 90);
		double y2_l = coord_y(loc[i].y, radius, loc[i].a + 90);

		double x1_r = coord_x(x, radius, a - 90);
		double y1_r = coord_y(y, radius, a - 90);

		double x2_r = coord_x(loc[i].x, radius, loc[i].a - 90);
		double y2_r = coord_y(loc[i].y, radius, loc[i].a - 90);
		
		if (check(x1_l, y1_l, x2_l, y2_l, k, a, loc[i].k, loc[i].a)) return false;/*left-left*/
		if (check(x1_l, y1_l, x2_r, y2_r, k, a, loc[i].k, loc[i].a)) return false; /*left-right*/
		if (check(x1_r, y1_r, x2_l, y2_l, k, a, loc[i].k, loc[i].a)) return false; /*right-left*/
		if (check(x1_r, y1_r, x2_r, y2_r, k, a, loc[i].k, loc[i].a)) return false; /*right-right*/
		
		
		if (belong(x1_l, y1_l, loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(x1_r, y1_r, loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(coord_x(x1_l, k, a), coord_y(y1_l, k, a), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;
		if (belong(coord_x(x1_r, k, a), coord_y(y1_r, k, a), loc[i].x, loc[i].y, loc[i].k, loc[i].a)) return false;


		if (belong(x2_l, y2_l, x, y, k, a)) return false;
		if (belong(x2_r, y2_r, x, y, k, a)) return false;
		if (belong(coord_x(x2_l, loc[i].k, loc[i].a), coord_y(y2_l, loc[i].k, loc[i].a), x, y, k, a)) return false;
		if (belong(coord_x(x2_r, loc[i].k, loc[i].a), coord_y(y2_r, loc[i].k, loc[i].a), x, y, k, a)) return false;

		


		//if (!check(x2_l, y2_l, x1_l, y1_l, radius * 2, loc[i].a - 90, k, a)) return false; /*left-bottom*/
		//
		//if (!check(coord_x(x2_l, loc[i].k, loc[i].a), coord_y(y2_l, loc[i].k, loc[i].a), 
		//				x1_l, y1_l, radius * 2, loc[i].a - 90, k, a)) return false; /*left - top*/


		///*bottom -> left, bottom, right, top*/
		//if (!check(x1_l, y1_l, x2_l, y2_l, radius * 2, a - 90, loc[i].k, loc[i].a)) return false; /*bottom-left*/
		//if (!check(x1_l, y1_l, x2_l, y2_l, radius * 2, a - 90, radius * 2, loc[i].a - 90)) return false; /*bottom-bottom*/
		//if (!check(x1_l, y1_l, x2_r, y2_r, radius * 2, a - 90, loc[i].k, loc[i].a)) return false; /*bottom-right*/
		//if (!check(x1_l, y1_l, coord_x(x2_l, loc[i].k, loc[i].a), coord_y(y2_l, loc[i].k, loc[i].a),
		//	radius * 2, a - 90, radius * 2, loc[i].a - 90)) return false; /*bottom-top*/



		///*right -> left, bottom, right, top*/
		//
		//if (!check(x2_l, y2_l, x1_r, y1_r, radius * 2, loc[i].a - 90, k, a)) return false; /*right - bottom*/
		//
		//if (!check(coord_x(x2_l, loc[i].k, loc[i].a), coord_y(y2_l, loc[i].k, loc[i].a), x1_r, y1_r, radius * 2, loc[i].a - 90, k, a)) return false; /*right - top*/


		///*top -> left, bottom, right, top*/
		//if (!check(coord_x(x1_l, k, a), coord_y(y1_l, k, a), x2_l, y2_l,
		//	radius * 2, a - 90, loc[i].k, loc[i].a)) return false; /*top-left*/
		//if (!check(coord_x(x1_l, k, a), coord_y(y1_l, k, a), x2_l, y2_l,
		//	radius * 2, a - 90, radius * 2, loc[i].a - 90)) return false; /*top-bottom*/
		//if (!check(coord_x(x1_l, k, a), coord_y(y1_l, k, a), x2_r, y2_r,
		//	radius * 2, a - 90, loc[i].k, loc[i].a)) return false; /*top-right*/
		//if (!check(coord_x(x1_l, k, a), coord_y(y1_l, k, a), coord_x(x2_l, loc[i].k, loc[i].a), coord_y(y2_l, loc[i].k, loc[i].a),
		//	radius * 2, a - 90, radius * 2, loc[i].a - 90)) return false; /*top-top*/


	}
	return true;
}

bool test(double x, double y, double k, int a) //true - удачное расположение
{
	if (!all_test(x, y, k, a, cnt) || !all_test(x, y, k, a, cnt_trans)) return false;

	if (coord_x(x, k, a) < 0) if (!all_test(x + L, y, k, a, cnt) || !all_test(x + L, y, k, a, cnt_trans)) return false;
	if (coord_x(x, k, a) > L) if (!all_test(x - L, y, k, a, cnt) || !all_test(x - L, y, k, a, cnt_trans)) return false;
	if (coord_y(y, k, a) < 0) if (!all_test(x, y + L, k, a, cnt) || !all_test(x, y + L, k, a, cnt_trans)) return false;
	if (coord_y(y, k, a) > L) if (!all_test(x, y - L, k, a, cnt) || !all_test(x, y - L, k, a, cnt_trans)) return false;

	flag = true; //удачное расположение + рисуем
	return true;

}
double bm()
{
	if (ready)
	{
		ready = false;
		dd << (second * devi + mean) << endl;
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
	dd << (r * v * devi + mean) << endl;
	return (r * v * devi + mean);
}
string toStr(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}



void t(double x, double y, double k, int a)
{
	if (belong(coord_x(x, radius, a + 90), coord_y(y, radius, a + 90))) f[0] = true;
	if (belong(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90))) f[1] = true;
	if (belong(coord_x(coord_x(x, radius, a + 90), k, a), coord_y(coord_y(y, radius, a + 90), k, a))) f[2] = true;
	if (belong(coord_x(coord_x(x, radius, a - 90), k, a), coord_y(coord_y(y, radius, a - 90), k, a))) f[3] = true;
}
void trans(double x, double y, double k, int a) 
{
	if ((coord_x(x, radius, a + 90) < 0 && !f[0]) || (coord_x(x, radius, a - 90) < 0 && !f[1]) ||
		(coord_x(coord_x(x, radius, a + 90), k, a) < 0 && !f[2]) || (coord_x(coord_x(x, radius, a - 90), k, a) < 0 &&f[3]))
	{
		cnt_trans.push_back(CNT(x + L, y, a, k));
		draw_CNT(x + L, y, k, a);
		t(x + L, y, k, a);
		trans(x + L, y, k, a);
		return;
	}

	if ((coord_x(x, radius, a + 90) > L && !f[0]) || (coord_x(x, radius, a - 90) > L && !f[1]) ||
		(coord_x(coord_x(x, radius, a + 90), k, a) > L && !f[2]) || (coord_x(coord_x(x, radius, a - 90), k, a) > L && !f[3]))
	{
		cnt_trans.push_back(CNT(x - L, y, a, k));
		draw_CNT(x - L, y, k, a);
		t(x - L, y, k, a);
		trans(x - L, y, k, a);
		return;
	}
	if ((coord_y(y, radius, a + 90) < 0 && !f[0]) || (coord_y(y, radius, a - 90) < 0 && !f[1]) ||
		(coord_y(coord_y(y, radius, a + 90), k, a) < 0 && !f[2]) || (coord_y(coord_y(y, radius, a - 90), k, a) < 0 && !f[3]))
	{
		cnt_trans.push_back(CNT(x, y + L, a, k));
		draw_CNT(x, y + L, k, a);
		t(x, y + L, k, a);
		trans(x, y + L, k, a);
		return;
	}
	if ((coord_y(y, radius, a + 90) > L && !f[0]) || (coord_y(y, radius, a - 90) > L && !f[1]) ||
		(coord_y(coord_y(y, radius, a + 90), k, a) > L && !f[2]) || (coord_y(coord_y(y, radius, a - 90), k, a) > L && !f[3]))
	{
		cnt_trans.push_back(CNT(x, y - L, a, k));
		draw_CNT(x, y - L, k, a);
		t(x, y - L, k, a);
		trans(x, y - L, k, a);
		return;
	}
}
void packaging()
{
	double x, y, k;
	int a, kol = 0;
	for (int i = 0; i < n; i++)
	{
		flag = false;
		kol = 0;
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
				n--;
				break;
			}
			kol++;
		} while (!test(x, y, k, a));

		if (flag)
		{
			cnt.push_back(CNT(x, y, a, k));
			draw_CNT(x, y, k, a);
			for (int i = 0; i < 4; i++)
				f[i] = false;

			t(x, y, k, a);

			for (int i = 0; i < 4; i++)
				if(!f[i]) trans(x, y, k, a);

			double x1 = coord_x(x, radius, a + 90);
			double y1 = coord_y(y, radius, a + 90);

			double x2 = coord_x(x, radius, a - 90);
			double y2 = coord_y(y, radius, a - 90);

			double x3 = coord_x(x2, k, a);
			double y3 = coord_y(y2, k, a);

			double x4 = coord_x(x1, k, a);
			double y4 = coord_y(y1, k, a);

			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
			//file << x1 << "|" << y1 << "|" << x2 << "|" << y2 << "|" << x3 << "|" << y3 << "|" << x4 << "|" << y4 << "|" << endl;
		}
	}
}


void main()
{
	f = new bool[4];
	setlocale(LC_ALL, "rus");
	cout << "Размер квадрата: ";
	cin >> L;
	cout << "Средняя длина трубки: ";
	cin >> mean;
	cout << "Радиус трубки: ";
	cin >> radius;
	cout << "Количество трубок: ";
	cin >> n;
	

	file << "Размер квадрата: " << L << endl;
	file << "Средняя длина трубки: " << mean << endl;
	file << "Количество трубок: " << n << endl;

	GraphInConsole();
	devi = mean * 0.1 / 3.0;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	double start = clock();
	packaging();
	double finish = clock();
	cout << "meow! " << endl;


	file.close();
	raspr.close();
	delete[]f;
	cout << "Упаковано " << n << " за " << (finish-start)/CLOCKS_PER_SEC << "c." << endl;

}


