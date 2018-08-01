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
#define M_PI 3.14159265358979323846

//************************************************************************
int GraphInConsole()
{
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);
	MoveToEx(hDC, 50, 130, NULL); // -
	LineTo(hDC, 50 + L, 130);

	MoveToEx(hDC, 50, 130 + L, NULL); // -
	LineTo(hDC, 50 + L, 130 + L);

	MoveToEx(hDC, 50 + L, 130, NULL); // |
	LineTo(hDC, 50 + L, 130 + L);

	MoveToEx(hDC, 50, 130, NULL); // |
	LineTo(hDC, 50, 130 + L);

	return 0;
}

bool parall(double a1, double a2, double b1, double b2)
{
	return (a1*b2 == a2*b1); //true - пр€мые параллельны
}
void intersect(double a1, double a2, double b1, double b2, double c1, double c2, double& x, double& y) //точка пересечени€ пр€мых
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

void draw_CNT(double x, double y, double k, double a)
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

	MoveToEx(hDC, x1_l + 50, y1_l + 130, NULL);
	LineTo(hDC, x2_l + 50, y2_l + 130);

	MoveToEx(hDC, x1_r + 50, y1_r + 130, NULL);
	LineTo(hDC, x2_r + 50, y2_r + 130);


	MoveToEx(hDC, x1_l + 50, y1_l + 130, NULL);
	LineTo(hDC, x1_r + 50, y1_r + 130);

	MoveToEx(hDC, x2_l + 50, y2_l + 130, NULL);
	LineTo(hDC, x2_r + 50, y2_r + 130);


}
double d(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
bool min_distance(double x1, double y1, double x2, double y2, double x, double y) //true - допустимое рассто€ние
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
	if (A > 90 || C > 90) return true; // можно не провер€ть
	double p = (d(x1, y1, x2, y2) + d(x1, y1, x, y) + d(x2, y2, x, y)) / 2.0;
	if ((2.0 / d(x1, y1, x2, y2))*sqrt(p*(p - d(x1, y1, x2, y2))*(p - d(x1, y1, x, y))*(p - d(x2, y2, x, y))) <= radius * 0.154) return false;
	return true;
}


void equation(double x1, double y1, double x2, double y2, double &a, double &b, double &c) // уравнение пр€мой 
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

	/*
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 0, 0));
	SelectObject(hDC, Pen);


	MoveToEx(hDC, x1 + 50+1, y1 + 130 + 1, NULL);
	LineTo(hDC, x2 + 50 + 1, y2 + 130 + 1);

	MoveToEx(hDC, x3 + 50 + 1, y3 + 130 + 1, NULL);
	LineTo(hDC, x4 + 50 + 1, y4 + 130 + 1);
	*/

	double k1 = (y2 - y1) / (x2 - x1);
	double k2 = (y4 - y1) / (x4 - x1);

	if (((k1*(x - x1) - (y - y1))*(k1*(x - x3) - (y - y3)) < 0 && (k2*(x - x1) - (y - y1))*(k2*(x - x2) - (y - y2)) < 0) ||
		!min_distance(x1, y1, x2, y2, x, y) ||
		!min_distance(x2, y2, x3, y3, x, y) ||
		!min_distance(x3, y3, x4, y4, x, y) ||
		!min_distance(x4, y4, x1, y1, x, y)) return true;


	return false;
}

bool check(double x1, double y1, double x2, double y2, double k1, double al1, double k2, double al2) //true - пересекаютс€
{
	double a1, a2, b1, b2, c1, c2;
	double x, y; //(x,y) - точка пересечени€ пр€мых
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
		intersect(a1, a2, b1, b2, c1, c2, x, y); //(x,y) - точка пересечени€ пр€мых

		double k_1 = (y2 - y1) / (x2 - x1);
		double k_2 = (coord_y(y2, k2, al2) - y1) / (coord_x(x2, k2, al2) - x1);

		if ((k_1*(x - x1) - (y - y1))*(k_1*(x - coord_x(x1, k1, al1)) - (y - coord_y(y1, k1, al1))) < 0 && 
			(k_2*(x - x1) - (y - y1))*(k_2*(x - x2) - (y - y2)) < 0)
			return true;
	}
	return false;
}

bool all_test(double x, double y, double k, double a, vector<CNT>loc) //true - удачное расположение
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

bool coincides(double x, double y, double x_real, double y_real) //true - трубка с такими координатами уже добавлена
{
	if (x == x_real && y == y_real) return true; 
	for (int i = 1; i < 4 && cnt_trans.size() >= i; i++)
		if (cnt_trans[cnt_trans.size() - i].x == x && cnt_trans[cnt_trans.size() - i].y == y) return true;
	return false;
}

void trans(double x, double y, double k, double a, double x_add, double y_add, double k_add, double a_add, int flag, double x_real, double y_real)
{
	if (flag == 3) return; //мы можем перенести к-мер максимум 3 раза
	if (!coincides(x_add, y_add - L, x_real, y_real) && check(x, y, -L/3.0, L, k, a, 5.0 * L/3.0, 0))
	{
		
		cnt_trans.push_back(CNT(x_add, y_add - L, a_add, k_add));
		draw_CNT(x_add, y_add - L, k_add, a_add);
		trans(x, y - L, k, a, x_add, y_add - L, k_add, a_add, flag+1, x_real, y_real);
	}

	if (!coincides(x_add, y_add + L, x_real, y_real) && check(x, y, -L / 3.0, 0, k, a, 5.0 * L / 3.0, 0))
	{
		cnt_trans.push_back(CNT(x_add, y_add + L, a_add, k_add));
		draw_CNT(x_add, y_add + L, k_add, a_add);
		trans(x, y + L, k, a, x_add, y_add + L, k_add, a_add, flag + 1, x_real, y_real);
	}

	if (!coincides(x_add - L, y_add, x_real, y_real) && check(x, y, L, -L / 3.0, k, a, 5.0 * L / 3.0, 90))
	{
		cnt_trans.push_back(CNT(x_add - L, y_add, a_add, k_add));
		draw_CNT(x_add - L, y_add, k_add, a_add);
		trans(x - L, y, k, a, x_add - L, y_add, k_add, a_add, flag + 1, x_real, y_real);
	}

	if (!coincides(x_add + L, y_add, x_real, y_real) && check(x, y, 0, - L / 3.0, k, a, 5.0 * L / 3.0, 90))
	{
		cnt_trans.push_back(CNT(x_add + L, y_add, a_add, k_add));
		draw_CNT(x_add + L, y_add, k_add, a_add);
		trans(x + L, y, k, a, x_add + L, y_add, k_add, a_add, flag + 1, x_real, y_real);
	}
	/*
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(2, 255, 255));
	SelectObject(hDC, Pen);


	MoveToEx(hDC, -L / 3.0 + 10, L + 100, NULL);
	LineTo(hDC, coord_x(-L / 3.0, 5.0 * L / 3.0, 0) + 10, coord_y(L, 5.0 * L / 3.0, 0) + 100);

	MoveToEx(hDC, L + 10, -L / 3.0 + 100, NULL);
	LineTo(hDC, coord_x( L, 5.0 * L / 3.0, M_PI/2.0) + 10, coord_y(-L / 3.0, 5.0 * L / 3.0, M_PI / 2.0) + 100);
	*/

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
			

			double x1 = coord_x(x, radius, a + 90);
			double y1 = coord_y(y, radius, a + 90);

			double x2 = coord_x(x, radius, a - 90);
			double y2 = coord_y(y, radius, a - 90);

			double x3 = coord_x(x2, k, a);
			double y3 = coord_y(y2, k, a);

			double x4 = coord_x(x1, k, a);
			double y4 = coord_y(y1, k, a);
			 
			if (!belong(x1, y1) || !belong(x2, y2) || !belong(x3, y3) || !belong(x4, y4))
			{
				trans(x1, y1, radius*2, a - 90, x, y, k, a,0, x, y);
				trans(x1, y1, k, a, x, y, k, a,0, x, y);
				trans(x4, y4, radius*2, a - 90, x, y, k, a,0, x, y);
				trans(x2, y2, k, a, x, y, k, a,0, x, y);
			}
			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
		}
	}
}


void main()
{
	f = new bool[4];
	setlocale(LC_ALL, "rus");
	cout << "–азмер квадрата: ";
	cin >> L;
	cout << "—редн€€ длина трубки: ";
	cin >> mean;
	cout << "–адиус трубки: ";
	cin >> radius;
	cout << " оличество трубок: ";
	cin >> n;
	

	file << "–азмер квадрата: " << L << endl;
	file << "—редн€€ длина трубки: " << mean << endl;
	file << " оличество трубок: " << n << endl;

	GraphInConsole();
	devi = mean * 0.1 / 3.0;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	double start = clock();
	packaging();
	double finish = clock();
	cout << "meow! " << endl;
	for (int i = 0; i < cnt_trans.size(); i++)
	{
		dd << cnt_trans[i].x << endl;
	}
	dd << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl;

	for (int i = 0; i < cnt_trans.size(); i++)
	{
		dd << cnt_trans[i].y << endl;
	}
	file.close();
	raspr.close();
	delete[]f;
	cout << "”паковано " << n << " за " << (finish-start)/CLOCKS_PER_SEC << "c." << endl;

}


