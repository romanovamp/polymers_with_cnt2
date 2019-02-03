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

int colour1=0, colour2=0, colour3=0;
vector <CNT> cnt(0);
vector <CNT> cntTrans(0);
vector <CNTInfo> cntInfo(0);
vector <CNTInfo> cntTransInfo(0);
ofstream file, raspr, dd, aa;
bool flag, fla=false;
bool *transFlag;
MtRng64 mt;
bool ready = false, firstTest;
double second = 0.0;
#define M_PI 3.14159265358979323846
int changeX = 350, changeY=50;
vector <vector<CNTInfo>> clustersInfo(0);
int kolClus = 0;
double mCh = 0;
vector <CNTInfo> allCnt(0);
CNTInfo cI;
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

void drawCNT(CNTInfo loc, int i, int j,int k)
{
	
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen;
	Pen = CreatePen(PS_SOLID, 2, RGB(i, j, k));
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
bool vzaim(CNTInfo cntNew, CNTInfo locInfo, double sv) //true - не пересекаются
{
	
	cI = CNTInfo(coordX(locInfo.x, sv, locInfo.a - 180), coordY(locInfo.y, sv, locInfo.a - 180), locInfo.k + 2 * sv, locInfo.a, radius + sv);

	//if(fla)drawCNT(cI,220,220,220);
	if (belong(cntNew.x1, cntNew.y1, cI)) return false;
	if (belong(cntNew.x2, cntNew.y2, cI)) return false;
	if (belong(cntNew.x3, cntNew.y3, cI)) return false;
	if (belong(cntNew.x4, cntNew.y4, cI)) return false;
	return true;
}

bool vzaim(CNTInfo cntNew, CNTInfo locInfo, double sv1, double sv2) //true - не пересекаются
{

	CNTInfo cI1 = CNTInfo(coordX(locInfo.x, sv1, locInfo.a - 180), coordY(locInfo.y, sv1, locInfo.a - 180), locInfo.k + 2 * sv1, locInfo.a, radius + sv1);
	CNTInfo cI2 = CNTInfo(coordX(cntNew.x, sv2, cntNew.a - 180), coordY(cntNew.y, sv2, cntNew.a - 180), cntNew.k + 2 * sv2, cntNew.a, radius + sv2);

	//if(fla)drawCNT(cI,220,220,220);
	if (belong(cI1.x1, cI1.y1, cI2)) return false;
	if (belong(cI1.x2, cI1.y2, cI2)) return false;
	if (belong(cI1.x3, cI1.y3, cI2)) return false;
	if (belong(cI1.x4, cI1.y4, cI2)) return false;

	return true;
}

bool vzaim(CNTInfo cntNew, CNTInfo locInfo) //true - подходит
{
	//if(fla)drawCNT(cI,220,220,220); 
	if (belong(cntNew.x1, cntNew.y1, locInfo)) return false;
	if (belong(cntNew.x2, cntNew.y2, locInfo)) return false;
	if (belong(cntNew.x3, cntNew.y3, locInfo)) return false;
	if (belong(cntNew.x4, cntNew.y4, locInfo)) return false;
	return true;
}
bool allTest(CNTInfo cntNew, vector<CNT>loc, vector <CNTInfo> locInfo) //true - удачное расположение
{
	for (int i = 0; i < loc.size(); i++)
	{
		if ((sqrt(pow(cntNew.k, 2) + pow(radius, 2)) + sqrt(pow(loc[i].k, 2) + pow(radius, 2))) < d(cntNew.x, cntNew.y, loc[i].x, loc[i].y)) continue;

		//double c_c = 0.154*radius;
		
		//if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), locInfo[i].x1, locInfo[i].y1, k, a, loc[i].k, loc[i].a)) return false; /*right-left*/
		//if (check(coord_x(x, radius, a - 90), coord_y(y, radius, a - 90), locInfo[i].x2, locInfo[i].y2, k, a, loc[i].k, loc[i].a)) return false; /*right-right*/
		
		/*от точки до прямой*/
		
		if (!vzaim(cntNew, locInfo[i])) return false;

		if (!vzaim(locInfo[i], cntNew)) return false;
		
		if (check(cntNew.x1, cntNew.y1, locInfo[i].x1, locInfo[i].y1, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false;/*left-left*/
		if (check(cntNew.x1, cntNew.y1, locInfo[i].x2, locInfo[i].y2, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false; /*left-right*/

	}
	return true;
}
bool test(double x, double y, double k, int a) //true - удачное расположение
{
	CNTInfo cIM = CNTInfo(x, y, k, a, radius);

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
	if (firstTest)drawCNT(cntTransInfo[(cntTransInfo.size() - 1)], 225,0,0);
}
void packaging()
{
	double x, y, k;
	int a, kol = 0;
	double S = 0;
	int fdf = 0;
	for(int i=0; i<n; i++)
	{
		//if (i % 1000 == 0) cout << i << endl;
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
				cout << " -n " << endl;
				//?????
				n = n - 1;
				return;
			}
			kol++;
		} while (!test(x, y, k, a));

		if (flag)
		{
			cnt.push_back(CNT(x, y, a, k));
			cntInfo.push_back(CNTInfo(x, y, k, a, radius));

			if(firstTest)drawCNT(cntInfo[(cntInfo.size()-1)], 225,225,225);

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

void sortA(double** A, int n, int m)
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

vector<CNTInfo> sortV(vector<CNTInfo> A)
{

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = i + 1; j < A.size(); j++)
			if (A[i].x >= A[j].x)
			{
				cI = A[i];
				A[i] = A[j];
				A[j] = cI;
			}
		A[i].clus = false;
	}
	return A;
}

vector<CNTInfo> sortV2(vector<CNTInfo> A)
{
	
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = i + 1; j < A.size(); j++)
			if (A[i].y >= A[j].y)
			{
				cI = A[i];
				A[i] = A[j];
				A[j] = cI;
			}
		A[i].clus = false;
	}
	return A;
}

void cluster(int i,int numClusters)
{
	
	for(; i < allCnt.size();i++)
	{
		//переделать условие, сделать условие Если растояние <= с_с, то точно кластер
		//пока что делаем так, что если <= mCh то тоже кластер 
		for (int j = 0; j < clustersInfo[numClusters].size(); j++)
		{
			if(!allCnt[i].clus)
				//if ((sqrt(pow(allCnt[i].k, 2) + pow(radius, 2)) + sqrt(pow(clustersInfo[numClusters][j].k, 2) + pow(radius, 2))) > d(allCnt[i].x, allCnt[i].y, clustersInfo[numClusters][j].x, clustersInfo[numClusters][j].y)) continue;

				if (!vzaim(allCnt[i], clustersInfo[numClusters][j], mCh, mCh))
				{
					//добавляем в кластер
					clustersInfo[numClusters].push_back(allCnt[i]);

					if (firstTest) drawCNT(allCnt[i], colour1, colour2, colour3);
					allCnt[i].clus = true;
					//Sleep(20);
				}
			
		}
	}
}

bool px(CNTInfo c)
{
	if (c.x1 < 0 + mCh ||
		c.x2 < 0 + mCh ||
		c.x3 < 0 + mCh ||
		c.x4 < 0 + mCh) return true;
	return false;
}

bool provx(int i)
{
	bool min = false, max = false;
	for (int j = 0; j < clustersInfo[i].size(); j++)
	{
		if (clustersInfo[i][j].x1 < 0 + mCh ||
			clustersInfo[i][j].x2 < 0 + mCh ||
			clustersInfo[i][j].x3 < 0 + mCh ||
			clustersInfo[i][j].x4 < 0 + mCh) min = true;

		if (clustersInfo[i][j].x1 > L - mCh ||
			clustersInfo[i][j].x2 > L - mCh ||
			clustersInfo[i][j].x3 > L - mCh ||
			clustersInfo[i][j].x4 > L - mCh) max = true;
	}
	if (min&&max) return true;
	return false;
}
bool provy(int i)
{
	bool min = false, max = false;
	for (int j = 0; j < clustersInfo[i].size(); j++)
	{
		if (clustersInfo[i][j].y1 < 0 + mCh ||
			clustersInfo[i][j].y2 < 0 + mCh ||
			clustersInfo[i][j].y3 < 0 + mCh ||
			clustersInfo[i][j].y4 < 0 + mCh) min = true;

		if (clustersInfo[i][j].y1 > L - mCh ||
			clustersInfo[i][j].y2 > L - mCh ||
			clustersInfo[i][j].y3 > L - mCh ||
			clustersInfo[i][j].y4 > L - mCh) max = true;
	}
	if (min&&max) return true;
	return false;
}


void clusters()
{

	fla = true;
	
	allCnt.insert(allCnt.end(), cntInfo.begin(), cntInfo.end());
	allCnt.insert(allCnt.end(), cntTransInfo.begin(), cntTransInfo.end());
	
	allCnt = sortV(allCnt);

	int numClusters = 0;

	for(int i=0; i<allCnt.size(); i++)
	{
		if (allCnt[i].clus) continue;

		clustersInfo.push_back(vector<CNTInfo>(0));
		if (firstTest)
		{
			colour1 = (int)(mt.getReal1() * 200);
			colour2 = (int)(mt.getReal1() * 220);
			colour3 = (int)(mt.getReal1() * 220);

			drawCNT(allCnt[i], colour1, colour2, colour3);
		}
		clustersInfo[numClusters].push_back(allCnt[i]);
		//if (!px(allCnt[i])) return;
		allCnt[i].clus = true;
		//clustersInfo[numClusters][(clustersInfo[numClusters].size() - 1)] = allCnt[0];
		//TempV.pop_back();

		//clustersInfo[numClusters].push_back(allCnt[0]);
		//allCnt.erase(allCnt.begin());
		cluster(i + 1, numClusters);
		if (provx(numClusters))
		{
			kolClus++;
			return;
		}
		numClusters++;
	}

	clustersInfo.clear();
	allCnt = sortV2(allCnt);
	numClusters = 0;

	for (int i = 0; i<allCnt.size(); i++)
	{
		if (allCnt[i].clus) continue;

		clustersInfo.push_back(vector<CNTInfo>(0));
		if (firstTest)
		{
			colour1 = (int)(mt.getReal1() * 200);
			colour2 = (int)(mt.getReal1() * 220);
			colour3 = (int)(mt.getReal1() * 220);

			drawCNT(allCnt[i], colour1, colour2, colour3);
		}
		clustersInfo[numClusters].push_back(CNTInfo(allCnt[i]));
		allCnt[i].clus = true;
		//clustersInfo[numClusters][(clustersInfo[numClusters].size() - 1)] = allCnt[0];
		//TempV.pop_back();

		//clustersInfo[numClusters].push_back(allCnt[0]);
		//allCnt.erase(allCnt.begin());
		cluster(i + 1, numClusters);
		if (provx(numClusters))
		{
			kolClus++;
			return;
		}
		numClusters++;
	}
	clustersInfo.clear();
	allCnt.clear();
}

void main()
{
	CreateDirectoryW(L"files", NULL);
	file.open("./files/coordinates.txt");
	raspr.open("./files/rasp_s.txt");
	dd.open("./files/raspr_a.txt");
	aa.open("./files/test.txt");
	double sh = 0.1;

	//ofstream aa("a.txt"), kk("k.txt");
	setlocale(LC_ALL, "rus");
	cout << "Размер квадрата: ";
	cin >> L;
	cout << "Средняя длина трубки: ";
	cin >> mean;
	cout << "Радиус трубки: ";
	cin >> radius;
	cout << "Начальная плотность: ";
	cin >> p;
	cout << "с шагом: ";
	cin >> sh;
	cout << "Количество испытаний: ";
	cin >> N;
	cout << "Межфазный слой: ";
	cin >> mCh;

	file << "Размер квадрата: " << L << endl;
	file << "Средняя длина трубки: " << mean << endl;
	file << "Радиус: " << radius << endl;
	//file << "Плотность: " << p << endl;
	file << "Количество испытаний: " << N << endl;
	file << "Межчастичный слой: " << mCh << endl;

	mCh = mCh * radius;

	GraphInConsole();
	devi = mean * 0.1;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	transFlag = new bool[8];
	bool*pr = new bool[3];
	pr[0] = false;
	pr[1] = false;
	pr[2] = false;

	for (p; p <= 0.5; p = p + sh)
	{	
		kolClus = 0;
		n = p * L * L / (mean * 2 * radius);
		cout << n << endl;
		for (int i = 0; i < N; i++)
		{
			//100 раз запускаем с p0, считаем кол-во попаданий/100, записываем 
			// делаем так, пока 3 раза подряд не получим 100/100
			
			if (i == 0) firstTest = true;
			else firstTest = false;
			file << "********************************************************************" << endl;
			file << "Испытание" << (i + 1) << endl;
			file << "Плотность: " << p << endl;
			file << "********************************************************************" << endl;
			cout << (i + 1) << "Исп: ";
			double start = clock();
			packaging();
			double finish = clock();
			cout << "уп. за " << (finish - start) / CLOCKS_PER_SEC << "c." << endl;

		//	if (firstTest)
		//		cout << "Примерное время: " << N * (((finish - start) / CLOCKS_PER_SEC) / 60) / 60 << "ч" << endl;

			clusters();

			cout <<"Перколяционный кластер найден: " << kolClus << endl;
			cnt.clear();
			cntTrans.clear();
			cntInfo.clear();
			cntTransInfo.clear();
			clustersInfo.clear();
			allCnt.clear();
		}
		double e = (double)kolClus / N;
		cout << endl << e << endl;
		aa  << e << endl;
		if (e == 1 && !pr[0])pr[0] = true;
		else if (e == 1 && !pr[1]) pr[1] = true; 
			else if (e == 1 && !pr[2]) pr[2] = true; 
		if (pr[0] && pr[1] && pr[2]) break;
	}
	aa.close();
	file.close();
	raspr.close();
	dd.close();
	cout << "Упаковка завершена" << endl;
	cin >> p; 


}


