#include <fstream>
#include <iostream>
#include <cmath>
#include "mersennetwister.h"
#include <Windows.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include "cnt.h"
#define M_PI 3.14159265358979323846

using namespace std;

//************************************************************************
int L, N, n;
double mean, devi, radius = 1;

int colour1 = 0, colour2 = 0, colour3 = 0;
ofstream file, raspr, flog, aa;
bool flag, fla=false;
bool *transFlag;
MtRng64 mt;
bool ready = false;
double second = 0.0;
int changeX = 350, changeY=50;

double mF = 0;
vector <cnt> cntList(0);
cnt cI1, cI2;
bool *visited, draw = false;


//************************************************************************
void GraphInConsole()
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
	
	DeleteObject(Pen);
	ReleaseDC(GetConsoleWindow(), hDC);

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

cnt cntWithMF(cnt c, int mf)
{
	return cnt(coordX(c.x, mf, c.a - 180), coordY(c.y, mf, c.a - 180), c.k + mf * 2, c.a, radius + mf);
}

void drawCNT(cnt loc, int i, int j, int k)
{
	if (draw)
	{
		HDC hDC = GetDC(GetConsoleWindow());
		HPEN Pen = CreatePen(PS_SOLID, 2, RGB(i, j, k));
		SelectObject(hDC, Pen);

		MoveToEx(hDC, loc.x1 + changeX, loc.y1 + changeY, NULL);
		LineTo(hDC, loc.x4 + changeX, loc.y4 + changeY);

		MoveToEx(hDC, loc.x2 + changeX, loc.y2 + changeY, NULL);
		LineTo(hDC, loc.x3 + changeX, loc.y3 + changeY);

		MoveToEx(hDC, loc.x1 + changeX, loc.y1 + changeY, NULL);
		LineTo(hDC, loc.x2 + changeX, loc.y2 + changeY);

		MoveToEx(hDC, loc.x4 + changeX, loc.y4 + changeY, NULL);
		LineTo(hDC, loc.x3 + changeX, loc.y3 + changeY);

		DeleteObject(Pen);
		ReleaseDC(GetConsoleWindow(), hDC);
	}

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
bool belong(double x, double y, cnt c) //true - точка (x, y) лежит внутри трубки или недопустимо близко
{
	if ((c.lA()*x + c.lB()*y + c.lC())*(c.rA()*x + c.rB()*y + c.rC()) > 0 &&
		(c.bA()*x + c.bB()*y + c.bC())*(c.tA()*x + c.tB()*y + c.tC()) > 0) return true;
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

bool vzaim(cnt c1, cnt c2, double mF) //true - не пересекаются
{
	cI1 = cntWithMF(c1, mF);
	cI2 = cntWithMF(c2, mF);


	if (belong(cI1.x1, cI1.y1, cI2)) return false;
	if (belong(cI1.x2, cI1.y2, cI2)) return false;
	if (belong(cI1.x3, cI1.y3, cI2)) return false;
	if (belong(cI1.x4, cI1.y4, cI2)) return false;

	if (belong(coordX(cI1.x, cI1.k / 2, cI1.a), coordY(cI1.y, cI1.k / 2, cI1.a), cI2)) return false;
	if (belong(coordX(cI2.x, cI2.k / 2, cI2.a), coordY(cI2.y, cI2.k / 2, cI2.a), cI1)) return false;

	return true;
}

bool vzaim(cnt cntNew, cnt locInfo) // true - не пересекаются
{
	if (belong(cntNew.x1, cntNew.y1, locInfo)) return false;
	if (belong(cntNew.x2, cntNew.y2, locInfo)) return false;
	if (belong(cntNew.x3, cntNew.y3, locInfo)) return false;
	if (belong(cntNew.x4, cntNew.y4, locInfo)) return false;

	return true;
}
bool allTest(cnt cntNew, vector<cnt>loc) //true - удачное расположение
{
	for (int i = 0; i < loc.size(); i++)
	{
		if ((loc[i].k+radius)*2.0 < d(cntNew.x, cntNew.y, loc[i].x, loc[i].y)) continue;

		if (!vzaim(cntNew, loc[i])) return false;
		if (!vzaim(loc[i], cntNew)) return false;

		if (check(cntNew.x1, cntNew.y1, loc[i].x1, loc[i].y1, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false;/*left-left*/
		if (check(cntNew.x1, cntNew.y1, loc[i].x2, loc[i].y2, cntNew.k, cntNew.a, loc[i].k, loc[i].a)) return false; /*left-right*/

	}
	return true;
}

bool test(double x, double y, double k, int a) //true - удачное расположение
{

	cnt cIM = cnt(x, y, k, a, radius);

	if (!allTest(cIM, cntList)) return false;

	if (cIM.x1 < 0 || cIM.x2 < 0 || cIM.x3 < 0 || cIM.x4 < 0)
	{
		if (!allTest(cnt(x + L, y, k, a, radius), cntList)) return false;
		transFlag[0] = true;
	}
	if (cIM.x1 > L || cIM.x2 > L || cIM.x3 > L || cIM.x4 > L)
	{
		if (!allTest(cnt(x - L, y, k, a, radius), cntList)) return false;
		transFlag[1] = true;
	}
	if (cIM.y1 < 0 || cIM.y2 < 0 || cIM.y3 < 0 || cIM.y4 < 0)
	{
		if (!allTest(cnt(x, y + L, k, a, radius), cntList)) return false;
		transFlag[2] = true;
	}
	if (cIM.y1 > L || cIM.y2 > L || cIM.y3 > L || cIM.y4 > L)
	{
		if (!allTest(cnt(x, y - L, k, a, radius), cntList)) return false;
		transFlag[3] = true;
	}

	if (cIM.x1 < 0 && cIM.y1 < 0 || cIM.x2 < 0 && cIM.y2 < 0 || cIM.x3 < 0 && cIM.y3 < 0 || cIM.x4 < 0 && cIM.y4 < 0)
	{
		if (!allTest(cnt(x + L, y + L, k, a, radius), cntList)) return false;
		transFlag[4] = true;
	}
	if (cIM.x1 > L && cIM.y1 < 0 || cIM.x2 > L && cIM.y2 < 0 || cIM.x3 > L && cIM.y3 < 0 || cIM.x4 > L && cIM.y4 < 0)
	{
		if (!allTest(cnt(x - L, y + L, k, a, radius), cntList)) return false;
		transFlag[5] = true;
	}

	if (cIM.x1 > L && cIM.y1 > L || cIM.x2 > L && cIM.y2 > L || cIM.x3 > L && cIM.y3 > L || cIM.x4 > L && cIM.y4 > L)
	{
		if (!allTest(cnt(x - L, y - L, k, a, radius), cntList)) return false;
		transFlag[6] = true;
	}
	if (cIM.x1 < 0 && cIM.y1 > L || cIM.x2 < 0 && cIM.y2 > L || cIM.x3 < 0 && cIM.y3 > L || cIM.x4 < 0 && cIM.y4 > L)
	{
		if (!allTest(cnt(x + L, y - L, k, a, radius), cntList)) return false;
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

void createMatrix(int i, int **m) //true - удачное расположение
{
	for (int j = 0; j < cntList.size(); j++)
	{
		if ((cntList[i].k + cntList[j].k + 2.0 * mF + 2.0) < d(cntList[j].x, cntList[j].y, cntList[i].x, cntList[i].y)) continue;

		if (!vzaim(cntList[i], cntList[j], mF))
		{
			m[i][j] = 1;
			m[j][i] = 1;
		}
		if (!vzaim(cntList[j], cntList[i], mF))
		{
			m[i][j] = 1;
			m[j][i] = 1;
		}
	}
}

void addTransCNT(double x, double y, int a, double k, int idParent, int **m)
{

	//parent = 0 - не дочерний, не родитель
	//parent = 1 - родитель
	//parent = 2 - дочерний

	cntList.push_back(cnt(x, y, k, a, radius, idParent, 2));
	cntList[idParent].parent = 1;

	int i = cntList.size() - 1;

	drawCNT(cntList[i], 225, 0, 0);
	createMatrix(i, m);

	m[idParent][i] = 2;
	m[i][idParent] = 2;
}
void packaging(int **m)
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
				cout << " n-- " << endl;
				n--;
				return;
			}
			kol++;
		} while (!test(x, y, k, a));

		if (flag)
		{
			cntList.push_back(cnt(x, y, k, a, radius));
			int id = cntList.size() - 1;
			createMatrix(id, m);

			drawCNT(cntList[id], 225, 225, 225);

			if (transFlag[0]) addTransCNT(x + L, y, a, k, id, m);
			if (transFlag[1]) addTransCNT(x - L, y, a, k, id, m);
			if (transFlag[2]) addTransCNT(x, y + L, a, k, id, m);
			if (transFlag[3]) addTransCNT(x, y - L, a, k, id, m);
			if (transFlag[4]) addTransCNT(x + L, y + L, a, k, id, m);
			if (transFlag[5]) addTransCNT(x - L, y + L, a, k, id, m);
			if (transFlag[6]) addTransCNT(x - L, y - L, a, k, id, m);
			if (transFlag[7]) addTransCNT(x + L, y - L, a, k, id, m);


			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << setw(7) << a << "|" << endl;
			S += k * radius * 2.0;
		}
	}
	file << "Реальная плотность: " << S / (L*L) << endl;
	//flog << "реальная p: " << S / (L*L) << " ";
}


void DFS(int st, int colour1, int colour2, int colour3, int clusters, int **m) // кластеризация
{
	visited[st] = true;
	drawCNT(cntList[st], colour1, colour2, colour3);
	cntList[st].idClus = clusters;

	for (int r = 0; r < cntList.size(); r++)
	{
		if ((m[st][r] != 0) && (!visited[r]))
			DFS(r, colour1, colour2, colour3, clusters, m);
	}
}
bool px0(cnt c)
{
	if (c.x1 < 0 + mF || c.x2 < 0 + mF || c.x3 < 0 + mF || c.x4 < 0 + mF) return true;
	return false;
}
bool pxL(cnt c)
{
	if (c.x1 > L - mF || c.x2 > L - mF || c.x3 > L - mF || c.x4 > L - mF) return true;
	return false;
}
bool py0(cnt c)
{
	if (c.y1 < 0 + mF || c.y2 < 0 + mF || c.y3 < 0 + mF || c.y4 < 0 + mF) return true;
	return false;
}
bool pyL(cnt c)
{
	if (c.y1 > L - mF || c.y2 > L - mF || c.y3 > L - mF || c.y4 > L - mF) return true;
	return false;
}

bool ifParent(int st, int r)
{
	if (st != cntList[r].idParent && r != cntList[st].idParent) return true;
	return false;
}

void DFS(int st, bool &pc, bool x, int **m) //x = true - по оси x,  = false - по оси y
{
	visited[st] = true;
	if (x)
	{
		if (pxL(cntList[st])) pc = true;
	}
	else
	{
		if (pyL(cntList[st])) pc = true;
	}

	drawCNT(cntWithMF(cntList[st], 1), 225, 225, 225);
	for (int r = 0; r < cntList.size(); r++)
		if (m[st][r] !=0 && !visited[r] && ifParent(st, r))
			DFS(r, pc, x, m);
}

void percolationClusters(vector <int> &pClus, vector <int> locVector, bool coord, int **m) //coord = true - по x, coord = false - по y
{
	for (int i = 0; i < locVector.size(); i++)
	{
		bool f = false;
		for (int j = 0; j < pClus.size(); j++)
			if (cntList[locVector[i]].idClus == cntList[pClus[j]].idClus) f = true; //этот перколяционный кластер уже найден
		if (f) continue;

		for (int j = 0; j < cntList.size(); j++)
			visited[j] = false;

		bool pc = false;
		DFS(locVector[i], pc, coord, m);
		if (pc) // если кластер перколяционный, то добавляем в список id и рисуем
		{
			for (int j = 0; j < cntList.size(); j++)
				if (cntList[j].idClus == cntList[locVector[i]].idClus) drawCNT(cntWithMF(cntList[j], 1), 225, 0, 0);
			pClus.push_back(locVector[i]);
		}
	}
}
void main()
{
	CreateDirectoryW(L"files", NULL);
	file.open("./files/coordinates.txt");
	//raspr.open("./files/rasp_s.txt");
	flog.open("./files/log.txt");
	aa.open("./files/значения вероятностей.txt");

	setlocale(LC_ALL, "rus");
	cout << "Размер квадрата: ";
	cin >> L;
	cout << "Средняя длина трубки: ";
	cin >> mean; 
	cout << "Количество испытаний: ";
	cin >> N;
	/*cout << "Концентрация: ";
	cin >> p;*/
	cout << "Проницаемый слой: ";
	cin >> mF;
	cout << "Рисовать(+/-): ";
	char strDraw;
	cin >> strDraw;

	if (strDraw == '+') draw = true;
	else draw = false;

	L *= radius;
	mean *= radius;
	mF *= radius;

	file << "Размер квадрата: " << L << endl;
	file << "Средняя длина трубки: " << mean << endl;
	file << "Радиус: " << radius << endl;
	file << "Количество испытаний: " << N << endl;
	file << "Проницаемый слой: " << mF << endl;

	devi = mean * 0.1;
	file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;

	transFlag = new bool[8];
	bool*pr = new bool[3];
	for (int i = 0; i < 3; i++)
		pr[i] = false;

	bool stop = false;
	aa << setw(4) << "p" << setw(5) << "вер." << endl;
	int pCl = 0;
	for (double p = 0; p <= 0.20; p +=0.005)
	{	
		pr[0] = pr[1];
		pr[1] = pr[2];

		n = p * L * L / (mean * 2 * radius); //количество к-меров 
		cout<<"Количество: " << n << " p = " << p <<  endl;
		double fullstart = clock();
		int sizeM = n*1.5;

		int **m = new int*[sizeM];
		for (int w = 0; w < sizeM; w++)
			m[w] = new int[sizeM];
		
		pCl = 0;

		for (int i = 0; i < N; i++)
		{
			double start = clock();
			//N раз запускаем с p, считаем кол-во попаданий/N, записываем 
			// делаем так, пока 3 раза подряд не получим 1
			
			file << "********************************************************************" << endl;
			file << "Испытание " << (i + 1) << endl;
			file << "Плотность: " << p << endl;
			file << "Количество: " << n << endl;
			file << "********************************************************************" << endl;

			//успаковка 
			for (int w = 0; w < sizeM; w++)
				for (int q = 0; q < sizeM; q++)
					m[w][q] = 0;

			packaging(m);

			//поиск всех кластеров 
			int clusters = 0;
			visited = new bool[cntList.size()];
			for (int w = 0; w<cntList.size(); w++)
				visited[w] = false;
			for (int q = 0; q < cntList.size(); q++)
			{
				if (cntList[q].idClus == 0)
				{
					clusters++;
					colour1 = (int)(mt.getReal1() * 220);
					colour2 = (int)(mt.getReal1() * 220);
					colour3 = (int)(mt.getReal1() * 220);
					DFS(q, colour1, colour2, colour3, clusters, m);
				}
			}

			//поиск перколяционных кластеров
			vector <int> x0(0), y0(0), pClusters(0);
			for (int w = 0; w < cntList.size(); w++)
			{
				if (px0(cntList[w])) x0.push_back(w);
				if (py0(cntList[w])) y0.push_back(w);
			}

			percolationClusters(pClusters, x0, true, m);
			percolationClusters(pClusters, y0, false, m);

			//flog <<"Всего кластеров найдено: " << clusters << endl;
			//flog << "Из них перколяционных: " << pClusters.size() << endl;
			cout << i+1 << "исп. perClus=" << pClusters.size() << " ";
			if (pClusters.size() != 0) pCl++;

			cntList.clear();
			pClusters.clear();
			x0.clear();
			y0.clear();

			if (draw)
			{
				HWND hwnd = GetConsoleWindow(); //Берём ориентир на консольное окно (В нём будем рисовать)
				RECT WinCoord = {}; //Массив координат окна  
				GetWindowRect(hwnd, &WinCoord); //Узнаём координаты
				HDC dc = GetDC(hwnd);

				HBRUSH brush = CreateSolidBrush(RGB(0, 0, 0)); //Создаём кисть определённого стиля и цвета
				SelectObject(dc, brush); //Выбираем кисть
				Rectangle(dc, changeX - 100, changeY - 10, WinCoord.right, 2 * L); //Рисуем новый прямоугольник, который будет небом
				DeleteObject(brush); //Очищаем память от созданной, но уже ненужной кисти
				ReleaseDC(hwnd, dc);
				GraphInConsole();
			}

			double finish = clock();
			cout << (finish - start) / CLOCKS_PER_SEC << "sec"<<endl;
		}
 
		double e = (double)pCl / N;
		if (e == 1.0) pr[2] = true;
		aa << setw(4) << p << setw(5) << e << endl;
		cout << "Вероятность: " << e << endl;
		for (int i = 0; i < sizeM; i++)
			delete[]m[i];
		delete[]m;
		double fullfinish = clock();
		//cout<<(fullfinish - fullstart) / CLOCKS_PER_SEC/60 << "m"<<endl;
		if (pr[0] && pr[1] && pr[2]) break;
	}
	cntList.clear();
	aa.close();
	file.close();
	raspr.close();
	flog.close();
	cout << "Упаковка завершена" << endl;
	system("pause");
}


