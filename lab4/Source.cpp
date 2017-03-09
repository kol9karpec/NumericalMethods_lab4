#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "TSLE.h"

#define ALPHA 0.210
#define BETHA -1.965
#define A 0.675
#define B 1.649
#define C -1.836
#define D 0.983

using namespace std;

double f1I_Iter(double x1, double x2);
double f2I_Iter(double x1, double x2);
double f1I(double x1, double x2);
double f2I(double x1, double x2);

double f1II(double x1, double x2);
double f2II(double x1, double x2);

double norma(double * x, int order);
double ** jakobiFII(double x1, double x2, int order);

typedef double(pfn)(double x1, double x2);
typedef double**(pFJakobi)(double x1, double x2, int order);

double * simpleIteration(double * x0,pfn ** fArr,pfn ** FArr,
	int order, double eps, char ** iterations);	
double * newtonsMethod(double * x0, pFJakobi * J, pfn ** FArr,
	int order, double eps, char ** iterations);


int main(int argc, const char * argv[])
{
	char * inFile = "input.txt";
	char * outFile = "output.txt";

	ifstream in(inFile);
	ofstream out(outFile);

	int order;
	in >> order;

	double * initApprx = new double[order];
	for (int i = 0; i < order; i++)
		in >> initApprx[i];

	double eps;
	in >> eps;

	char ** iterations = new char*[1];

	pfn ** fArr = new pfn*[2];
	fArr[0] = f1I_Iter;
	fArr[1] = f2I_Iter;

	pfn ** FArr = new pfn*[2];
	FArr[0] = f1I;
	FArr[1] = f2I;

	//Simple iterations
	double * root = simpleIteration(initApprx, fArr, FArr, order, eps, iterations);

	out << "Simple iterations" << endl;
	out << "Root: x = (";
	for (int i = 0; i < order; i++)
	{
		out << setprecision(-(int)(log10(eps / 10))) << root[i];
		if (i != order - 1) out << " , ";
	}

	out << ")";
	out << *iterations << endl;		
	delete[] root;

	//Newtons method

	FArr[0] = f1II;
	FArr[1] = f2II;
	root = newtonsMethod(initApprx,jakobiFII,FArr,order,eps,iterations);

	out << "------------------------------------------------" << endl;
	out << "Newtons method" << endl;
	out << "Root: x = (";
	for (int i = 0; i < order; i++)
	{
		out << setprecision(-(int)(log10(eps / 10))) << root[i];
		if (i != order - 1) out << " , ";
	}

	out << ")";
	out << *iterations << endl;

	delete[] initApprx;
	delete[] * iterations;
	delete[] iterations;
	delete[] root;
	delete[] fArr;
	delete[] FArr;

	in.close();
	out.close();


	return 0;
}

double f1I_Iter(double x1, double x2)
{
	return (D - sin(x2 + BETHA));
}

double f2I_Iter(double x1, double x2)
{
	return (C/B - cos(x1 + ALPHA)*(1 / B));
}

double f1I(double x1, double x2)
{
	return (cos(x1+ALPHA)+B*x2-C);
}

double f2I(double x1, double x2)
{
	return (x1+sin(x2+BETHA)-D);
}

double f1II(double x1, double x2)
{
	return (sin(x1+x2)+C*x1-D);
}

double f2II(double x1, double x2)
{
	return (pow(x1,2.0)+pow(x2,2.0)-1);
}

double norma(double * x, int order)
{
	double sum = 0;
	for (int i = 0; i < order; i++)
		sum += pow(x[i], 2.0);

	return sqrt(sum);
}

double ** jakobiFII(double x1, double x2, int order)
{
	double ** res = new double*[order];
	for (int i = 0; i < order; i++)
		res[i] = new double[order];

	res[0][0] = cos(x1 + x2) + C;
	res[0][1] = cos(x1 + x2);
	res[1][0] = 2 * x1;
	res[1][1] = 2 * x2;

	return res;
}

double * simpleIteration(double * x0, pfn ** fArr, pfn ** FArr,
	int order, double eps, char ** iterations)
{
	double * res = new double[order];
	for (int i = 0; i < order; i++)
		res[i] = x0[i];

	int iter = 1;

	*iterations = new char[50000];
	strcpy(*iterations,"");

	while (true)
	{
		double * xPrev = new double[order];
		for (int i = 0; i < order; i++)
			xPrev[i] = res[i];

		for (int i = 0; i < order; i++)
			res[i] = fArr[i](res[0], res[1]);

		double * FRes = new double[order];
		for (int i = 0; i < order; i++)
			FRes[i] = FArr[i](res[0], res[1]);

		double normaF = norma(FRes, order);

		//strcat(*iterations, "Root [simple iteration method]: x=(");
		char * buff = new char[255];
		strcat(*iterations, "\n\nIteration #");
		strcat(*iterations, _itoa(iter, buff, 10));
		strcat(*iterations, "\nxk = (");
		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", res[i]);
			strcat(*iterations, buff);
			if (i != order - 1) strcat(*iterations, " , ");
		}
		strcat(*iterations, ")\nResiduals vector: (");
		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", FRes[i]);
			strcat(*iterations, buff);
			if (i != order - 1) strcat(*iterations, " , ");
		}
		strcat(*iterations, ")\nNorma: ");
		strcpy(buff, "");
		sprintf(buff, "%2.6f", normaF);
		strcat(*iterations, buff);
		delete[] buff;

		double * diffV = new double[order];
		for (int i = 0; i < order; i++)
			diffV[i] = res[i] - xPrev[i];
		double normaDiffV = norma(diffV, order);

		delete[] xPrev;
		delete[] FRes;

		iter++;

		if ((normaF < eps) && (normaDiffV < eps)) break;
	}

	return res;
}

double * newtonsMethod(double * x0, pFJakobi * J, pfn ** FArr,
	int order, double eps, char ** iterations)
{
	double * res = new double[order];
	for (int i = 0; i < order; i++)
		res[i] = x0[i];

	int iter = 1;

	*iterations = new char[50000];
	strcpy(*iterations, "");

	while (true)
	{
		double * xPrev = new double[order];
		for (int i = 0; i < order; i++)
			xPrev[i] = res[i];

		double ** jMatr = J(res[0], res[1], order);

		TSLE dXSystem(jMatr,order);
		double ** b = new double*;
		*b = new double[order];
		for (int i = 0; i < order; i++)
			(*b)[i] = FArr[i](res[0], res[1]);

		double det = 0;
		double ** dX = dXSystem.gauss(b, 1, det, NULL);

		for (int i = 0; i < order; i++)
			res[i] += dX[i][0];

		double * FRes = new double[order];
		for (int i = 0; i < order; i++)
			FRes[i] = FArr[i](res[0], res[1]);

		double normaF = norma(FRes, order);

		//strcat(*iterations, "Root [simple iteration method]: x=(");
		char * buff = new char[255];
		strcat(*iterations, "\n\nIteration #");
		strcat(*iterations, _itoa(iter, buff, 10));
		strcat(*iterations, "\nxk = (");
		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", res[i]);
			strcat(*iterations, buff);
			if (i != order - 1) strcat(*iterations, " , ");
		}
		strcat(*iterations, ")\nResiduals vector: (");
		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", FRes[i]);
			strcat(*iterations, buff);
			if (i != order - 1) strcat(*iterations, " , ");
		}
		strcat(*iterations, ")\nNorma: ");
		strcpy(buff, "");
		sprintf(buff, "%2.6f", normaF);
		strcat(*iterations, buff);
		delete[] buff;

		double * diffV = new double[order];
		for (int i = 0; i < order; i++)
			diffV[i] = res[i] - xPrev[i];
		double normaDiffV = norma(diffV, order);

		delete[] xPrev;
		delete[] FRes;

		iter++;

		if ((normaF < eps) && (normaDiffV < eps)) break;
	}

}


