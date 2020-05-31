// ConsoleApplication7.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include<iostream>
#include<fstream>
#include<math.h>

double res;
using namespace std;
void reverse(double **a, int n) {
	double **b = new double*[n];
	for (int i = 0; i < n; i++) {
		b[i] = new double[n];
		b[i][i] = 1.;
	}
	double **ab = new double*[n];
	for (int i = 0; i < n; i++)
		ab[i] = new double[2 * n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ab[i][j] = a[i][j];
			ab[i][j + n] = b[i][j];
		}
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < 2 * n; i++)
			ab[k][i] = ab[k][i] / a[k][k];
		for (int i = k + 1; i < n; i++) {
			double K = ab[i][k] / ab[k][k];
			for (int j = 0; j < 2 * n; j++)
				ab[i][j] -= ab[k][j] * K;
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = ab[i][j];
	}
	for (int k = n - 1; k > -1; k--) {
		for (int i = 2 * n - 1; i > -1; i--)
			ab[k][i] = ab[k][i] / a[k][k];
		for (int i = k - 1; i > -1; i--) {
			double K = ab[i][k] / ab[k][k];
			for (int j = 2 * n - 1; j > -1; j--)
				ab[i][j] -= ab[k][j] * K;
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i][j] = ab[i][j + n];
}
double cond(double **a, double *x, int n) {
	double cond = 1.;
	double axn = 0.;
	double xn = 0.;
	for (int i = 0; i < n; i++) {
		xn += x[i] * x[i];
		for (int j = 0; j < n; j++)
			axn += pow(a[i][j] * x[j], 2);
	}
	cond *= sqrt(axn) / xn;
	double invaxn = 0.;
	reverse(a, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			invaxn += pow(a[i][j] * x[j], 2);
	cond *= sqrt(invaxn);
	return cond;
}
double norm_of_res(double **a, double *x, double *y, int n) {
	double r[100];
	for (int i = 0; i < n; i++) {
		double s = 0;
		for (int j = 0; j < n; j++)
			s += a[i][j] * x[j];
		r[i] = s;
	}
	double res = 0;

	for (int i = 0; i < n; i++) {
		res += (y[i] - r[i]) * (y[i] - r[i]);
	}
	return sqrt(res);
}
void gauss(double **a, double *b, int m, double *x) {
	for (int k = 1; k < m; k++) {
		for (int j = k; j < m; j++) {
			double h = a[j][k - 1] / a[k - 1][k - 1];
			for (int i = 0; i < m; i++)
				a[j][i] -= h * a[k - 1][i];
			b[j] -= h * b[k - 1];
		}
		for (int i = m - 1; i >= 0; i--) {
			x[i] = b[i] / a[i][i];
			for (int c = m - 1; c > i; c--)
				x[i] -= a[i][c] * x[c] / a[i][i];
		}
	}
}
void PVR(double **a, double *y, int n, int &count, double *x) {
	double *xn;
	xn = new double[10000];
	double norm;
	double w = 1.5;
	double eps = 1e-7;
	for (int i = 0; i < n; i++) {
		xn[i] = 0;
		x[i] = xn[i];
	}
	do {
		count++;
		norm = 0;
		for (int i = 0; i < n; i++) {
			x[i] = y[i];
			for (int j = 0; j < n; j++) {
				if (i != j)
					x[i] = x[i] - a[i][j] * x[j];
			}
			x[i] /= a[i][i];
			x[i] = w * x[i] + (1 - w)*xn[i];
			if (fabs(x[i] - xn[i]) > norm)
				norm = fabs(x[i] - xn[i]);
			xn[i] = x[i];
		}
	} while (norm > eps);
}
void generate(double **a, double *y, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				a[i][j] = 10.;
			else if (j < i + 2)
				a[i][j] = j / (i + j + 2);
			else
				a[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++)
		y[i] = i + 3;
}
int main() {
	const int SIZE = 100;
	double **a = new double*[SIZE];
	for (int i = 0; i < SIZE; i++)
		a[i] = new double[SIZE];
	double *f = new double[SIZE];
	generate(a, f, SIZE);
	double x[SIZE] = {};
	gauss(a, f, SIZE, x);
	ofstream data1;
	data1.open("ans1.dat");
	for (int i = 0; i < SIZE; i++)
		data1 << x[i] << endl;
	data1.close();
	cout << "Gauss: " << endl;
	cout << "	norm of residual: ";
	cout << norm_of_res(a, x, f, SIZE) << endl;
	cout << "	number of condition: ";
	cout << cond(a, x, SIZE) << endl;
	cout << "============================" << endl;
	generate(a, f, SIZE);
	int count = 0;
	PVR(a, f, SIZE, count, x);
	ofstream data2;
	data2.open("ans2.dat");
	for (int i = 0; i < SIZE; i++)
		data2 << x[i] << endl;
	data2.close();
	cout << "PVR: " << endl;
	cout << "	iter: ";
	cout << count << endl;
	cout << "	norm of residual: ";
	cout << norm_of_res(a, x, f, SIZE) << endl;
	cout << "	number of condition: ";
	cout << cond(a, x, SIZE) << endl;
	delete[] f;
	delete[] a;
	return 0;
}