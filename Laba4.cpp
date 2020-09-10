#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
using namespace std;

const double eps = 0.001;

/////////////////////////////////////////////////////////////////////////////////////////////////

double Residual_vect(double** A, double* f, double* metod, int n)
{
	double* x0;
	double rate = 0;
	x0 = new double[n];
	for (int i = 0; i < n; i++)
		x0[i] = 0.0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			x0[i] += A[i][j] * metod[j];
		}
		x0[i] = f[i] - x0[i];
	}
	for (int i = 0; i < n; i++)
	{
		rate += (x0[i] * x0[i]);
	}
	return sqrt(rate);
}


/////////////////////////////////////////////////////////////////////////////////////////////////


//Метод Гауса
double* GausMet(double** A, double* F, int N)
{
	double* x, Max;
	int k_count, ind;
	x = new double[N];
	k_count = 0;
	while (k_count < N)
	{
		Max = abs(A[k_count][k_count]);
		ind = k_count;
		for (int i = k_count + 1; i < N; i++)
		{
			if (abs(A[i][k_count]) > Max)
			{
				Max = abs(A[i][k_count]);
				ind = i;
			}
		}
		for (int j = 0; j < N; j++)
		{
			double temp = A[k_count][j];
			A[k_count][j] = A[ind][j];
			A[ind][j] = temp;
		}
		double temp = F[k_count];
		F[k_count] = F[ind];
		F[ind] = temp;

		for (int i = k_count; i < N; i++)
		{
			double temp = A[i][k_count];
			if (abs(temp) < eps)
				continue;
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] / temp;
			F[i] = F[i] / temp;
			if (i == k_count)
				continue;
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] - A[k_count][j];
			F[i] = F[i] - F[k_count];
		}
		k_count++;
	}
	for (k_count = N - 1; k_count >= 0; k_count--)
	{
		x[k_count] = F[k_count];
		for (int i = 0; i < k_count; i++)
			F[i] = F[i] - A[i][k_count] * x[k_count];
	}
	return x;
}


/////////////////////////////////////////////////////////////////////////////////////////////////

//метод ПВР
double* PVR(double** A, double* F, int n)
{
	int k = 0;
	const double eps = 0.00001;
	double w, norma;
	double* xn = 0;
	double* x;
	x = new double[n];
	xn = new double[n];
	w = 1.000001;
	for (int i = 0; i < 10; i++)
	{
		xn[i] = 0;
		x[i] = xn[i];
	}
	do
	{
		k++;
		norma = 0;
		for (int i = 0; i < n; i++)
		{
			x[i] = F[i];
			for (int j = 0; j < n; j++)
			{
				if (i != j)
					x[i] = x[i] - A[i][j] * x[j];
			}
			x[i] /= A[i][i];

			x[i] = w * x[i] + (1 - w) * xn[i];

			if (fabs(x[i] - xn[i]) > norma)
				norma = fabs(x[i] - xn[i]);
			xn[i] = x[i];
		}
	} while (norma > eps);
	cout << "Number of iterations: "; cout << k << " \n";

	return x;
}

/////////////////////////////////////////////////////////////////////////////////////////////////


//очистить выделенную память
void clear(double** arr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] arr[i];
	delete[] arr;
}
//создать копию массива
double** clone(double** arr, int n)
{
	double** newArr = new double*[n];
	for (int row = 0; row < n; row++)
	{
		newArr[row] = new double[n];
		for (int col = 0; col < n; col++)
			newArr[row][col] = arr[row][col];
	}
	return newArr;
}
//матричное умножение матриц
double** matrix_multi(double** A, double** B, int n)
{
	double** result = new double*[n];
	//заполнение нулями
	for (int row = 0; row < n; row++)
	{
		result[row] = new double[n];
		for (int col = 0; col < n; col++)
		{
			result[row][col] = 0.0;
		}
	}

	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			for (int j = 0; j < n; j++)
			{
				result[row][col] += A[row][j] * B[j][col];
			}
		}
	}
	return result;
}
//умножение матрицы на число
void scalar_multi(double** m, int n, double a)
{
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
		{
			m[row][col] *= a;
		}
}
//вычисление суммы двух квадратных матриц
void sum(double** A, double** B, int n)
{
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			A[row][col] += B[row][col];
}

//вычисление определителя
double det(double** matrix, int n) //квадратная матрица размера n*n
{
	double** B = clone(matrix, n);
	//приведение матрицы к верхнетреугольному виду
	for (int step = 0; step < n - 1; step++)
		for (int row = step + 1; row < n; row++)
		{
			double coeff = -B[row][step] / B[step][step]; //метод Гаусса
			for (int col = step; col < n; col++)
				B[row][col] += B[step][col] * coeff;
		}
	//Рассчитать определитель как произведение элементов главной диагонали
	double Det = 1;
	for (int i = 0; i < n; i++)
		Det *= B[i][i];
	//Очистить память
	clear(B, n);
	return Det;
}

void invert(double** M, int n)
{
	//Исходная матрица, динамический двухмерный массив
	double** A = new double*[n];
	for (int row = 0; row < n; row++)
	{
		A[row] = new double[n];
		for (int col = 0; col < n; col++)
			A[row][col] = M[row][col];
	}


	double N1 = 0, Ninf = 0; //норма матрицы по столбцам и по строкам
	double** A0 = clone(A, n);       //инициализация начального приближения
	for (size_t row = 0; row < n; row++) {
		double colsum = 0, rowsum = 0;
		for (size_t col = 0; col < n; col++) {
			rowsum += fabs(A0[row][col]);
			colsum += fabs(A0[col][row]);
		}
		N1 = max(colsum, N1);
		Ninf = max(rowsum, Ninf);
	}
	//транспонирование
	for (size_t row = 0; row < n - 1; row++) {
		for (size_t col = row + 1; col < n; col++)
			swap(A0[col][row], A0[row][col]);
	}
	scalar_multi(A0, n, (1 / (N1 * Ninf))); //нормирование матрицы
											//инициализация удвоенной единичной матрицы нужного размера
	double** E2 = new double*[n];
	for (int row = 0; row < n; row++)
	{
		E2[row] = new double[n];
		for (int col = 0; col < n; col++)
		{
			if (row == col)
				E2[row][col] = 2;
			else
				E2[row][col] = 0;
		}
	}
	double** inv = clone(A0, n); 
	double EPS = 0.001;   //погрешность
	if (det(A, n) != 0) //если матрица не вырождена
	{
		while (fabs(det(matrix_multi(A, inv, n), n) - 1) >= EPS) 
		{
			double** prev = clone(inv, n); 
			inv = matrix_multi(A, prev, n);  
			scalar_multi(inv, n, -1);        
			sum(inv, E2, n);                  
			inv = matrix_multi(prev, inv, n); 
			clear(prev, n);
		}

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				M[i][j] = inv[i][j];
	}
	else
		cout << "Impossible\n";
	clear(A, n);
	clear(E2, n);
}

////////////////////////////////////////////////////////

int main()
{
	double ** Coef, ** CoefN, ** CoefG, ** Coefpvr, *fG, *fpvr, *f, *Gaus_Ans, *PVR_Ans, Max_1 = -256, Max_2 = -256, Sum = 0;
	int n = 100;
	int b = 10;
	Coef = new double*[n];
	CoefG = new double*[n];
	CoefN = new double*[n];
	Coefpvr = new double*[n];
	f = new double[n];
	fG = new double[n];
	fpvr = new double[n];
	for (int i = 0; i < n; i++)
	{
		Coef[i] = new double[n];
		CoefN[i] = new double[n];
		CoefG[i] = new double[n];
		Coefpvr[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			if (j > i + 1)
				Coef[i][j] = 0;
			else if (!(i == j))
				Coef[i][j] = (1.0 * (j + 1) / ((i + 1) + (j + 1)));
			else if (i == j)
				Coef[i][j] = b;

			CoefG[i][j] = Coef[i][j];
			Coefpvr[i][j] = Coef[i][j];
			CoefN[i][j] = Coef[i][j];
		}
	}
	for (int i = 0; i < n; i++)
	{
		f[i] = ((i + 1) + 2);
		fG[i] = f[i];
		fpvr[i] = f[i];
	}



	
	Gaus_Ans = GausMet(CoefG, fG, n);

	ofstream gauss("ans1.dat");
	for (int i = 0; i < n; i++)
	{
		gauss << Gaus_Ans[i] << endl;
	}
	gauss.close();

	PVR_Ans = PVR(Coefpvr, fpvr, n);

	ofstream pvr("ans2.dat");
	for (int i = 0; i < n; i++)
	{
		pvr << PVR_Ans[i] << endl;
	}
	pvr.close();

	for (int i = 0; i < n; i++)
	{
		cout << "x(PVR)[" << i + 1 << "] = " << PVR_Ans[i] << "\t\t" << "x(GAUSS)[" << i + 1 << "] = " << Gaus_Ans[i] << endl;
	}
	/////////////////////////////////////////////////////////////////////////////////
	cout << "Norm of the discrepancy vector for Gauss metod: " << Residual_vect(CoefN, f, Gaus_Ans, n) << endl;
	cout << "Norm of the discrepancy vector for PVR metod: " << Residual_vect(CoefN, f, PVR_Ans, n) << endl;

	/////////////////////////////////////////////////////////////////////////////////

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
			Sum += abs(Coef[i][j]);
		if (Max_1 < Sum)
			Max_1 = Sum;
		Sum = 0;
	}
	invert(Coef, n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
			Sum += abs(Coef[i][j]);
		if (Max_2 < Sum)
			Max_2 = Sum;
		Sum = 0;
	}
	cout << "condition number of the matrix: " << Max_1 * Max_2 << endl;



	return 0;
}
