#include <iostream>
#include <fstream>
#include <cmath>


void OutputVect(double* A, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << A[i];
        std::cout << ' ';
    }
    std::cout << std::endl;
}

double NormaVectora1(double* b, int n)
{
    double answer = 0;
    for (int i = 0; i < n; i++)
    {
        answer += abs(b[i]);
    }
    return answer;
}


double NormaMat1(double** A, int n)
{
    double maks = 0;
    double DopCount;

    for (int j = 0; j < n; j++)
    {
        DopCount = 0;
        for (int i = 0; i < n; i++)
        {
            DopCount += abs(A[i][j]);
        }
        if (maks < DopCount)
        {
            maks = DopCount;
        }
    }

    return maks;
}

double NormaVectoraInf(double* b, int n)
{
    double maks = 0;
    for (int i = 0; i < n; i++)
    {
        if (maks < abs(b[i]))
        {
            maks = abs(b[i]);
        }
    }
    return maks;
}

double NormaMatInf(double** A, int n)
{
    double maks = 0;
    double DopCount;

    for (int i = 0; i < n; i++)
    {
        DopCount = 0;
        for (int j = 0; j < n; j++)
        {
            DopCount += abs(A[i][j]);
        }
        if (maks < DopCount)
        {
            maks = DopCount;
        }
    }

    return maks;
}

void MultiplyMatrixToVector(double** Matrix, double* vector, double* result, int n) {
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) {
            s += Matrix[i][j] * vector[j];
        }
        result[i] = s;
    }
}

void SumVect(double* VectA, double* VectB, int n) {
    for (int i = 0; i < n; i++)
    {
        VectA[i] += VectB[i];
    }
}

void DiffVect(double* VectA, double* VectB, int n)
{
    for (int i = 0; i < n; i++)
    {
        VectA[i] -= VectB[i];
    }
}

void Copy(double* data, double* newData, int n)
{
    for (int i = 0; i < n; i++)
    {
        data[i] = newData[i];
    }
}

void ToNull(double** OldM, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            OldM[i][j] = 0;
        }
    }
}
void ToNull(double* OldV, int n)
{
    for (int i = 0; i < n; i++)
    {
        OldV[i] = 0;
    }
}

void ToOne(double** OldM, int n)
{
    ToNull(OldM, n);
    for (int j = 0; j < n; j++)
    {
        OldM[j][j] = 1;
    }
}

void Output(double** A, int n) {
    //std::cout << "Матрица:\n" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << A[i][j];
            std::cout << ' ';
        }
        std::cout << std::endl;
    }
}

void SimpleIter(double** A, double** C, double* b, double* y, double* x, double* x_old, int n, double* diff)
{
    double eps = 0.001;
    //Начальное приближение x^0 любое
    ToNull(x, n);
    double tao = 1;
    //C = -(tao*A-E)
    //y = tao*b
    ToOne(C, n);
    while (NormaMat1(C, n) >= 1)
    {
        tao -= 0.0001;
        for (int i = 0; i < n; i++)
        {
            y[i] = tao * b[i];
            for (int j = 0; j < n; j++)
            {
                C[i][j] = -tao * A[i][j];
                if (i == j)
                    C[i][j]++;
            }
        }
        std::cout << NormaMat1(C, n) << " " << NormaMatInf(C, n) << " " << tao << std::endl;
    }
    //x^(k+1)=C*x^k+y
    Copy(x_old, x, n);
    MultiplyMatrixToVector(C, x_old, x, n);
    SumVect(x, y, n);
    Copy(diff, x, n);
    OutputVect(x_old, n);
    DiffVect(diff, x_old, n);
    OutputVect(diff, n);
    while (NormaVectora1(diff, n) > eps)
    {
        Copy(x_old, x, n);
        MultiplyMatrixToVector(C, x_old, x, n);
        SumVect(x, y, n);
        Copy(diff, x, n);
        DiffVect(diff, x_old, n);
    }

}

int main() {
    std::ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");
    int n;
    //std::cout << "Введите размер системы: ";
    fin >> n;
    //std::cout << "Введите матрицу системы:\n";
    //Переделать матрицы под локальное выделение памяти?
    double** C;
    C = new double* [n];

    double** A;
    A = new double* [n];

    double* b;
    b = new double[n];

    double* y;
    y = new double[n];

    double* x;
    x = new double[n];

    double* x_old;
    x_old = new double[n];

    double* diff;
    diff = new double[n];

    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        C[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> A[i][j];
        }
    }

    //std::cout << "Введите свободные коэффициенты:\n";
    for (int i = 0; i < n; i++) {
        fin >> b[i];
    }

    std::cout << "Решение системы из файла: " << std::endl;
    SimpleIter(A, C, b, y, x, x_old, n, diff);
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }

    for (int i = 0; i < n; i++) {
        delete[] A[i];
        delete[] C[i];
    }
    delete[] C;
    delete[] A;
    delete[] b;
    delete[] y;
    delete[] x;
    delete[] x_old;
    delete[] diff;
    fin.close();
    return 0;
}