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

void LDU(double** A, int n, double** L, double** D, double** U)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j < i)
            {
                L[i][j] = A[i][j];
            }
            else if (j > i)
            {
                U[i][j] = A[i][j];
            }
            else
            {
                D[i][j] = A[i][j];
            }
        }
    }
}

void Jacobi(double** A, double** C, double* b, double* y, double* x, double* x_old, int n, double* diff, double** L, double** D, double** U)
{
    int iterations = 0;
    double eps = 0.001;
    //Начальное приближение x^0 любое
    ToNull(x, n);

    //ToOne(C, n);

    for (int i = 0; i < n; i++)
    {
        y[i] = b[i] / A[i][i];
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                C[i][j] = -A[i][j] / A[i][i];
            }
            else
                C[i][j] = 0;
        }
    }

    //x^(k+1)=C*x^k+y
    std::cout << "Матрица С и вектор y: " << std::endl;
    Output(C, n);
    OutputVect(y, n);

    ToNull(L, n);
    ToNull(U, n);
    LDU(C, n, L, D, U);
    std::cout << "Матрица L: " << std::endl;
    Output(L, n);
    std::cout << "Матрица U: " << std::endl;
    Output(U, n);

    std::cout << "Норма матрицы С: " << NormaMat1(C, n) << std::endl;
    Copy(x_old, x, n);
    MultiplyMatrixToVector(C, x_old, x, n);
    SumVect(x, y, n);
    Copy(diff, x, n);
    DiffVect(diff, x_old, n);
    while (NormaVectora1(diff, n) > eps)
    {
        Copy(x_old, x, n);
        MultiplyMatrixToVector(C, x_old, x, n);
        SumVect(x, y, n);
        Copy(diff, x, n);
        DiffVect(diff, x_old, n);
        iterations++;
    }
    std::cout << "Количество иттераций: " << iterations << "\n";
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

    double** L;
    L = new double* [n];
    double** D;
    D = new double* [n];
    double** U;
    U = new double* [n];

    for (int i = 0; i < n; i++)
    {
        C[i] = new double[n];
        A[i] = new double[n];
        L[i] = new double[n];
        D[i] = new double[n];
        U[i] = new double[n];
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
    Jacobi(A, C, b, y, x, x_old, n, diff, L, D, U);
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }
    double* b1;
    b1 = new double[n];

    MultiplyMatrixToVector(A, x, b1, n);
    DiffVect(b1, b, n);
    std::cout << "Норма невязки: " << NormaVectora1(b1, n);
    
    for (int i = 0; i < n; i++)
    {
        delete[] C[i];
        delete[] A[i];
        delete[] L[i];
        delete[] D[i];
        delete[] U[i];
    }

    delete[] L;
    delete[] D;
    delete[] U;
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