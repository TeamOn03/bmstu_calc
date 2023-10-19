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

void matrix_diag_poz_sign(double** matrix, int n, double* b) {
    for (int i = 0; i < n; i++) {
        if (matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
                matrix[i][j] = (-1) * matrix[i][j];
            b[i] *= -1;
        }
    }
}

void SimpleIter(double** A, double** C, double* b, double* y, double* x, double* x_old, int n, double* diff)
{
    int iterations = 0;
    double eps = 1e-9;
    //Начальное приближение x^0 любое
    ToNull(x, n);
    double tao = 0.0001;
    matrix_diag_poz_sign(A, n, b);
    //C = -(tao*A-E)
    //y = tao*b
    ToOne(C, n);
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
    //x^(k+1)=C*x^k+y
    std::cout << "Матрица С и вектор y: " << std::endl;
    Output(C, n);
    OutputVect(y, n);
    std::cout << "Норма матрицы С: " << NormaMat1(C, n) << std::endl;
    Copy(x_old, x, n);
    for (int i = 0; i < n; i++)
    {
        x[i] = tao * b[i];
    }
    Copy(diff, x, n);
    DiffVect(diff, x_old, n);
    while (NormaVectora1(diff, n) > eps)
    {
        Copy(x_old, x, n);
        for (int i = 0; i < n; i++)
        {
            x[i] = 0;
            for (int j = 0; j < n; j++)
            {
                x[i] += -tao * A[i][j] * x_old[j];
                if (i == j)
                    x[i] += x_old[i];
            }
            x[i] += tao * b[i];
        }
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
    SimpleIter(A, C, b, y, x, x_old, n, diff);
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }
    double* b1;
    b1 = new double[n];

    MultiplyMatrixToVector(A, x, b1, n);
    DiffVect(b1, b, n);
    std::cout << "Норма невязки: " << NormaVectora1(b1, n);

    delete[] b1;

    for (int i = 0; i < n; i++)
    {
        delete[] C[i];
        delete[] A[i];
        delete[] L[i];
        delete[] D[i];
        delete[] U[i];
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