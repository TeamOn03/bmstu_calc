#include <iostream>
#include <fstream>
#include <cmath>

void transpose(double** matrix, int n) //либо int matrix[][5], либо int (*matrix)[5]
{
    double t;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            t = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = t;
        }
    }
}

void NormVid(int j, int n, double** Mat, double** QMat)
{
    double max = 0;
    int k = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Mat[i][j]) > max)
        {
            k = i;
            max = abs(Mat[i][j]);
        }

    }
    double* temp;
    temp = Mat[j];
    Mat[j] = Mat[k];
    Mat[k] = temp;
    temp = QMat[j];
    QMat[j] = QMat[k];
    QMat[k] = temp;
}

void OutputVect(double* x, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i];
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

void SumMat(double** MatA, double** MatB, double** Result, int n) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Result[i][j] = MatA[i][j] + MatB[i][j];
        }
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

void Copy(double** data, double** newData, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            data[i][j] = newData[i][j];
        }
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

void ToOne(double* OldV, int n)
{
    ToNull(OldV, n);
    for (int j = 0; j < n; j++)
    {
        OldV[j] = 1;
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

double* SumVectWithResult(double* VectA, double* VectB, int n) {
    double* result = new double[n];
    for (int i = 0; i < n; i++) {
        result[i] = VectA[i] + VectB[i];
    }
    return result;
}

void MultiplyMatrix(double** aMatrix, double** bMatrix, double** product, int n)
{
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            product[row][col] = 0;
            for (int inner = 0; inner < n; inner++) {
                product[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
            }
        }
        //std::cout << "\n";
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

void LDU(double* a, double* b, double* c, double* d, int n, double** L, double** D, double** U)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j < i)
            {
                L[i][j] = a[j];
            }
            else if (j > i)
            {
                U[i][j] = c[j];
            }
            else
            {
                D[i][j] = b[j];
            }
        }
    }
}

void MatrixOnCoef(double** Mat, double coef, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Mat[i][j] *= -1;
        }
    }
}

void Zeidel(double** A, double* b, double eps, double* x, double* x_old, double* diff, int n, double** L, double** D, double** U, double** C)
{
    double* y = new double[n];
    int iterations = 0;
    double omega = 1;
    ///*
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
    std::cout << "Матрица С и вектор y: " << std::endl;
    Output(C, n);
    OutputVect(y, n);
    std::cout << "Норма матрицы С: " << NormaMat1(C, n) << std::endl;
    //*/
    do
    {
        ToNull(diff, n);
        Copy(x_old, x, n);
        for (int i = 0; i < n; i++)
        {
            double temp1 = 0, temp2 = 0;

            for (int j = 0; j < n; j++)
            {
                if (j < i)
                {
                    temp1 += A[i][j] / A[i][i] * x[j];
                }
                if (j > i)
                {
                    temp2 += A[i][j] / A[i][i] * x[j];
                }
            }
            x[i] = -omega * temp1 - omega * temp2 + omega * b[i] / A[i][i] + (1 - omega) * x[i];
            diff[i] = x_old[i] - x[i];
            //difference[i] = sol[i] - solution[i];
        }
        iterations++;
    } while (NormaVectora1(diff, n) > eps);
    OutputVect(x, n);

    std::cout << "Количество иттераций: " << iterations << "\n";
    delete[] y;
}


void Zeidel4Vect(double* a, double* b, doublAe* c, double* d, int n, double* x, double* x_old, double* diff, double** L, double** D, double** U, double** C)
{
    double* y = new double[n];
    double omega = 1;
    double eps = 1e-10;
    for (int i = 0; i < n; i++) {
        x_old[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        x[i] = 1;
    }

    for (int i = 0; i < n; i++)
    {
        y[i] = d[i] / b[i];
        for (int j = 0; j < n; j++)
        {
            if (j > i)
            {
                C[i][j] = -c[j] / b[i];
            }
            if (j < i)
            {
                C[i][j] = -a[j] / b[i];
            }
            else
                C[i][j] = 0;
        }
    }
    std::cout << "Норма матрицы С: " << NormaMat1(C, n) << std::endl;
    std::cout << "Матрица С и вектор y: " << std::endl;
    //Output(C, n);
    //OutputVect(y, n);

    Copy(diff, x, n);
    DiffVect(diff, x_old, n);
    double normD = NormaVectora1(diff, n);
    while (normD > eps)
    {
        Copy(x_old, x, n);
        for (int i = 0; i < n; i++) {
            double temp = 0;
            if (i == 0)
                temp += c[i] * x[i + 1];
            if (i == n - 1)
                temp += a[i - 1] * x[i - 1];
            if ((i != 0) && (i != n - 1))
                temp += a[i - 1] * x[i - 1] + c[i] * x[i + 1];
            x[i] = omega * (d[i] - temp) / b[i] + (1 - omega) * x[i];
        }
        Copy(diff, x, n);
        DiffVect(diff, x_old, n);
        normD = NormaVectora1(diff, n);
    }
    delete[] y;
}

int main() {
    std::ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");
    int n;
    fin >> n;
    //std::cout << "Введите размер системы: ";
    //std::cout << "Введите матрицу системы:\n";
    //Переделать матрицы под локальное выделение памяти?
    double* a;
    a = new double[n - 1];

    double* b;
    b = new double[n];

    double* c;
    c = new double[n - 1];

    double* d;
    d = new double[n];

    double* x;
    x = new double[n];

    double* x_old;
    x_old = new double[n];

    double* diff;
    diff = new double[n];

    double** A;
    A = new double* [n];

    double** C;
    C = new double* [n];

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

    ToOne(x, n);
    std::cout << "Решение системы из файла: " << std::endl;
    Zeidel(A, b, 1e-10, x, x_old, diff, n, L, D, U, C);
    double* b1;
    b1 = new double[n];

    MultiplyMatrixToVector(A, x, b1, n);
    DiffVect(b1, b, n);
    std::cout << "Норма невязки: " << NormaVectora1(b1, n);


    n = 216;
    A = new double* [n];

    C = new double* [n];

    L = new double* [n];

    D = new double* [n];

    U = new double* [n];

    for (int i = 0; i < n; i++)
    {
        C[i] = new double[n];
        A[i] = new double[n];
        L[i] = new double[n];
        D[i] = new double[n];
        U[i] = new double[n];
    }
    a = new double[n - 1];

    b = new double[n];

    c = new double[n - 1];

    d = new double[n];

    x = new double[n];

    x_old = new double[n];

    diff = new double[n];
    for (int i = 0; i < n - 1; i++) {
        a[i] = 1;
        c[i] = 1;
        b[i] = 4;
        d[i] = 10 - 2 * ((i + 1) % 2);
    }
    b[n - 1] = 4;
    d[0] = 6;
    d[n - 1] = 9 - 3 * (n % 2);

    std::cout << "Решение заданной системы: " << std::endl;
    Zeidel4Vect(a, b, c, d, n, x, x_old, diff, L, D, U, C);
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }
    //b1 = new double[n];

    //MultiplyMatrixToVector(A, x, b1, n);
    //DiffVect(b1, b, n);
    //std::cout << "Норма невязки: " << NormaVectora1(b1, n);

    delete[] b1;

    for (int i = 0; i < n; i++)
    {
        delete[] C[i];
        delete[] A[i];
        delete[] L[i];
        delete[] D[i];
        delete[] U[i];
    }

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] x;
    delete[] x_old;
    delete[] diff;
    delete[] C;
    delete[] A;
    delete[] L;
    delete[] D;
    delete[] U;
    fin.close();
    return 0;
}