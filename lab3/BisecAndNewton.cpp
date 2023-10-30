#include <iostream>
#include <fstream>
#include <math.h>


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

void MultiplyMatrix(double** aMatrix, double** bMatrix, double** Result, int n)
{
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            Result[row][col] = 0;
            for (int inner = 0; inner < n; inner++) {
                Result[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
            }
        }
        //std::cout << "\n";
    }
}

void MatrixOnCoef(double** Mat, double coef, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Mat[i][j] *= coef;
        }
    }
}

void VectOnCoef(double* Vect, double coef, int n)
{
    for (int i = 0; i < n; i++)
    {
        Vect[i] *= coef;
    }
}

void QR(double** A, double** R, double** Q, double** Qn, int n, double* b, double* bn, double* x)
{
    double eps = 0.0001;
    ToOne(Q, n);
    ToOne(Qn, n);
    double cij;
    double sij;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (abs(A[j][j]) < eps)
                NormVid(j, n, A, Q);
            Copy(R, A, n);
            Copy(Qn, Q, n);
            cij = A[j][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]));
            sij = A[i][j] / (sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]));

            for (int k = 0; k < n; k++)
            {
                R[j][k] = cij * A[j][k] + sij * A[i][k];
                R[i][k] = -sij * A[j][k] + cij * A[i][k];
            }
            for (int k = 0; k < n; k++)
            {
                Qn[j][k] = cij * Q[j][k] + sij * Q[i][k];
                Qn[i][k] = -sij * Q[j][k] + cij * Q[i][k];
            }
            Copy(A, R, n);
            Copy(Q, Qn, n);
        }
    }

    MultiplyMatrixToVector(Q, b, bn, n);
    transpose(Q, n);
    MultiplyMatrix(Q, R, A, n);

    double s;

    //Обратный проход
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += R[i][j] * x[j];
        }
        x[i] = (bn[i] - s) / R[i][i];
    }
}

void Inverse(double** Ainv, double** A, int n)
{
    double* e = new double[n];
    double* bn = new double[n];
    double** R = new double* [n];
    double** Q = new double* [n];
    double** Qn = new double* [n];
    for (int i = 0; i < n; i++)
    {
        R[i] = new double[n];
        Q[i] = new double[n];
        Qn[i] = new double[n];
    }

    for (int i = 0; i < n; i++)
    {
        ToNull(e, n);
        e[i] = 1;
        QR(A, R, Q, Qn, n, e, bn, Ainv[i]);
    }
    transpose(Ainv, n);
    for (int i = 0; i < n; i++)
    {
        delete[] R[i];
        delete[] Q[i];
        delete[] Qn[i];
    }
    delete[] e;
    delete[] bn;
    delete[] R;
    delete[] Q;
    delete[] Qn;
}

void Copy(double data[2], double newData[2])
{
    for (int i = 0; i < 2; i++)
    {
        data[i] = newData[i];
    }
}

double func(double x)
{
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double deriv(double x)
{
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double Bisec(double a, double b)//double(*func)(double))
{
    double acc = 1e-10;
    double x = (a + b) / 2;
    double eps = (b - a) / 2;

    if (func(a) * func(x) <= 0)
    {
        b = x;
    }
    if (func(b) * func(x) <= 0)
    {
        a = x;
    }
    while ((b - a) > acc)
    {
        //std::cout << a << ":a  b:" << b << std::endl;
        x = (a + b) / 2;
        eps = (b - a) / 2;

        if (func(a) * func(x) <= 0)
        {
            b = x;
        }
        //std::cout << func(b) * func(x) << std::endl;
        if (func(b) * func(x) <= 0)
        {
            a = x;
        }
    }
    return x;
}

//Метод Ньютона
//Аналитическое+численное вычисление частных производных

double Newton(double a, double b)//double(*func)(double))
{
    double acc = 1e-3;
    double x = (func(a) * b - func(b) * a) / (func(a) - func(b));

    double memory = x;
    x = x - func(x) * acc / (func(x + acc) - func(x));

    while (abs(x - memory) > acc)
    {
        memory = x;
        x = x - func(x) * acc / (func(x + acc) - func(x));

    }
    return x;
}

double* Newton(double a[2], double b[2])//double(*func)(double))
{
    double acc = 1e-3;
    double x[2];
    x[0] = (func(a[0]) * b[0] - func(b[0]) * a[0]) / (func(a[0]) - func(b[0]));
    x[1] = (func(a[0]) * b[0] - func(b[0]) * a[0]) / (func(a[0]) - func(b[0]));

    double memory[2];
    Copy(memory, x);
    x[0] = x[0] - func(x[0]) * acc / (func(x[0] + acc) - func(x[0]));
    x[1] = x[1] - func(x[1]) * acc / (func(x[1] + acc) - func(x[1]));

    while (abs(x - memory) > acc)
    {
        Copy(memory, x);
        x[0] = x[0] - func(x[0]) * acc / (func(x[0] + acc) - func(x[0]));
        x[1] = x[1] - func(x[1]) * acc / (func(x[1] + acc) - func(x[1]));

    }
    return x;
}



int main() {
    std::ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");

    double NulFunc = Newton(0.56, 0.74);

    std::cout << NulFunc;

    return 0;
}
