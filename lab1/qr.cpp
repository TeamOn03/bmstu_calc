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

void MultiplyMatrix(double** aMatrix, double** bMatrix, double** product, int n)
{
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            product[row][col] = 0;
            for (int inner = 0; inner < n; inner++) {
                product[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
            }
        }
        std::cout << "\n";
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

void Output(double** A, double* b, int n) {
    //std::cout << "Матрица:\n" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << A[i][j];
            std::cout << ' ';
        }
        std::cout << "   ";
        std::cout << b[i];
        std::cout << std::endl;
    }
}


void NormVid(int j, int n, double** Mat, double** QMat)
{
    double maks = 0;
    int k = 0;
    for (int i = j; i < n; i++)
    {
        if (maks < std::max(maks, abs(Mat[i][j])))
        {
            k = i;
            maks = std::max(maks, abs(Mat[i][j]));
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

double NormaVectora1(double* b, int n)
{
    double answer = 0;
    for (int i = 0; i < n; i++)
    {
        answer += b[i];
    }
    return answer;
}

double* DopCount;

double NormaMat1(double** A, int n)
{
    double maks = 0;
    int k = 0;
    DopCount = new double[n];

    for (int i = 0; i < n; i++)
    {
        DopCount[i] = 0;
    }
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            DopCount[j] += A[i][j];
        }
        if (maks < DopCount[j])
        {
            maks = DopCount[j];
            k = j;
        }
    }

    return DopCount[k];
}

double NormaVectoraInf(double* b, int n)
{
    double maks = 0;
    int k = 0;
    for (int i = 0; i < n; i++)
    {
        if (maks < b[i])
        {
            maks = b[i];
            k = i;
        }
    }
    return b[k];
}

double NormaMatInf(double** A, int n, double* DopCount)
{
    double maks = 0;
    int k = 0;

    for (int i = 0; i < n; i++)
    {
        DopCount[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            DopCount[i] += A[i][j];
        }
        if (maks < DopCount[i])
        {
            maks = DopCount[i];
            k = i;
        }
    }

    return DopCount[k];
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
    //bn = T * b;
    for (int i = 0; i < n; i++)
    {
        bn[i] = 0;
        for (int j = 0; j < n; j++)
        {
            bn[i] += Q[i][j] * b[j];
        }
    }

    MultiplyMatrix(Q, R, A, n);

    transpose(Q, n);
    MultiplyMatrixToVector(Q, b, bn, n);

    double s;

    //Обратный проход
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += R[i][j] * x[j];
        }
        x[i] = (bn[i] - s) / R[i][i];
    }


    //std::cout << "\nОтвет (QR):\n";
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }
}

void Inverse(double** Ainv, double** A, double** R, double** Q, double** Qn, int n, double* e, double* bn)
{
    for (int i = 0; i < n; i++)
    {
        ToNull(e, n);
        e[i] = 1;
        QR(A, R, Q, Qn, n, e, bn, Ainv[i]);
    }
    transpose(Ainv, n);
}


int main() {
    std::ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");
    double eps = 0.0001;
    int n;
    double s;
    std::cout << "Введите размер системы: ";
    fin >> n;
    std::cout << "Введите матрицу системы:\n";
    //Переделать матрицы под локальное выделение памяти?
    double** A;
    A = new double* [n];

    double** Ainv;
    Ainv = new double* [n];

    double** Qn;
    Qn = new double* [n];

    double** Q;
    Q = new double* [n];

    double** R;
    R = new double* [n];

    double* b;
    b = new double[n];

    double* bn;
    bn = new double[n];

    double* x;
    x = new double[n];

    double* temp;
    temp = new double[n];

    double* DopCount;
    DopCount = new double[n];

    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        Ainv[i] = new double[n];
        Q[i] = new double[n];
        R[i] = new double[n];
        Qn[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> A[i][j];
        }
    }

    std::cout << "Введите свободные коэффициенты:\n";
    for (int i = 0; i < n; i++) {
        fin >> b[i];
    }

    std::cout << "Решение системы:";
    QR(A, R, Q, Qn, n, b, bn, x);

    //3
    double* b1;
    b1 = new double[n];

    double bdop;

    for (int i = 0; i < n; i++)
    {
        bdop = 0;
        for (int j = 0; j < n; j++)
        {
            bdop += A[i][j] * x[j];
        }
        b1[i] = bdop;
        for (int j = 0; j < n; j++)
        {
            b1[i] = b[i] - b1[i];
        }
    }
    std::cout << "Норма невязки:" << (NormaVectora1(b1, n));
    delete[] b1;
    /////////////////////////////////////////////////////////

    //4,6
    double* e;
    e = new double[n];


    Inverse(Ainv, A, R, Q, Qn, n, e, bn);
    std::cout << "Обратная матрица:" << std::endl;
    Output(Ainv, x, n);

    std::cout << "Число обусловленности: " << NormaMatInf(A, n, DopCount) * NormaMatInf(Ainv, n, DopCount);
    delete[] e;
    //+оценка снизу с помощью полученных решений

    //5
    double vozm = 0.01;
    for (int i = 0; i < n; i++) {
        b[i] = b[i] + vozm;
    }

    double* x1;
    x1 = new double[n];

    QR(A, R, Q, Qn, n, b, bn, x1);

    for (int i = 0; i < n; i++)
    {
        x1[i] = abs(x[i]) - abs(x1[i]);
        std::cout << x1[i] << std::endl;
    }

    std::cout << "Норма вектора расхождений решений" << NormaVectora1(x1, n);

    delete[] x1;
    ////////////////////////////////////////////////////////////////////////////

    std::cout << "Единичная:";
    MultiplyMatrix(A, Ainv, Q, n);
    Output(Q, b, n);

    for (int i = 0; i < n; i++) {
        delete[] A[i];
        delete[] Ainv[i];
        delete[] Q[i];
        delete[] R[i];
        delete[] Qn[i];
    }
    delete[] Q;
    delete[] Qn;
    delete[] R;
    delete[] A;
    delete[] Ainv;
    delete[] b;
    delete[] bn;
    delete[] x;
    delete[] temp;
    delete[] DopCount;
    fin.close();
    return 0;
}