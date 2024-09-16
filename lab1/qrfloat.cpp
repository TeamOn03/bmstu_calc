#include <iostream>
#include <fstream>
#include <cmath>


void transpose(float** matrix, int n) //либо int matrix[][5], либо int (*matrix)[5]
{
    float t;
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

void MultiplyMatrix(float** aMatrix, float** bMatrix, float** product, int n)
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

void Copy(float** data, float** newData, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            data[i][j] = newData[i][j];
        }
    }
}

void ToNull(float** OldM, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            OldM[i][j] = 0;
        }
    }
}
void ToNull(float* OldV, int n)
{
    for (int i = 0; i < n; i++)
    {
        OldV[i] = 0;
    }
}

void ToOne(float** OldM, int n)
{
    ToNull(OldM, n);
    for (int j = 0; j < n; j++)
    {
        OldM[j][j] = 1;
    }
}

void Output(float** A, float* b, int n) {
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


void NormVid(int j, int n, float** Mat, float** QMat)
{
    float max = 0;
    int k = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Mat[i][j]) > max)
        {
            k = i;
            max = abs(Mat[i][j]);
        }

    }
    float* temp;
    temp = Mat[j];
    Mat[j] = Mat[k];
    Mat[k] = temp;
    temp = QMat[j];
    QMat[j] = QMat[k];
    QMat[k] = temp;
}

float NormaVectora1(float* b, int n)
{
    float answer = 0;
    for (int i = 0; i < n; i++)
    {
        answer += abs(b[i]);
    }
    return answer;
}

float* DopCount;

float NormaMat1(float** A, int n, float* DopCount)
{
    float maks = 0;
    DopCount = new float[n];

    for (int i = 0; i < n; i++)
    {
        DopCount[i] = 0;
    }
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            DopCount[j] += abs(A[i][j]);
        }
        if (maks < DopCount[j])
        {
            maks = DopCount[j];
        }
    }

    return maks;
}

float NormaVectoraInf(float* b, int n)
{
    float maks = 0;
    for (int i = 0; i < n; i++)
    {
        if (maks < abs(b[i]))
        {
            maks = abs(b[i]);
        }
    }
    return maks;
}

float NormaMatInf(float** A, int n, float* DopCount)
{
    float maks = 0;

    for (int i = 0; i < n; i++)
    {
        DopCount[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            DopCount[i] += abs(A[i][j]);
        }
        if (maks < DopCount[i])
        {
            maks = DopCount[i];
        }
    }

    return maks;
}

void MultiplyMatrixToVector(float** Matrix, float* vector, float* result, int n) {
    for (int i = 0; i < n; i++) {
        float s = 0;
        for (int j = 0; j < n; j++) {
            s += Matrix[i][j] * vector[j];
        }
        result[i] = s;
    }
}

void QR(float** A, float** R, float** Q, float** Qn, int n, float* b, float* bn, float* x)
{
    float eps = 0.0001;
    ToOne(Q, n);
    ToOne(Qn, n);
    float cij;
    float sij;
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

    float s;

    //Обратный проход
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += R[i][j] * x[j];
        }
        x[i] = (bn[i] - s) / R[i][i];
    }
}

void Inverse(float** Ainv, float** A, float** R, float** Q, float** Qn, int n, float* e, float* bn)
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
    int n;
    //std::cout << "Введите размер системы: ";
    fin >> n;
    //std::cout << "Введите матрицу системы:\n";
    //Переделать матрицы под локальное выделение памяти?
    float** A;
    A = new float* [n];

    float** Ainv;
    Ainv = new float* [n];

    float** Qn;
    Qn = new float* [n];

    float** Q;
    Q = new float* [n];

    float** R;
    R = new float* [n];

    float* b;
    b = new float[n];

    float* bn;
    bn = new float[n];

    float* x;
    x = new float[n];

    float* temp;
    temp = new float[n];

    float* DopCount;
    DopCount = new float[n];

    for (int i = 0; i < n; i++) {
        A[i] = new float[n];
        Ainv[i] = new float[n];
        Q[i] = new float[n];
        R[i] = new float[n];
        Qn[i] = new float[n];
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
    QR(A, R, Q, Qn, n, b, bn, x);
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }

    //3
    float* b1;
    b1 = new float[n];

    float bdop;

    for (int i = 0; i < n; i++)
    {
        b1[i] = b[i];
        for (int j = 0; j < n; j++)
        {
            b1[i] -= A[i][j] * x[j];
        }
    }
    std::cout << std::endl << "Норма невязки: " << (NormaVectora1(b1, n));
    float* x1;
    x1 = new float[n];

    x1[0] = 1;
    x1[1] = 1000;
    x1[2] = -20;
    x1[3] = 3;

    for (int i = 0; i < n; i++)
    {
        x1[i] = abs(abs(x1[i]) - abs(x[i]));
    }
    std::cout << std::endl << "Норма рассхождения численных решений от теоретических значений: " << (NormaVectora1(x1, n));
    /////////////////////////////////////////////////////////

    //4,6
    float* e;
    e = new float[n];


    Inverse(Ainv, A, R, Q, Qn, n, e, bn);
    std::cout << std::endl << "Обратная матрица:" << std::endl;
    Output(Ainv, x, n);

    std::cout << "Число обусловленности: " << NormaMatInf(A, n, DopCount) * NormaMatInf(Ainv, n, DopCount) << std::endl;
    delete[] e;
    //+оценка снизу с помощью полученных решений

    //5
    float vozm = 0.01;
    for (int i = 0; i < n; i++) {
        b[i] = b[i] + vozm;
    }

    QR(A, R, Q, Qn, n, b, bn, x1);
    std::cout << std::endl << "Решение системы с возмущенной правой частью:" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << x1[i] << '\n';
    }

    std::cout << std::endl;


    for (int i = 0; i < n; i++) {
        b[i] = b[i] - vozm;
    }
    for (int i = 0; i < n; i++) {
        b1[i] = vozm;
    }
    for (int i = 0; i < n; i++)
    {
        x1[i] = abs(x[i]) - abs(x1[i]);
        std::cout << x1[i] << std::endl;
    }

    std::cout << "Норма вектора расхождений решений: " << NormaVectora1(x1, n) << std::endl;
    std::cout << "Оценка снизу числа обусловленности: " << NormaVectoraInf(x1, n) / NormaVectoraInf(x, n) * NormaVectoraInf(b, n) / NormaVectoraInf(b1, n) << std::endl;

    delete[] b1;
    delete[] x1;
    ////////////////////////////////////////////////////////////////////////////

    std::cout << "Единичная:" << std::endl;
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