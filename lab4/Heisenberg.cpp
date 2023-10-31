#include <iostream>
#include <fstream>
#include <math.h>

void swap_ij(double** matrix, int i, int j, int n)
{
    double temp;
    for (int k = 0; k < n; k++)
    {
        temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
}
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
    //std::cout << bMatrix[0][0] << std::endl;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            Result[row][col] = 0;
            for (int inner = 0; inner < n; inner++)
            {
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

void sum_i2j_with_c(double** matrix, int i, int j, double c, int n) {
    for (int k = 0; k < n; k++) {
        matrix[j][k] = matrix[j][k] + c * matrix[i][k];
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

void QR(double** A, double** R, double** Q, double** Qn, int n)
{
    double eps = 0.0001;
    ToOne(Q, n);
    ToOne(Qn, n);
    double cij;
    double sij;
    double** Amem = new double* [n];
    for (int i = 0; i < n; i++)
    {
        Amem[i] = new double[n];
    }
    Copy(Amem, A, n);
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
    Copy(A, Amem, n);
    for (int i = 0; i < n; i++)
    {
        delete[] Amem[i];
    }
    delete[] Amem;
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

void Tmat(double** T, int k, int l, double alpha, double beta, int n)
{
    ToOne(T, n);
    T[k][k] = alpha;
    T[k][l] = beta;
    T[l][k] = -beta;
    T[l][l] = alpha;
}

void QRAlgorithm(double** A, int n)
{
    double** R = new double* [n];
    double** Q = new double* [n];
    double** Qn = new double* [n];
    for (int i = 0; i < n; i++)
    {
        R[i] = new double[n];
        Q[i] = new double[n];
        Qn[i] = new double[n];
    }

    double eps = 1e-3;
    //while (abs(A[n-1][n-1])>eps) {
    for (int i = 0; i < 1; i++)
    {
        QR(A, R, Q, Qn, n);
        MultiplyMatrix(A, Q, Qn, n);
        transpose(Q, n);
        MultiplyMatrix(Q, Qn, A, n);
        //Output(A, n);

        //Output(Qn, n);
        //Output(A, n);
    }
    //}

    for (int i = 0; i < n; i++)
    {
        delete[] R[i];
        delete[] Q[i];
        delete[] Qn[i];
    }
    delete[] R;
    delete[] Q;
    delete[] Qn;
}

void Heisenberg(double** A, int n)
{
    double alpha;
    double beta;
    double c;
    double** T;
    double** Result;
    T = new double* [n];
    Result = new double* [n];
    // int iter = 0;
    for (int i = 0; i < n; i++)
    {
        T[i] = new double[n];
        Result[i] = new double[n];
    }

    /*for (int k = 1; k < n - 1; k++) {
        //Обнуляем столбик
        //Находим максимальный элемент в столбце
        int max_ind = k;
        for (int x = k + 1; x < n; x++) {
            if (fabs(A[x][k-1]) > fabs(A[max_ind][k-1])) {
                max_ind = x;
            }
        }
        if (A[max_ind][k-1] == 0)
        {
            std::cout << "Матрица несовместная";
            return;
        }
        //Ставим его на k место
        swap_ij(A, k, max_ind, n);
        for (int x = k + 1; x < n; x++) {
            sum_i2j_with_c(A, k, x, -A[x][k-1] / A[k][k-1], n);
            //std::cout << k << " " << x+1 << std::endl;
        }
    }
    Output(A, n)*////////////////////////////////////////////////////////////////////Работает но не то
    for (int k = 1; k < n - 1; k++)
    {
        for (int l = k + 1; l < n; l++)
        {
            //iter++;
            alpha = A[k][k - 1] / sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1]);
            beta = A[l][k - 1] / sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1]);
            Tmat(T, k, l, alpha, beta, n);
            MultiplyMatrix(T, A, Result, n);
            transpose(T, n);
            MultiplyMatrix(Result, T, A, n);
            //c = alpha * alpha + beta * beta;
            //for (int i = 0; i < n; i++)
            //{
            //    A[k][i] = alpha * A[k][i] + beta * A[l][i];
            //    A[l][i] = -beta * A[k][i] + alpha * A[l][i];
            //}
        }
    }
    Output(A, n);
    //std::cout << iter << " ";
    for (int i = 0; i < n; i++)
    {
        delete[] T[i];
        delete[] Result[i];
    }
    delete[] T;
    delete[] Result;
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

    int n;
    fin >> n;

    double** A;
    A = new double* [n];

    double* b;
    b = new double[n];

    for (int i = 0; i < n; i++)
    {
        A[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> A[i][j];
        }
    }
    //Output(A, n);
    for (int i = 0; i < n; i++) {
        fin >> b[i];
    }

    Heisenberg(A, n);
    QRAlgorithm(A, n);
    Output(A, n);

    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }

    delete[] b;
    delete[] A;

    return 0;
}
