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

double NormaVectora2(double* b, int n)
{
    double answer = 0;
    for (int i = 0; i < n; i++)
    {
        answer += (b[i]) * (b[i]);
    }
    return sqrt(answer);
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
    double** Amem = new double* [n];
    for (int i = 0; i < n; i++)
    {
        Amem[i] = new double[n];
    }
    //std::cout << bMatrix[0][0] << std::endl;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            Amem[row][col] = 0;
            for (int inner = 0; inner < n; inner++)
            {
                Amem[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
            }
        }
        //std::cout << "\n";
    }
    Copy(Result, Amem, n);
    for (int i = 0; i < n; i++)
    {
        delete[] Amem[i];
    }
    delete[] Amem;
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
    //Output(A, n);
    double s;

    //Обратный проход
    //OutputVect(x, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (std::isnan(R[i][j]))
            {
                R[i][j] = 1e-30;
            }
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        if (std::isnan(bn[i]))
        {
            bn[i] = 1e-50;
        }
        for (int j = i + 1; j < n; j++) {
            if (std::isnan(R[i][j]))
            {
                R[i][j] = 1e-30;
            }
            s += R[i][j] * x[j];
        }
        if (abs(R[i][i]) < 1e-60 or std::isnan(R[i][i]))
        {
            if (R[i][i] < 0)
                R[i][i] = -1e-30;
            else
                R[i][i] = 1e-30;
        }
        x[i] = (bn[i] - s) / R[i][i];
        //std::cout << x[i];
    }
    //OutputVect(x, n);
}

void QR(double** A, double** R, double** Q, int n)
{
    double eps = 0.0001;
    ToOne(Q, n);
    double cij;
    double sij;
    double** Amem = new double* [n];
    double** Qn = new double* [n];
    for (int i = 0; i < n; i++)
    {
        Amem[i] = new double[n];
        Qn[i] = new double[n];
    }
    ToOne(Qn, n);
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
    transpose(Q, n);
    Copy(A, Amem, n);
    for (int i = 0; i < n; i++)
    {
        delete[] Amem[i];
        delete[] Qn[i];
    }
    delete[] Amem;
    delete[] Qn;
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
    //QR(A, R, Q, Qn, n, e, bn, Ainv[0]);

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


void Tmat(double** T, int k, int l, double alpha, double beta, int n)
{
    ToOne(T, n);
    T[k][k] = alpha;
    T[k][l] = beta;
    T[l][k] = -beta;
    T[l][l] = alpha;
}

void getL(double** A, double** L, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            L[j][i] = A[j][i];
        }
    }
}

double shift(double** A, int n)
{
    double** AlE = new double* [n];
    double** E = new double* [n];
    for (int i = 0; i < n; i++)
    {
        E[i] = new double[n];
        AlE[i] = new double[n];
    }
    ToOne(E, n);
    double mem = A[n - 1][n - 1];

    MatrixOnCoef(E, -mem, n);
    SumMat(A, E, AlE, n);

    Copy(A, AlE, n);

    for (int i = 0; i < n; i++)
    {
        delete[] E[i];
        delete[] AlE[i];
    }
    delete[] E;
    delete[] AlE;

    return mem;
}

void unshift(double** A, double mem, int n)
{
    double** AlE = new double* [n];
    double** E = new double* [n];
    for (int i = 0; i < n; i++)
    {
        E[i] = new double[n];
        AlE[i] = new double[n];
    }
    ToOne(E, n);
    MatrixOnCoef(E, mem, n);
    SumMat(A, E, AlE, n);

    Copy(A, AlE, n);

    for (int i = 0; i < n; i++)
    {
        delete[] E[i];
        delete[] AlE[i];
    }
    delete[] E;
    delete[] AlE;
}

void QRAlgorithm(double** A, int n)
{
    double* nuls = new double[n - 1];
    double** R = new double* [n];
    double** Q = new double* [n];
    double** L = new double* [n];
    for (int i = 0; i < n; i++)
    {
        R[i] = new double[n];
        Q[i] = new double[n];
        L[i] = new double[n];
    }

    int iter = 0;
    double mem = 0;
    bool flag = 0;

    ToNull(R, n);
    ToNull(Q, n);

    double eps = 1e-6;
    //double eps1 = 1e-2;
    ToNull(L, n);
    getL(A, L, n);
    //double nmem = n;
    while (NormaMat1(L, n) > eps)
        //while (n > 1)
    {
        iter++;
        if (flag == 0)
        {
            mem += shift(A, n);
            flag = 1;
            //Output(A, n);
            //std::cout << mem << std::endl;
        }
        QR(A, R, Q, n);
        MultiplyMatrix(R, Q, A, n);
        getL(A, L, n);
        for (int i = 0; i < n - 1; i++)
        {
            nuls[i] = A[n - 1][i];
        }
        //std::cout << NormaVectora1(nuls, n - 1) << std::endl;
        if (NormaVectora1(nuls, n - 1) < eps)
        {
            if (flag == 1)
            {
                flag = 0;
                unshift(A, mem, n);
                mem = 0;
            }
            //Output(A, n);
            n = n - 1;

        }
    }
    std::cout << iter << std::endl;
    for (int i = 0; i < n; i++)
    {
        delete[] L[i];
        delete[] R[i];
        delete[] Q[i];
    }
    delete[] L;
    delete[] R;
    delete[] Q;
    delete[] nuls;
}

void Heisenberg(double** A, int n)
{
    double alpha;
    double beta;
    //double c;
    //double** T;
    double** Result;
    //T = new double* [n];
    Result = new double* [n];
    for (int i = 0; i < n; i++)
    {
        //T[i] = new double[n];
        Result[i] = new double[n];
    }

    double Aki;
    double Ali;
    double Ail;
    double Aik;
    for (int k = 1; k < n - 1; k++)
    {
        for (int l = k + 1; l < n; l++)
        {
            alpha = A[k][k - 1] / sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1]);
            beta = A[l][k - 1] / sqrt(A[k][k - 1] * A[k][k - 1] + A[l][k - 1] * A[l][k - 1]);
            /*Tmat(T, k, l, alpha, beta, n);
            MultiplyMatrix(T, A, Result, n);
            transpose(T, n);
            MultiplyMatrix(Result, T, A, n);*/
            for (int i = 0; i < n; i++)
            {
                Aki = A[k][i];
                Ali = A[l][i];
                A[k][i] = alpha * Aki + beta * Ali;
                A[l][i] = alpha * Ali - beta * Aki;
            }
            for (int i = 0; i < n; i++)
            {
                Aik = A[i][k];
                Ail = A[i][l];
                A[i][k] = alpha * Aik + beta * Ail;
                A[i][l] = alpha * Ail - beta * Aik;
            }
            //Output(A, n);
        }
    }
    Output(A, n);
    for (int i = 0; i < n; i++)
    {
        //delete[] T[i];
        delete[] Result[i];
    }
    //delete[] T;
    delete[] Result;
}
void SimpleIter(double** A, double* b, double* x, int n)
{
    double* x_old = new double[n];
    double* diff = new double[n];
    int iterations = 0;
    double eps = 1e-2;
    //Начальное приближение x^0 любое
    ToNull(x, n);
    double tao = 1e-4;
    matrix_diag_poz_sign(A, n, b);

    //std::cout << "Норма матрицы С: " << NormaMat1(C, n) << std::endl;
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
    //std::cout << "Количество иттераций: " << iterations << "\n";
    delete[] diff;
    delete[] x_old;
}

void Zeidel(double** abcd, int n, double* x)
{
    int iterations = 0;
    double omega = 1;
    double eps = 1e-10;

    double* diff = new double[n];
    double* x_old = new double[n];

    for (int i = 0; i < n; i++) {
        x_old[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        x[i] = 1;
    }


    Copy(diff, x, n);
    DiffVect(diff, x_old, n);
    double normD = NormaVectora1(diff, n);
    while (normD > eps)
    {
        Copy(x_old, x, n);
        for (int i = 0; i < n; i++) {
            double temp = 0;
            if (i == 0)
                temp += abcd[2][i] * x[i + 1];
            if (i == n - 1)
                temp += abcd[0][i - 1] * x[i - 1];
            if ((i != 0) && (i != n - 1))
                temp += abcd[0][i - 1] * x[i - 1] + abcd[2][i] * x[i + 1];
            x[i] = omega * (abcd[3][i] - temp) / abcd[1][i] + (1 - omega) * x[i];
        }
        Copy(diff, x, n);
        DiffVect(diff, x_old, n);
        normD = NormaVectora1(diff, n);
        iterations++;
    }
    std::cout << "Количество итераций: " << iterations << "\n";

    delete[] x_old;
    delete[] diff;
}

void ReverseIter(double** AlE, double* x, int n)
{
    double* y = new double[n];
    //ToNull(x, n);


    double** AlEmem = new double* [n];
    for (int i = 0; i < n; i++)
    {
        AlEmem[i] = new double[n];
    }
    double eps = 1e-1;
    double* estime = new double[n];
    double* estimeMem = new double[n];
    MultiplyMatrixToVector(AlE, x, estime, n);
    for (int i = 0; i < n; i++)
    {
        x[i] = estime[i] / NormaVectora1(estime, n);
    }
    //std::cout << NormaVectora2(x, n);
    ToOne(estimeMem, n);

    //Inverse(AlEmem, AlE, n);
    while (NormaVectora1(estime, n) > eps)
    {
        //MultiplyMatrixToVector(AlEmem, x, y, n);
        SimpleIter(AlE, x, y, n);
        for (int i = 0; i < n; i++)
        {
            x[i] = y[i] / NormaVectora1(y, n);
        }
        MultiplyMatrixToVector(AlE, x, estime, n);
        //std::cout << NormaVectora1(estime, n) << std::endl;
    }
    //std::cout << NormaVectoraInf(x, n) << " " << NormaVectora1(x, n) << std::endl;
    for (int i = 0; i < n; i++)
    {
        delete[] AlEmem[i];
    }
    delete[] AlEmem;
    delete[] estime;
    delete[] estimeMem;
    delete[] y;
}

double MultiplyVectors(double* a, double* b, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

double Reley(double** A, double* x, int n)
{
    double lymbda = 0;
    double* xextra = new double[n];
    double* estime = new double[n];

    double** AlE = new double* [n];
    double** AlEinv = new double* [n];
    double** E = new double* [n];
    for (int i = 0; i < n; i++)
    {
        E[i] = new double[n];
        AlE[i] = new double[n];
        AlEinv[i] = new double[n];
    }
    double eps = 1e-3;
    ToOne(estime, n);

    /* [0] = -0.8644;
    x[1] = -0.0033;
    x[2] = -0.2445;
    x[3] = -0.4391;

    std::cout << NormaVectora1(x, n) << std::endl;*/
    //ToNull(AlEinv, 0);

    while (NormaVectora2(estime, n) > eps)
    {
        //MatrixOnCoef(A, 1.55, n);
        MultiplyMatrixToVector(A, x, xextra, n);
        lymbda = MultiplyVectors(xextra, x, n);
        //MatrixOnCoef(A, 1 / 1.55, n);

        ToOne(E, n);
        MatrixOnCoef(E, -lymbda, n);
        SumMat(A, E, AlE, n);
        Inverse(AlEinv, AlE, n);
        //Output(AlEinv, n);
        MultiplyMatrixToVector(AlEinv, x, xextra, n);
        for (int i = 0; i < n; i++)
        {
            x[i] = xextra[i] / NormaVectora2(xextra, n);
        }
        MultiplyMatrixToVector(AlE, x, estime, n);

        //std::cout << NormaVectora1(estime, n) << std::endl;
    }

    for (int i = 0; i < n; i++)
    {
        delete[] E[i];
        delete[] AlE[i];
        delete[] AlEinv[i];
    }
    delete[] estime;
    delete[] E;
    delete[] AlE;
    delete[] AlEinv;
    delete[] xextra;
    return lymbda;
}

/*В отдельный класс записываем, сохраняем поля abcd потом на каждом иксе подсталвяем коэффициенты просто
double** abcd = new double* [4];
        for (int i = 0; i < 4; i++)
        {
            abcd[i] = new double[n + 1];
        }

        double* bv = new double[n + 1];
        double** Syst = new double* [n + 1];
        for (int i = 0; i < n + 1; i++)
        {
            Syst[i] = new double[n + 1];
        }

        for (int i = 0; i < n; i++)
        {
            abcd[0][i] = y[i];
        }

        ToNull(Syst, n);
        Syst[0][0] = 1;
        bv[0] = 0;
        for (int i = 1; i < n; i++)
        {
            Syst[i][i - 1] = this->x[i] - this->x[i - 1];
            Syst[i][i] = 2 * (this->x[i + 1] - this->x[i - 1]);
            Syst[i][i + 1] = this->x[i + 1] - this->x[i];
            bv[i] = 3 * ((y[i + 1] - y[i]) / (this->x[i + 1] - this->x[i]) + (y[i] - y[i - 1]) / (this->x[i] - this->x[i - 1]));
        }
        Syst[int(n)][int(n)] = 1;
        bv[int(n)] = 0;

        SimpleIter(Syst, bv, abcd[2], n + 1);*/

class Lagranj
{
private:
    double* x;
    double* y;
    double n;
public:
    Lagranj(int n)
    {
        this->n = n;
        x = new double[n + 1];
        y = new double[n + 1];
    }
    double func(double x)
    {
        return 5;
        //return x * x;
        //return 1/(1+25*x*x);
        //return 1/atan(1+10*x*x);
        //return pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
    }
    void even(double a, double b)
    {
        for (int i = 0; i <= n; i++)
        {
            x[i] = a + (b - a) / n * i;
            y[i] = func(x[i]);
        }
        //OutputVect(x, n+1);
    }
    void cheb(double a, double b)
    {
        for (int i = 0; i <= n; i++)
        {
            x[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * acos(-1) / (2 * (n + 1)));
            y[i] = func(x[i]);
        }
        //OutputVect(x, n+1);
    }
    double Ck(int k, double x)
    {
        double result = 1;
        for (int i = 0; i <= n; i++)
        {
            if (k != i)
                result *= (x - this->x[i]) / (this->x[k] - this->x[i]);
        }
        return result;
    }
    /*double Si(int k, double x)
    {
        return ;
    }*/
    double LInterpol(double x)
    {
        double result = 0;
        for (int k = 0; k <= n; k++)
        {
            result += Ck(k, x) * y[k];
        }
        return result;
    }
    /*double KubSInterpol(double x)
    {


        double result = 0;
        for (int k = 0; k < n; k++)
        {
            result += a + b * (x - this->x[k - 1]) + c * pow((x - this->x[k - 1]), 2) + d * pow((x - this->x[k - 1]), 3);
        }

        for (int i = 0; i < n + 1; i++)
        {
            delete[] Syst[i];
            delete[] abcd;
        }
        delete[] abcd;
        delete[] Syst;
        delete[] bv;

        return result;
    }*/
    ~Lagranj()
    {
        delete[] x;
        delete[] y;
    }
};


class Spline
{
private:
    bool flagCheb = false;
    double* x;
    double* y;
    int n;
    double** abcd = new double* [4];
    double** Syst = new double* [4];
public:
    Spline(int n)
    {
        this->n = n;
        x = new double[n + 1];
        y = new double[n + 1];
        for (int i = 0; i < 4; i++)
        {
            abcd[i] = new double[n + 2];
        }
    }
    double func(double x)
    {
        //return 5;
        //return x;
        //return x * x;
        //return exp(x);
        //return 1/(1+25*x*x);
        //return 1/atan(1+10*x*x);
        return pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
    }
    void even(double a, double b)
    {
        for (int i = 0; i <= n; i++)
        {
            x[i] = a + (b - a) / n * i;
            y[i] = func(x[i]);
        }
        //OutputVect(x, n+1);
    }
    void cheb(double a, double b)
    {
        for (int i = 0; i <= n; i++)
        {
            x[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * acos(-1) / (2 * (n + 1)));
            y[i] = func(x[i]);
        }
        flagCheb = true;
        //OutputVect(x, n + 1);
    }
    void FindAbcd()
    {
        for (int i = 0; i < 4; i++)
        {
            Syst[i] = new double[n + 1];
            ToNull(abcd[i], n + 2);
            ToNull(Syst[i], n + 1);
        }

        for (int i = 1; i < n + 1; i++)
        {
            abcd[0][i] = this->y[i - 1];
        }

        Syst[1][0] = 1;
        Syst[3][0] = 0;
        Syst[1][n] = 1;
        Syst[3][n] = 0;
        for (int i = 1; i < n; i++)
        {
            Syst[0][i - 1] = this->x[i] - this->x[i - 1];
            Syst[1][i] = 2 * (this->x[i + 1] - this->x[i - 1]);
            Syst[2][i - 1] = this->x[i + 1] - this->x[i];
            Syst[3][i] = 3 * ((y[i + 1] - y[i]) / (this->x[i + 1] - this->x[i]) - (y[i] - y[i - 1]) / (this->x[i] - this->x[i - 1]));
        }
        Zeidel(Syst, n + 1, abcd[2]);
        //OutputVect(abcd[2], n+1);

        for (int i = 1; i < n + 1; i++)
        {
            abcd[1][i] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]) - (abcd[2][i + 1] + 2 * abcd[2][i]) * (x[i] - x[i - 1]) / 3;
            abcd[3][i] = (abcd[2][i + 1] - abcd[2][i]) / (3 * (x[i] - x[i - 1]));
        }

        /*OutputVect(abcd[0], n + 2);
        OutputVect(abcd[1], n + 2);
        OutputVect(abcd[2], n + 2);
        OutputVect(abcd[3], n + 2);*/


        //delete[] Syst;
    }
    double KubSInterpol(double x)
    {
        double result = 0;
        int k = 1;
        if (flagCheb)
            for (int i = 1; i <= n + 1; i++)
            {
                if ((x < this->x[i - 1]) and (x >= this->x[i]))
                {
                    k = i;
                    break;
                }
            }
        else
            for (int i = 1; i <= n + 1; i++)
            {
                if ((x >= this->x[i - 1]) and (x < this->x[i]))
                {
                    k = i;
                    break;
                }
            }
        //for (int k = 1; k < n; k++)
        //{
        result += abcd[0][k] + abcd[1][k] * (x - this->x[k - 1]) + abcd[2][k] * pow((x - this->x[k - 1]), 2) + abcd[3][k] * pow((x - this->x[k - 1]), 3);
        //}

        return result;
    }
    ~Spline()
    {
        for (int i = 0; i < 4; i++)
        {
            delete[] abcd[i];
            delete[] Syst[i];
        }
        /*for (int i = 0; i < 4; i++)
        {
            delete[] Syst[i];
        }*/

        delete[] Syst;
        delete[] x;
        delete[] y;
        delete[] abcd;
    }
};

int main() {
    std::ofstream fout;
    fout.open("func.txt");
    setlocale(LC_ALL, "Russian");


    int n = 50;
    double a = -1;
    double b = 1;
    //fin >> n >> a >> b;

    double step = 0.01;

    //double* y = new double[(b-a)/step];

    Lagranj xsqr(n);
    //Spline xsqr(n);
    xsqr.even(a, b);
    //xsqr.FindAbcd();
    for (double x = a; x <= b; x += step)
    {
        //y[int((x-a)/step)] = xsqr.Interpol(x);
        fout << x << " " << xsqr.LInterpol(x) << std::endl;
        //fout << x << " " << xsqr.KubSInterpol(x) << std::endl;
    }

    //OutputVect(y, (b - a) / step);
    fout.close();
    //delete[] y;
    //return 0;
    system("Pause");
}



/*int main() {
    std::ifstream fin("test2.txt");
    setlocale(LC_ALL, "Russian");

    int n;
    fin >> n;

    double** A;
    A = new double* [n];
    double** Amem;
    Amem = new double* [n];
    double** E;
    E = new double* [n];

    double* x;
    x = new double[n];
    double* evalues;
    evalues = new double[n];

    for (int i = 0; i < n; i++)
    {
        Amem[i] = new double[n];
        A[i] = new double[n];
        E[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> A[i][j];
        }
    }
    //Inverse(E, A, n);///////////////////////////////////////////////
    Copy(Amem, A, n);

    //Heisenberg(A, n);
    Output(A, n);
    QRAlgorithm(A, n);
    Output(A, n);
    for (int i = 0; i < n; i++)
    {
        evalues[i] = A[i][i];
    }
    double lymbda;

    Copy(A, Amem, n);
    ToOne(E, n);
    for (int i = 0; i < n; i++)
    {
        Copy(A, Amem, n);
        lymbda = evalues[i];


        MatrixOnCoef(E, -lymbda, n);
        SumMat(A, E, A, n);
        ToOne(x, n);
        //x[i] = 1.11;
        ReverseIter(A, x, n);
        ToOne(E, n);



        std::cout << std::endl;
        OutputVect(x, n);

        Copy(A, Amem, n);
        lymbda = Reley(A, x, n);
        std::cout << lymbda << std::endl;
        /*if (i == 2)
        {
            x[0] = 4.82e-6;
            x[1] = 4.53e-6;
            x[2] = 0.52732;
            x[3] = 0.850651;
        }
        OutputVect(x, n);
        evalues[i] = lymbda;
    }
    /*x[0] = 12e-70;
    x[1] =23e-70;
    x[2] = 43e-70;
    x[3] = 1.000000312;
    std::cout << std::endl;
    lymbda = 1;
    OutputVect(x, n);
    std::cout << lymbda << std::endl;




    double norma = NormaVectora1(evalues, n);
    evalues[0] = 1;
    evalues[1] = 1;
    evalues[2] = 1;
    evalues[3] = 1;
    std::cout << std::endl;
    std::cout << abs(NormaVectora1(evalues, n) - norma) << std::endl;

    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
        delete[] Amem[i];
        delete[] E[i];
    }

    delete[] evalues;
    delete[] E;
    delete[] x;
    delete[] A;
    delete[] Amem;

    return 0;
}*/
