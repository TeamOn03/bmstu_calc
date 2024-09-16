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
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += R[i][j] * x[j];
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

void ReverseIter(double** AlE, double* x, int n)
{
    double* y = new double[n];
    ToOne(x, n);

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
        x[i] = estime[i] / NormaVectora2(y, n);
    }
    ToOne(estimeMem, n);

    Inverse(AlEmem, AlE, n);
    while (NormaVectora2(estime, n) > eps)
    {
        MultiplyMatrixToVector(AlEmem, x, y, n);
        for (int i = 0; i < n; i++)
        {
            x[i] = y[i] / NormaVectora2(y, n);
        }
        MultiplyMatrixToVector(AlE, x, estime, n);
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
        x = new double[n];
        y = new double[n];
    }
    void even(double a, double b)
    {
        for (int i = 0; i < n; i++)
        {
            x = a + (b - a) / i;
        }
    }
    double Ck(int k, double x)
    {
        double result = 0;
        for (int i = 0; i < n; i++)
        {
            result *= (x - x[i]) / (x[k] - x[j]);
        }
        return result;
    }
    ~Lagranj()
    {
        delete[] x;
        delete[] y;
    }
};

int main() {
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
    for (int i = 0; i < n-1; i++)
    {
        lymbda = evalues[i];


        MatrixOnCoef(E, -lymbda, n);
        SumMat(A, E, A, n);
        ReverseIter(A, x, n);

        std::cout << std::endl;
        OutputVect(x, n);

        MatrixOnCoef(E, -1, n);
        SumMat(A, E, A, n);
        MatrixOnCoef(E, 1 / lymbda, n);
        //lymbda = Reley(A, x, n);
        std::cout << lymbda << std::endl;
        //(x, n);
        evalues[i] = lymbda;
    }
    x[0]=12e-70;
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
}
