#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;


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

void ToOne(double** OldM, int n)
{
    ToNull(OldM, n);
    for (int j = 0; j < n; j++)
    {
        OldM[j][j] = 1;
    }
}

void Output(double** A, double* b, int n) {
    //cout << "Матрица:\n" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j];
            cout << ' ';
        }
        cout << "   ";
        cout << b[i];
        cout << endl;
    }
}


void Norm(int j, int n, double** Mat, double** QMat)
{
    double maks = 0;
    int k = 0;
    for (int i = j; i < n; i++)
    {
        if (maks < max(maks, abs(Mat[i][j])))
        {
            k = i;
            maks = max(maks, abs(Mat[i][j]));
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

void MultiplyMatrixToVector(double** Matrix, double* vector, double* result, int n){
    for(int i=0;i<n;i++){
        double s=0;
        for(int j=0;j<n;j++){
            s+=Matrix[i][j]*vector[j];
        }
        result[i]=s;
    }
}


int main() {
    ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");
    double eps = 0.0001;
    int n;
    double s;
    cout << "Введите размер системы: ";
    fin >> n;
    cout << "Введите матрицу системы:\n";
    //Переделать матрицы под локальное выделение памяти?
    double** A;
    A = new double* [n];

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

    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        Q[i] = new double[n];
        R[i] = new double[n];
        Qn[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> A[i][j];
        }
    }

    cout << "Введите свободные коэффициенты:\n";
    for (int i = 0; i < n; i++) {
        fin >> b[i];
    }

    ToOne(Q, n);
    ToOne(Qn, n);
    double cij;
    double sij;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (abs(A[j][j]) < eps)
                Norm(j, n, A, Q);
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

    transpose(Q, n);
    cout<<"\nМатрица Q: \n";
    Output(Q, b, n);
    cout<<"\nМатрица R: \n";
    Output(R, bn, n);

    MultiplyMatrix(Q, R, A, n);

    transpose(Q,n);
    MultiplyMatrixToVector(Q,b,bn,n);

    //Обратный проход
    for(int i=n-1;i>=0;i--){
        s=0;
        for(int j=i+1;j<n;j++){
            s+=R[i][j]*x[j];
        }
        x[i]=(bn[i]-s)/R[i][i];
    }


    cout<<"\nОтвет (QR):\n";
    for(int i=0;i<n;i++){
        cout<<x[i]<<'\n';
    }


    for (int i = 0; i < n; i++) {
        delete[] A[i];
        delete[] Q[i];
        delete[] R[i];
        delete[] Qn[i];
    }
    delete[] Q;
    delete[] Qn;
    delete[] R;
    delete[] A;
    delete[] b;
    delete[] bn;
    delete[] x;
    delete[] temp;
    fin.close();
    return 0;
}