#include <iostream>
#include <fstream>
#include <cmath>

void sum_i2j_with_c(double** matrix, int i, int j, double c, int n) {
    for (int k = 0; k < n; k++) {
        matrix[j][k] = matrix[j][k] + c * matrix[i][k];
    }
}

void swap_ij(double** matrix, int i, int j, int n) {
    double temp;
    for (int k = 0; k < n; k++) {
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
        answer += abs(b[i]);
    }
    return answer;
}

double* DopCount;

double NormaMat1(double** A, int n, double* DopCount)
{
    double maks = 0;
    DopCount = new double[n];

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

double NormaMatInf(double** A, int n, double* DopCount)
{
    double maks = 0;

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

void MultiplyMatrixToVector(double** Matrix, double* vector, double* result, int n) {
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) {
            s += Matrix[i][j] * vector[j];
        }
        result[i] = s;
    }
}

void Gauss(double** A, double* b, double* x, int n)
{
    //Прямой проход
    for (int k = 0; k < n; k++) {
        //Находим максимальный элемент в столбце
        int max_ind = k;
        for (int x = k + 1; x < n; x++) {
            if (fabs(A[x][k]) > fabs(A[max_ind][k])) {
                max_ind = x;
            }
        }
        if (A[max_ind][k] == 0) {
            std::cout << "Матрица несовместная";
            return;
        }
        //Ставим его на k место
        swap_ij(A, k, max_ind, n);
        double temp = b[k];
        b[k] = b[max_ind];
        b[max_ind] = temp;

        //Обнуляем столбик
        for (int x = k + 1; x < n; x++) {
            b[x] = b[x] + b[k] * (-A[x][k] / A[k][k]);
            sum_i2j_with_c(A, k, x, -A[x][k] / A[k][k], n);
        }
    }

    double s;

    //Output(A, b, n);
    //Обратный проход
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }

    //std::cout << "\nОтвет (Прямой метод Гаусса):\n";
    //for (int i = 0; i < n; i++) {
    //    std::cout << x[i] << '\n';
    //}
}

void Inverse(double** Ainv, double** A, int n, double* e)
{
    for (int i = 0; i < n; i++)
    {
        ToNull(e, n);
        e[i] = 1;
        Gauss(A, e, Ainv[i], n);
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
    double** A;
    A = new double* [n];
    double** E;
    E = new double* [n];
    double** Ainv;
    Ainv = new double* [n];
    double* b;
    b = new double[n];
    double* DopCount;
    DopCount = new double[n];
    double* x;
    x = new double[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new double[n];
        Ainv[i] = new double[n];
        E[i] = new double[n];
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

    Gauss(A, b, x, n);

    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }

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
        b1[i] = b[i] - bdop;
    }
    std::cout << std::endl << "Норма невязки: " << (NormaVectora1(b1, n));
    /////////////////////////////////////////////////////////

    //4,6
    double* e;
    e = new double[n];


    Inverse(Ainv, A, n, e);
    std::cout << std::endl << "Обратная матрица:" << std::endl;
    Output(Ainv, x, n);

    std::cout << "Число обусловленности: " << NormaMatInf(A, n, DopCount) * NormaMatInf(Ainv, n, DopCount) << std::endl;
    delete[] e;
    //+оценка снизу с помощью полученных решений

    //5
    double vozm = 0.01;
    for (int i = 0; i < n; i++) {
        b[i] = b[i] + vozm;
    }

    double* x1;
    x1 = new double[n];

    Gauss(A, b, x1, n);
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
    MultiplyMatrix(A, Ainv, E, n);
    Output(E, b, n);


    fin.close();
    for (int i = 0; i < n; i++) {
        delete[] A[i];
        delete[] Ainv[i];
        delete[] E[i];
    }
    delete[] DopCount;
    delete[] E;
    delete[] Ainv;
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}