#include <iostream>
#include <fstream>
#include <cmath>

void sum_i2j_with_c(float** matrix, int i, int j, float c, int n) {
    for (int k = 0; k < n; k++) {
        matrix[j][k] = matrix[j][k] + c * matrix[i][k];
    }
}

void swap_ij(float** matrix, int i, int j, int n) {
    float temp;
    for (int k = 0; k < n; k++) {
        temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
}

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

void Copy(float* data, float* newData, int n)
{
    for (int i = 0; i < n; i++)
    {
        data[i] = newData[i];
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
    float maks = 0;
    int k = 0;
    for (int i = j; i < n; i++)
    {
        if (maks < abs(Mat[i][j]))
        {
            k = i;
            maks = abs(Mat[i][j]);
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

void Gauss(float** A, float* b, float* x, int n, float** Acopy, float* bcopy)
{
    Copy(bcopy, b, n);
    Copy(Acopy, A, n);
    //Прямой проход
    for (int k = 0; k < n; k++) {
        //Находим максимальный элемент в столбце
        int max_ind = k; #include <iostream>
#include <fstream>
#include <cmath>


            void Output(double* A, int n)
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

        double DopSum(double** A, double* x, int low_board, int up_board, int i)
        {
            double result = 0;
            for (int j = low_board; j < up_board; j++)
            {
                result += A[i][j] / A[i][i] * x[j];
            }
            return result;
        }

        void Zeidel(double** A, double** C, double* b, double* y, double* x, double* x_old, int n, double* diff)
        {
            double eps = 1e-10;
            //Начальное приближение x^0 любое
            ToNull(x, n);
            double omega = 1;

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

            Copy(x_old, x, n);
            x[0] = (1 - omega) * x_old[0] - omega * DopSum(A, x_old, 1, n, 0) + omega * b[0] / A[0][0];

            x[1] = -omega * A[1][0] / A[1][1] * x[0] + (1 - omega) * x_old[1] - omega * DopSum(A, x_old, 2, n, 1) + omega * b[1] / A[1][1];

            for (int i = 2; i < n - 1; i++)
            {
                x[i] = -omega * DopSum(A, x, 0, i, i) + (1 - omega) * x_old[i] - omega * DopSum(A, x_old, i + 1, n, i) + omega * b[i] / A[i][i];
            }

            x[n - 1] = -omega * DopSum(A, x, 0, n - 1, n - 1) + (1 - omega) * x_old[n - 1] + omega * b[n - 1] / A[n - 1][n - 1];

            Copy(diff, x, n);
            DiffVect(diff, x_old, n);
            while (NormaVectora1(diff, n) > eps)
            {
                Copy(x_old, x, n);
                x[0] = (1 - omega) * x_old[0] - omega * DopSum(A, x_old, 1, n, 0) + omega * b[0] / A[0][0];

                x[1] = -omega * A[1][0] / A[1][1] * x[0] + (1 - omega) * x_old[1] - omega * DopSum(A, x_old, 2, n, 1) + omega * b[1] / A[1][1];

                for (int i = 2; i < n - 1; i++)
                {
                    x[i] = -omega * DopSum(A, x, 0, i, i) + (1 - omega) * x_old[i] - omega * DopSum(A, x_old, i + 1, n, i) + omega * b[i] / A[i][i];
                }

                x[n - 1] = -omega * DopSum(A, x, 0, n - 1, n - 1) + (1 - omega) * x_old[n - 1] + omega * b[n - 1] / A[n - 1][n - 1];

                Copy(diff, x, n);
                DiffVect(diff, x_old, n);
                Output(x, n);
            }

        }

        int main() {
            std::ifstream fin("testZeidel.txt");
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

            for (int i = 0; i < n; i++) {
                A[i] = new double[n];
                C[i] = new double[n];
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
            Zeidel(A, C, b, y, x, x_old, n, diff);
            Output(x, n);

            for (int i = 0; i < n; i++) {
                delete[] A[i];
                delete[] C[i];
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
        float temp = b[k];
        b[k] = b[max_ind];
        b[max_ind] = temp;

        //Обнуляем столбик
        for (int x = k + 1; x < n; x++) {
            b[x] = b[x] + b[k] * (-A[x][k] / A[k][k]);
            sum_i2j_with_c(A, k, x, -A[x][k] / A[k][k], n);
        }
    }

    float s;

    //Output(A, b, n);
    //Обратный проход
    for (int i = n - 1; i >= 0; i--) {
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }

    Copy(A, Acopy, n);
    Copy(b, bcopy, n);
    //std::cout << "\nОтвет (Прямой метод Гаусса):\n";
    //for (int i = 0; i < n; i++) {
    //    std::cout << x[i] << '\n';
    //}
}

void Inverse(float** Ainv, float** A, int n, float* e, float** Acopy, float* bcopy)
{
    for (int i = 0; i < n; i++)
    {
        ToNull(e, n);
        e[i] = 1;
        Gauss(A, e, Ainv[i], n, Acopy, bcopy);
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
    float** A;
    A = new float* [n];
    float** Acopy;
    Acopy = new float* [n];
    float** E;
    E = new float* [n];
    float** Ainv;
    Ainv = new float* [n];
    float* b;
    b = new float[n];
    float* bcopy;
    bcopy = new float[n];
    float* DopCount;
    DopCount = new float[n];
    float* x;
    x = new float[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new float[n];
        Ainv[i] = new float[n];
        E[i] = new float[n];
        Acopy[i] = new float[n];
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

    //Copy(E, A, n);
    Gauss(A, b, x, n, Acopy, bcopy);

    for (int i = 0; i < n; i++) {
        std::cout << x[i] << '\n';
    }

    //3
    float* b1;
    b1 = new float[n];

    float bdop;

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

    //3
    float* xReal;
    xReal = new float[n];
    xReal[0] = 1;
    xReal[1] = 1000;
    xReal[2] = -20;
    xReal[3] = 3;

    for (int i = 0; i < n; i++)
    {
        xReal[i] = abs(abs(xReal[i]) - abs(x[i]));
    }
    std::cout << std::endl << "Норма вектора расхождений от теоретических решений: " << (NormaVectora1(xReal, n));
    /////////////////////////////////////////////////////////

    //4,6
    float* e;
    e = new float[n];

    //Copy(A, E, n);
    Inverse(Ainv, A, n, e, Acopy, bcopy);
    std::cout << std::endl << "Обратная матрица:" << std::endl;
    Output(Ainv, x, n);

    std::cout << "Число обусловленности: " << NormaMatInf(A, n, DopCount) * NormaMatInf(Ainv, n, DopCount) << std::endl;
    delete[] e;
    /////////////////////////////////////////////////////////

    //5
    float vozm = 0.01;
    for (int i = 0; i < n; i++) {
        b[i] = b[i] + vozm;
    }

    float* x1;
    x1 = new float[n];

    //Copy(A, E, n);
    Gauss(A, b, x1, n, Acopy, bcopy);
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
        delete[] Acopy[i];
    }
    delete[] Acopy;
    delete[] bcopy;
    delete[] DopCount;
    delete[] E;
    delete[] Ainv;
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}