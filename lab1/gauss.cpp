#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

//ДЛЯ ТЕСТА
void sum_i2j_with_c(double** matrix, int i,int j, double c,int n){
    for(int k=0;k<n;k++){
        matrix[j][k]=matrix[j][k]+c*matrix[i][k];

    }
}

void Output(double** A,double* b,int n){
    cout << "Матрица:\n"<<endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                cout << A[i][j];
                cout<<' ';
            }
            cout<<"   ";
            cout<<b[i];
            cout << endl;
        }
}

void swap_ij(double** matrix,int n,int i,int j){
    double temp;
    for(int k=0;k<n;k++){
        temp=matrix[i][k];
        matrix[i][k]=matrix[j][k];
        matrix[j][k]=temp;
    }
}


int main() {
    ifstream fin("test.txt");
    setlocale(LC_ALL, "Russian");
    int n;
    double s;
    fin>>n;
    double** A;
    A=new double*[n];
    double* b;
    b=new double[n];
    double* x;
    x=new double[n];
    for(int i=0;i<n;i++)
        A[i]=new double[n];

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            fin>>A[i][j];
        }
    }

    Output(A,b,n);
    for(int i=0;i<n;i++){
        fin>>b[i];
    }
    //Прямой проход
    for(int k=0;k<n;k++){
        //Находим максимальный элемент в столбце
        int max_ind=k;
        for(int x=k+1;x<n;x++){
            if(fabs(A[x][k])>fabs(A[max_ind][k])){
                max_ind=x;
            }
        }
        if(A[max_ind][k]==0){
            cout<<"Матрица несовместная";
            return 1;
        }
        //Ставим его на k место
        swap_ij(A,n,k,max_ind);
        double temp=b[k];
        b[k]=b[max_ind];
        b[max_ind]=temp;

        //Обнуляем столбик
        for(int x=k+1;x<n;x++){
            b[x]=b[x]+b[k]*(-A[x][k]/A[k][k]);
            sum_i2j_with_c(A,k,x,-A[x][k]/A[k][k],n);
        }

    }

    Output(A,b,n);
    //Обратный проход
    for(int i=n-1;i>=0;i--){
        s=0;
        for(int j=i+1;j<n;j++){
            s+=A[i][j]*x[j];
        }
        x[i]=(b[i]-s)/A[i][i];
    }

    cout<<"\nОтвет (Прямой метод Гаусса):\n";
    for(int i=0;i<n;i++){
        cout<<x[i]<<'\n';
    }

    fin.close();
    for(int i=0;i<n;i++){
        delete [] A[i];
    }
    delete [] A;
    delete [] b;
    delete [] x;
    return 0;
}