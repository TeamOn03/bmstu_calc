#include <iostream>
#include <fstream>
using namespace std;

double skalar_product(double* a, double* b, int n){
    double ans=0;
    for(int i=0;i<n;i++){
        ans+=a[i]*b[i];
    }
    return n;
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

//eto govno ne rabotaet
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

void proj(double* b, double* a,double* temp,int n){
    double c1=skalar_product(a,b,n);
    double c2=skalar_product(b,b,n);
    for(int i=0;i<n;i++){
        temp[i]=b[i]*c1/c2;
    }
}

int main(){
    setlocale(LC_ALL, "Russian");
    int n;
    double s;
    cout<<"Введите размер системы: ";
    cin>>n;
    cout<<"Введите матрицу системы:\n";
    double** A;
    A=new double*[n];
    double** Q;
    Q=new double*[n];
    double** R;
    R=new double*[n];
    double* b;
    b=new double[n];
    double* x;
    x=new double[n];
    double* temp;
    temp=new double[n];


    for(int i=0;i<n;i++){
        A[i]=new double[n];
        Q[i]=new double[n];
        R[i]=new double[n];
    }

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>A[j][i];
        }
    }



    cout<<"Введите свободные коэффициенты:\n";
    for(int i=0;i<n;i++){
        cin>>b[i];
    }

    Output(A,b,n);
    //Грамм-Шмидт
    for(int i=0;i<n;i++){
        Q[i]=A[i];
        for(int k=0;k<i;k++){
            proj(Q[k],A[i],temp,n);
            cout<<"\nПроекция "<<k<<": ";
            for(int j=0;j<n;j++){//Q[i]=Q[i]-temp
                cout<<temp[j]<<" ";
                Q[i][j]=Q[i][j]-temp[j];
            }
        }
    }

    Output(Q,b,n);

    for(int i=0;i<n;i++){
        delete [] A[i];
        delete [] Q[i];
        delete [] R[i];
    }
    delete [] Q;
    delete [] R;
    delete [] A;
    delete [] b;
    delete [] x;
    delete [] temp;
    return 0;
}