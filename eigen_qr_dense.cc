#include <iostream>
#include <time.h> 
#include "toolkit/eigen-3.4.0/Eigen/Dense"
using namespace std;
using namespace Eigen;
VectorXd qr_AX_ex_b_solving() {
    clock_t startinit = clock();
    int n;
    cout<<"Введите размерность системы\n";
    cin>>n;
    MatrixXd A = MatrixXd::Ones(n,n);
    VectorXd b(n);
    for(int i=0;i<n;i++){
        b(i)=sin(i);
        for(int j=0;j<n;j++){
            if (i!=j && (i-j)!=0){
                A(i,j)=1.0/abs(i-j);
            }
        }
    }
    double eps=b.norm()*pow(10,-8);
    cout<<"Система сгенерирована,идет факторизация \n";
    /*cout<<"A=\n--------------------------------------------------\n"<<A<<endl;
    cout<<"b=\n--------------------------------------------------\n"<<b<<endl;*/
    HouseholderQR<MatrixXd> qr(n,n);
    qr.compute(A);
    MatrixXd R = qr.matrixQR().triangularView<Upper>();
    MatrixXd Q =  qr.householderQ();
    clock_t endinit = clock();
    clock_t start = clock();
    /*cout<<"Q=\n--------------------------------------------------\n"<<Q<<endl;
    //cout<<"R=\n--------------------------------------------------\n"<<R<<endl;
    cout<<"Q*R=\n--------------------------------------------------\n"<<Q*R<<endl;*/
    MatrixXd x=qr.solve(b);
    clock_t end = clock();
    //cout<<"x=\n--------------------------------------------------\n"<<x<<endl;
    VectorXd ax = A*x;
    //cout<<"A*x=\n--------------------------------------------------\n"<<ax<<endl;
    if((ax-b).norm()<eps){
        cout<<"Решение в норме\n";
    }
    else cout<<"pzdc";
    double timeinit = double(endinit - startinit) / CLOCKS_PER_SEC;
    double timeexec = double(end - start) / CLOCKS_PER_SEC;
    cout<<"Время инициализации = "<<timeinit<<" Сек"<<endl<<"Время исполнения = "<<timeexec<<" Сек"<<endl;
    return x;
}
int main(){
    qr_AX_ex_b_solving();
    return 0;
}
