#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <random>  //Libreria que genera numeros aleatorios reales y enteros


#define T 0.3         //Valor temperatura
#define N 100         //Dimension de la matriz
#define PMC 300       //Pasos Monte Carlo      


double CalcularP(double s[N][N], int n, int m);    //Funcion que calcula el valor de P 


using namespace std;

int main(){

    int i, j, k, n, m, f, b;
    double random, x,p, mag;
    double s[N][N];


    ofstream pos, datos3, magn;
    magn.open("magnetizacion_T0,3.txt");
    pos.open("posiciones.txt");
    datos3.open("T_0,3(2).txt");


    mt19937 semilla(time(NULL));
    uniform_int_distribution<int> random_entero(0,N-1);
    uniform_real_distribution<double> random_real(0.0,1.0);


    //Inicializamos la matriz desordenada
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            random=random_real(semilla);
            if(random<=0.5){
                s[i][j]=-1.0;
            }
            else{
                s[i][j]=+1.0;
            }
        }
    }


    //Volcamos los datos iniciales en el fichero

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            datos3<<s[i][j]<<"\t";
        }
        datos3<<endl;
    }

    datos3<<endl;




    //Inicializamos la magnetizacion a 0
    mag=0.0;

    for(k=0;k<PMC;k++){

        for(i=0;i<N*N;i++){

            //Elegimos punto (n,m) aleatorio de la red

            n=random_entero(semilla);
            m=random_entero(semilla);

            //Evaluamos p=min[1,exp-(E/T)]

            p=CalcularP(s, n, m);

            //Generamos numero aleatorio uniforme en el intervalo [0,1], y si es <=p cambiamos el signo al elemento

            x=random_real(semilla);

            if(x<=p){
                s[n][m]=-s[n][m];
            }
            else    s[n][m]=s[n][m];

        }


        //Volcamos los datos de un PMC en el archivo
        for(b=0;b<N;b++){
            for(f=0;f<N;f++){
                datos3<<s[b][f]<<"\t";
            }
            datos3<<endl;
        }
        datos3<<endl;

        //Calculamos la magnetizacion de un PMC
        for(i=0;i<N-1;i++){
            for(j=0;j<N-1;j++){
                mag=mag+s[i][j];
            }
        }
        mag=abs(mag);
        mag=mag/(N*N);
        magn<<mag<<endl;

        mag=0.0;

        

    }


    pos.close();
    datos3.close();
    magn.close();


    return 0;
}


double CalcularP(double s[N][N], int n, int m){

    double E, p;

    //Se calcula E teniendo en cuenta las condiciones de contorno

    if((n==0)&&(m==0)){
        E=2.0*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
    }

    if((n==0)&&(m==N-1)){
        E=2.0*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
    }

    if((n==N-1)&&(m==0)){
        E=2.0*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }

    if((n==N-1)&&(m==N-1)){
        E=2.0*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }

    if((n==0)&&(m!=0)&&(m!=N-1)){
        E=2.0*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);
    }

    if((n==N-1)&&(m!=0)&&(m!=N-1)){
        E=2.0*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }

    if((m==0)&&(n!=0)&&(n!=N-1)){
        E=2.0*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }

    if((m==N-1)&&(n!=0)&&(n!=N-1)){
        E=2.0*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }

    else{
        E=2.0*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }


    p=min(1.0,exp(-E/T));


    return p;

}