#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <time.h>

using namespace std;

//From Numerical Recipes, third edition
struct Ranq1 { //Recommended generator for everyday use. The period is ~1.8E19, so it should not be used by an application that makes more than 1E12 calls.
	unsigned long long int v;
	Ranq1(unsigned long long int j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline unsigned long long int int64() {//Return 64-bit random integer
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); } //Return random double-precision floating value in the range 0. to 1
	inline unsigned int int32() { return (unsigned int)int64(); }  //Return 32-bit random integer
};

Ranq1 myran(time(0));

//Red de LxL
const int L=50;

//Tiempo final (si lo hago hasta un dado tiempo) o cantidad de avalanchas totales (si lo hago hasta cierta cantidad de avalanchas)
const int T=1E6;

//En la notación del enunciado: alpha
const double a=0.15;

int dinamica_lenta(double red[L+2][L+2], double &t, double &fuerzas) //recibe la red, el tiempo y un double para guardar la suma de fuerzas
{
    double f_max=0;
    int n;
    //recorro la red (sin contar los bordes) para encontrar el máximo (y aprovecho para calcular la suma de fuerzas)
    for(int i=1; i<L+1; i++)
    {
        for(int j=1; j<L+1; j++)
        {
            if(red[i][j]>f_max)
            {
                f_max=red[i][j];
                n=i*(L+2)+j;
            }
            fuerzas+=red[i][j];
        }
    }
    //avanzo el tiempo. Tomo k_1*V=1
    t+=(1-f_max);
    for(int i=1; i<L+1; i++)
    {
        for(int j=1; j<L+1; j++)
            red[i][j]+=(1-f_max);
    }
    return n;
}

void dinamica_rapida(double red[L+2][L+2], int &n, double &s)//recibe la red, el número correspondiente al sitio de red que llegó al umbral y un contador de cantidad de veces que se llama a la función (momento sísmico)
{
    int x, y;
    x=n/(L+2);
    y=n%(L+2);
    //Delta F (si llevo a cero la fuerza del sitio que se volvi� inestable: igualo a red[x][y])(si en cambio resto 1, igualo a 1)
    double dF=1;
    //disminuyo la fuerza sobre el sitio que llegó a su umbral
    red[x][y]-=dF;
    //aumento el momento sísmico en dF
    s+=dF;
    //aumento en (alpha*Delta F) la fuerza sobre los vecinos
    red[x-1][y]+=a*dF;
    red[x+1][y]+=a*dF;
    red[x][y-1]+=a*dF;
    red[x][y+1]+=a*dF;
    //Me fijo si alguno de los vecinos llegó al umbral
    if(red[x-1][y]>=1 && x-1>0)
    {
        int este=(x-1)*(L+2)+y;
        dinamica_rapida(red,este,s);
    }
    if(red[x+1][y]>=1 && x+1<(L+1))
    {
        int oeste=(x+1)*(L+2)+y;
        dinamica_rapida(red,oeste,s);
    }
    if(red[x][y-1]>=1 && y-1>0)
    {
        int norte=x*(L+2)+y-1;
        dinamica_rapida(red,norte,s);
    }
    if(red[x][y+1]>=1 && y+1<(L+1))
    {
        int sur=x*(L+2)+y+1;
        dinamica_rapida(red,sur,s);
    }
    return;
}

int main()
{
    //Tiempo
    double t=0;
    //cantidad de avalanchas
    unsigned int avalanchas=0;
    //suma de todas las fuerzas
    double fuerzas;
    //Red de sitios
    double red[L+2][L+2];// los 2 bordes en cada dirección son para sitios que nunca se descargan y sólo absorben (condiciones de borde abiertas)
    //archivo para guardar los valores de momento sísmico de cada terremoto
    ofstream datos_momento_sismico;
    datos_momento_sismico.open("65_datos_momento_sismico_L_50_a_15_T_a_1E6.dat");
    //archivo para guardar la suma de los f y los tiempos (lo comento cuando no quiero la evolución temporal)
    ofstream datos_fuerzas;
    datos_fuerzas.open("65_datos_fuerzas_L_50_a_15_T_a_1E6.dat");
    //inicializo los sitios de red con un double aleatorio entre 0 y 1
    for(int i=0; i<L+2; i++)
    {
        for(int j=0; j<L+2; j++)
        {
            red[i][j]=myran.doub();
        }
    }
    //Dejo transcurrir el tiempo
    while(avalanchas<T)
    {
        //cout << avalanchas << endl;
        fuerzas=0;
        datos_fuerzas << t << '\t';
        int n=dinamica_lenta(red,t,fuerzas);
        datos_fuerzas << fuerzas << '\n';
        //momento sísmico
        double s=0;
        dinamica_rapida(red,n,s);
        avalanchas++;
        //Guardo el momento sísmico
        datos_momento_sismico << t << '\t' << s << '\n';
    }
    datos_momento_sismico.close();
    datos_fuerzas.close();
    return 0;
}















