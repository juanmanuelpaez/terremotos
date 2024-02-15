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
const int T=1E5;

//En la notación del enunciado: alpha
const double a=0.2;

//sigma para la distribución gaussiana de los umbrales
const double sigma=0.3;

int dinamica_lenta(double red[L+2][L+2], double &t, double &fuerzas, double umbral[L+2][L+2]) //recibe la red, el tiempo y un double para guardar la suma de fuerzas
{
    double df_min=2; //pongo 2 de valor inicial porque sé que todos los df son menores que 2
    int n;
    //recorro la red (sin contar los bordes) para encontrar el máximo (y aprovecho para calcular la suma de fuerzas)
    for(int i=1; i<L+1; i++)
    {
        for(int j=1; j<L+1; j++)
        {
            if((umbral[i][j]-red[i][j])<df_min)
            {
                df_min=umbral[i][j]-red[i][j];
                n=i*(L+2)+j;
            }
            fuerzas+=red[i][j];
        }
    }
    //avanzo el tiempo. Tomo k_1*V=1
    t+=df_min;
    for(int i=1; i<L+1; i++)
    {
        for(int j=1; j<L+1; j++)
            red[i][j]+=df_min;
    }
    return n;
}

void dinamica_rapida(double red[L+2][L+2], int &n, double &s, double umbral[L+2][L+2])//recibe la red, el número correspondiente al sitio de red que llegó al umbral y un contador de cantidad de veces que se llama a la función (momento sísmico)
{
    int x, y;
    x=n/(L+2);
    y=n%(L+2);
    //Delta F (si llevo a cero la fuerza del sitio que se volvió inestable: igualo a red[x][y])(si en cambio resto 1, igualo a 1)
    double dF=red[x][y];
    //Aumento el momento sísmico en dF
    s+=dF;
    //disminuyo la fuerza sobre el sitio que llegó a su umbral (resto 1 o llevo a 0)
    red[x][y]-=dF;
    //aumento en (alpha*Delta F) la fuerza sobre los vecinos
    red[x-1][y]+=a*dF;
    red[x+1][y]+=a*dF;
    red[x][y-1]+=a*dF;
    red[x][y+1]+=a*dF;
    //cambio el umbral del sitio que se había vuelto inestable
    while(true)
    {
        double aleatorio_1=1+(myran.doub()-0.5)*6*sigma;
        double aleatorio_2=myran.doub();
        if(aleatorio_2<(exp(-(aleatorio_1-1)*(aleatorio_1-1)/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma)) && aleatorio_1>red[x][y])
        {
            umbral[x][y]=aleatorio_1;
            break;
        }
    }
    //Me fijo si alguno de los vecinos llegó al umbral
    if(red[x-1][y]>=umbral[x-1][y] && x-1>0)
    {
        int este=(x-1)*(L+2)+y;
        dinamica_rapida(red,este,s,umbral);
    }
    if(red[x+1][y]>=umbral[x+1][y] && x+1<(L+1))
    {
        int oeste=(x+1)*(L+2)+y;
        dinamica_rapida(red,oeste,s,umbral);
    }
    if(red[x][y-1]>=umbral[x][y-1] && y-1>0)
    {
        int norte=x*(L+2)+y-1;
        dinamica_rapida(red,norte,s,umbral);
    }
    if(red[x][y+1]>=umbral[x][y+1] && y+1<(L+1))
    {
        int sur=x*(L+2)+y+1;
        dinamica_rapida(red,sur,s,umbral);
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
    //Red de sitios (contiene la fuerza de cada sitio)
    double red[L+2][L+2];// los 2 bordes en cada dirección son para sitios que nunca se descargan y sólo absorben (condiciones de borde abiertas)
    //Red de fuerzas umbrales correspondientes a cada sitio (incluyo los bordes para no tener índices distintos)
    double umbral[L+2][L+2];
    //archivo para guardar los valores de momento sísmico de cada terremoto
    ofstream datos_momento_sismico;
    datos_momento_sismico.open("40_datos_momento_sismico_gauss_2_L_50_a_2_T_a_1E5.dat");
    //archivo para guardar la suma de los f y los tiempos (lo comento cuando no quiero la evolución temporal)
    ofstream datos_fuerzas;
    datos_fuerzas.open("40_datos_fuerzas_gauss_2_L_50_a_2_T_a_1E5.dat");
    //inicializo los sitios de red con un double aleatorio entre 0 y 1 y un umbral aleatorio con distribución gaussiana de ancho 0.3 centrada en 1
    for(int i=0; i<L+2; i++)
    {
        for(int j=0; j<L+2; j++)
        {
            red[i][j]=myran.doub();
            while(true)
            {
                double aleatorio_1=1+(myran.doub()-0.5)*6*sigma;//elijo un número alrededor de 1 con 3 sigmas de ancho
                double aleatorio_2=myran.doub();
                if(aleatorio_2<(exp(-(aleatorio_1-1)*(aleatorio_1-1)/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma)) && aleatorio_1>red[i][j])
                {
                    umbral[i][j]=aleatorio_1;
                    break;
                }
            }
        }
    }
    //Dejo transcurrir el tiempo
    while(avalanchas<T)
    {
        //cout << avalanchas << endl;
        fuerzas=0;
        datos_fuerzas << t << '\t';
        int n=dinamica_lenta(red,t,fuerzas,umbral);
        datos_fuerzas << fuerzas << '\n';
        //momento sísmico
        double s=0;
        dinamica_rapida(red,n,s,umbral);
        avalanchas++;
        //Guardo el momento sísmico
        datos_momento_sismico << t << '\t' << s << '\n';
    }
    datos_momento_sismico.close();
    datos_fuerzas.close();
    return 0;
}



