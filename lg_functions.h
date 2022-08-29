
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DIM 100

// dichiaro funzioni utilizzate
void scegli_problema(int*, float*, float*);
float f(int, float);
float calcolaPolyInt(float x, int N, float nodi[], int scelta);
float calcolaPolyIntLin(int h, float x, float nodi[], int scelta);
float calcolaPolyIntQuad(int h, float x, float nodi[], int scelta);
void interpolazioneGlobale(float a, float b, int N, int punti, int scelta, float nodi[], float A[N][N]);
void interpolazioneCompositaLineare(float a, float b, int N, int punti, int scelta, float nodi[]);
void interpolazioneCompositaQuadratica(float a, float b, int N, int punti, int scelta, float nodi[]);
float disegnaGrafico(float a, float b, int punti, int scelta);
float calcolaPolyIntNewton(float x, int N, float nodi[], float A[N][N], int scelta);
void calcolaDiffDivise(int N, float nodi[], int scelta, float A[N][N]);


// definizione

float f(int scelta, float x){
  float fx;
  switch(scelta)
    {
    case 1: fx=sin(x); break;
    case 2: fx=sin(3*x); break;
    case 3: fx=exp(x); break;
    case 4: fx=x*cos(x); break;
    case 5: fx=1/(1+x*x); break;
    case 6: fx=fabs(x); break;
    case 7: fx=fabs(x); break;
    case 8: fx=sin(x*x); break;
    case 9: fx=fabs(sin(x*x)); break;
    case 10: fx=fabs(x)*sin(x); break;
    case 11: fx=fabs(x)*sin(x); break;
    case 12: fx=pow(x,5)+3*pow(x,4)+2*pow(x,3)-pow(x,2)-5*x+1; break;
    case 13: fx=log(x); break;
    }
  return fx;
}

void interpolazioneGlobale(float a, float b, int N, int punti, int scelta, float nodi[], float A[N][N])
{
    float deltaX, dx, m, xi, yi, yi_2;
    long double maxRes, res;
    int i;
    FILE *file;
    deltaX = (b-a)/(N); //Divide l'intervallo da studiare per il numero di sottointervalli equispaziati
    dx = (b-a)/(punti); //Questo valore rappresenza la distanza tra due punti che debbano essere rappresentati nel grafico
    file = fopen("interpolazioneGlobale.txt", "wt");
    for (i=0; i<=punti; i++)
    {
        xi = a + (dx*i);
        yi = calcolaPolyInt(xi, N, nodi, scelta); //Valuta il polinomio interpolatore nel punto xi
        yi_2 = calcolaPolyIntNewton(xi, N, nodi, A, scelta);
        //printf("rv: %f iv:_%f\n",f(scelta, xi), yi);
        fprintf(file, "%f  %f \n", xi, yi_2);
    }
    fclose(file);
    maxRes = 0;
    for(i=0; i<N; i++)
    {
        m = nodi[i] + (deltaX/2); //m è il punto medio di ciascun sottointervallo
        res = fabs(f(scelta,m)-calcolaPolyIntNewton(m, N, nodi, A, scelta));
        //printf("res : %Lf, value: %.5f\n", res, calcolaPolyIntNewton(xi, N, nodi, A, scelta));
        if(res>maxRes){
            maxRes = res;
        }
    }
    if(maxRes == 0 && res != res){ // ossia se res è sempre nan o -inf e quindi il polinomio diverge
        maxRes = res;
    }
    if(maxRes!=maxRes){
        printf("\nPer l'interpolazione globale, la norma infinito del vettore degli errori calcolati nel punto medio di ogni intervallo e' troppo grande, poichè il grado del polinomio è elevato e perciò il grafico del polinomio interpolatore non è plottato\n");
    }
    else{
        printf("\nPer l'interpolazione globale, la norma infinito del vettore degli errori calcolati nel punto medio di ogni intervallo e' %Lf\n", maxRes);
    }
}


void interpolazioneCompositaLineare(float a, float b, int N, int punti, int scelta, float nodi[])
{
    float deltaX, dx, m, xj, yj, maxRes, res;
    int i, j, acc;
    FILE *file;
    deltaX = (b-a)/N; //Divide l'intervallo da studiare per il numero di sottointervalli
    acc = (int)(punti/N); //Questo valore rappresenta i punti che devono essere rappresentati in ogni sottointervallo
    dx = deltaX/acc; //Questo valore rappresenza la distanza tra due punti che debbano essere rappresentati nel grafico
    file = fopen("interpolazioneCompositaLineare.txt", "wt");
    for (i=0; i<N; i++)
        for(j=0; j<=acc; j++)
        {
            xj = nodi[i] + (j*dx);
            yj = calcolaPolyIntLin(i, xj, nodi, scelta);
            fprintf(file, "%f  %f\n", xj, yj);
        }
    fclose(file);
    maxRes = 0;
    for(j=0; j<N; j++)
    {
        m = nodi[j] + (deltaX/2);
        res = fabs(f(scelta,m)-calcolaPolyIntLin(j, m, nodi, scelta));
        if(res>maxRes)
            maxRes = res;
    }
    printf("\nPer l'interpolazione composita lineare, la norma infinito del vettore degli errori calcolati nel punto medio di ogni intervallo e' %f\n", maxRes);
}

void interpolazioneCompositaQuadratica(float a, float b, int N, int punti, int scelta, float nodi[])
{
    float deltaX, dx, m, xj, yj, maxRes, res;
    int i, j, acc;
    FILE *file;
    deltaX = (b-a)/N; //Divide l'intervallo da studiare per il numero di sottointervalli
    acc = (int)(punti/N); //Questo valore rappresenta i punti che devono essere rappresentati in ogni sottointervallo
    dx = deltaX/acc; //Questo valore rappresenza la distanza tra due punti che debbano essere rappresentati nel grafico
    file = fopen("interpolazioneCompositaQuadratica.txt", "wt");
    for (i=0; i<N; i++)
        for(j=0; j<=acc; j++)
        {
            xj = nodi[i] + (j*dx); // sono i punti in cui valuto il polinomio
            yj = calcolaPolyIntQuad(i, xj, nodi, scelta);
            fprintf(file, "%f  %f\n", xj, yj);
        }
    fclose(file);
    maxRes = 0;
    for(j=0; j<N; j++)
    {
        m = nodi[j] + (deltaX/4); //I nodi su cui è costruito il polinomio sono due volte più fitti, quando si usa l'interpolazione quadratica
        res = fabs(f(scelta,m)-calcolaPolyIntQuad(j, m, nodi, scelta)); //Calcolo l'errore sul punto medio della prima metà del j-esimo intervallo
        if(res>maxRes)
            maxRes = res;
        m = nodi[j] + (3*deltaX/4);
        res = fabs(f(scelta,m)-calcolaPolyIntQuad(j, m, nodi, scelta)); //Calcolo l'errore sul punto medio della seconda metà del j-esimo intervallo
        if(res>maxRes)
            maxRes = res;
    }
    printf("\nPer l'interpolazione composita quadratica, la norma infinito del vettore degli errori calcolati nel punto medio di ogni intervallo e' %f\n", maxRes);
}

float disegnaGrafico(float  a, float b, int punti, int scelta)
{
    int i;
    double ymax;
    float xi, yi, dx;
    FILE *file;
    dx = (b-a)/(punti-1); //Questo valore rappresenza la distanza tra due punti che debbano essere rappresentati nel grafico
    file = fopen("grafico.txt", "wt");
    for(i=0; i<punti; i++)
    {
        xi = a+i*dx;
	yi=f(scelta,xi);
	if(fabs(yi)>ymax) // possibile errore
	  ymax = fabs(yi);
        fprintf(file, "%f  %f\n", xi, yi);
    }
    fclose(file);
    return ymax;
}

float calcolaPolyInt(float x, int N, float nodi[], int scelta)
{
    float l_i, p;
    int i, k;
    p=0;
    for(k=0; k <= N; k++)
    {
        l_i = 1;
        for (i = 0; i <= N; i++)
            if(i!=k)
                l_i *= (x-nodi[i])/(nodi[k]-nodi[i]);
        l_i *= f(scelta,nodi[k]);
        p+= l_i;
    }
    return p;
}

float calcolaPolyIntLin(int h, float x, float nodi[], int scelta)
{
    return f(scelta,nodi[h])+((f(scelta,nodi[h+1])-f(scelta,nodi[h]))/(nodi[h+1]-nodi[h]))*(x-nodi[h]); //Calcola il valore del polinomio interpolatore a tratti considerato nell'intervallo [x[h], x[h+1]]
}

float calcolaPolyIntQuad(int h, float x, float nodi[], int scelta)
{
    float l_i, p, halfInt, xk, xi;
    int i, k;
    p=0;
    halfInt = (nodi[h+1]-nodi[h])/2; // diventa /k se uso composito di grado k
    for(k=0; k <= 2; k++) // uso 3 nodi e quindi una parabola
    {
        xk = nodi[h] + k* halfInt;
        l_i = 1; // polinomi di base per lagrange
        for (i = 0; i <= 2; i++)
            if(i!=k)
            {
                xi = nodi[h] + i*halfInt;
                l_i*=(x-xi)/(xk-xi);
            }
        l_i *= f(scelta,xk);
        p+= l_i;
    }
    return p;
}

void calcolaDiffDivise(int N, float nodi[], int scelta, float A[N][N]){
    for(int i = 0; i <= N; i++){
        for(int j = 0; j <= N; j++){
            A[i][j] = 0;
        }
    }

    for(int i = 0; i <= N; i++){
        for(int j = i; j <= N; j++){
            if(i != 0){
                A[j][i] = (A[j][i-1] - A[j-1][i-1]) / (nodi[j] - nodi[j-i]);
            }
            else{
                A[j][i] = f(scelta, nodi[j]);
            }
        }
    }
}

// calcola pol di Lagrange con differenze divise
// usando l'algoritmo di Horner

float calcolaPolyIntNewton(float x, int N, float nodi[], float A[N][N], int scelta){
    float p = A[N][N] * (x - nodi[N-1]) + A[N-1][N-1];// valore polinomio
    for(int j = 2; j <= N; j++){
        p = p * (x - nodi[N-j]) + A[N-j][N-j];
    }

    return p;
}