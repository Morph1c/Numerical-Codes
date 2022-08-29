#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include </home/andrea/Scrivania/esameAN/lg_functions.h>


int main(void)
{
    float a, b, deltaX, eps, ymax;
    int i, N, punti, scelta, type_interpolation;
    float nodi[DIM];
    FILE *file;
    float *pa, *pb;
	int *pscelta;
    pa = &a;
    pb = &b;
    pscelta = &scelta;
    scegli_problema(pscelta,pa,pb);
	printf("\nInserisci il numero di sottointervalli(i.e. grado del polinomio per globale o numero sottointervalli per composita):\n");
    scanf("%d", &N);
    //printf("\nInserisci la precisione della rappresentazione grafica (piu' e' bassa, piu' e' preciso il grafico ossia il passo di distanza tra i punti):\n");
    //scanf("%f", &eps);  //Il valore inserito dall'utente rappresenta la distanza tra due punti rappresentati nel grafico
    eps = 0.01;

    int node_type;
    printf("Scegli il tipo di griglia\n[0] griglia equispaziata\n[1]Nodi chebyscev\n: ");
    scanf("%d", &node_type);
    // crea nodi e decide numero di punti in cui valutare pol
    punti = (int)((b-a)/eps); //Il valore "punti" corrisponde al numero di punti da rappresentare nel grafico
    deltaX = (b-a)/(N); //Divide l'intervallo da studiare per il numero di sottointervalli equispaziati
    for(i=0; i<=N; i++){
        if(node_type == 0){
            nodi[i] = a + i*deltaX;
        }
        else{
            nodi[i]=a+(b-a)*((-1*cos((2*i+1)*M_PI/(2*N+2)))+1)/2;
        }
    }
    printf("\nInserisci il tipo di interpolazione:\n[1]Globale(Con polinomio newton)\n[2]composita lineare\n[3]Composita Quadratica\n: ");
    scanf("%d", &type_interpolation);
    ymax = disegnaGrafico(a, b, punti,scelta);
    //printf("ymax:_%.4f\n", ymax);
    file = fopen("comando.gp", "wt");
    switch (type_interpolation){
        case 1:{
            float A[N][N];
            calcolaDiffDivise(N, nodi, scelta, A);
            interpolazioneGlobale(a, b, N, punti, scelta, nodi, A);
            fprintf(file, "set yrange [%.4f:%.4f];\nplot \'grafico.txt\' with lines, \'interpolazioneGlobale.txt\' using 1:2 with lines", -2*ymax, 2*ymax);
            break;
        }

        case 2:{
            interpolazioneCompositaLineare(a, b, N, punti,scelta, nodi);
            fprintf(file, "set yrange [%.4f:%.4f];\nplot \'grafico.txt\' with lines, \'interpolazioneCompositaLineare.txt\' with lines", -2*ymax, 2*ymax);
            break;
        }

        case 3:{
            interpolazioneCompositaQuadratica(a, b, N, punti,scelta, nodi);
            fprintf(file, "set yrange [%.4f:%.4f];\nplot \'grafico.txt\' with lines, \'interpolazioneCompositaQuadratica.txt\' with lines", -2*ymax, 2*ymax);
            break;
        }
    }
    //ymax = disegnaGrafico(a, b, punti,scelta); //disegnaGrafico, oltre a disegnare il grafico della funzione, restituisce la massima ordinata, in valore assoluto, che la funzione assume in [a,b]
    //file = fopen("comando.gp", "wt");
    //fprintf(file, "set yrange [%f:%f];\nplot \'grafico.txt\' with lines, \'interpolazioneGlobale.txt\' with lines, \'interpolazioneCompositaLineare.txt\' with lines, \'interpolazioneCompositaQuadratica.txt\' with lines", -ymax, ymax);
    fclose(file);
    system("gnuplot -p comando.gp");
}

void scegli_problema(int *scelta, float *a, float *b){
  printf("\n\n Scegli la funzione da approssimare:");
  printf("\n [1] f(x) = sin(x), x in [0,2pi]");
  printf("\n [2] f(x) = sin(3x), x in [0,2pi]");
  printf("\n [3] f(x) = exp(x), x in [0,1]");
  printf("\n [4] f(x) = xcos(x), x in [-2pi,2pi]");
  printf("\n [5] f(x) = 1/(1+x^2), x in [-5,5]");
  printf("\n [6] f(x) = |x|, x in [0,5]");
  printf("\n [7] f(x) = |x|, x in [-3,2]");
  printf("\n [8] f(x) = sin(x^2), x in [-5,5]");
  printf("\n [9] f(x) = |sin(x^2)|, x in [-5,5]");
  printf("\n [10] f(x) = |x|sin(x), x in [0,5]");
  printf("\n [11] f(x) = |x|sin(x), x in [-5,5]");
  printf("\n [12] f(x) = x^5 + 3x^4 + 2x^3 -x^2 -5x + 1, x in [-5,5]");
  printf("\n [13] f(x) = log(x), x in [1,1.5]");
  printf("\n Scegli la funzione: ");
  scanf("%d",scelta);

  int int_pre;
  printf("Vuoi scegliere un'intervallo o usare quello preimpostato?(1 = si/ 0 = no): ");
  scanf("%d", &int_pre);
  while((*scelta<1)&&(*scelta>13))
    {
      printf("\n Il valore inserito non e' corretto.");
      printf("\n Scegli la funzione: ");
      scanf("%d",scelta);
    }

  if(int_pre == 1){
    printf("\nInserisci gli estremi: ");
    scanf("%f", a);
    scanf("%f", b);
  }

  else{
  switch (*scelta)
    {
    case 1: *a=0; *b=2*acos(-1.); break;
    case 2: *a=0; *b=2*acos(-1.); break;
    case 3: *a=0; *b=1; break;
    case 4: *a=-2*acos(-1.); *b=2*acos(-1.); break;
    case 5: *a=-5; *b=5; break;
    case 6: *a=0; *b=5; break;
    case 7: *a=-3; *b=2; break;
    case 8: *a=-5; *b=5; break;
    case 9: *a=-5; *b=5; break;
    case 10: *a=0; *b=5; break;
    case 11: *a=-5; *b=5; break;
    case 12: *a=-5; *b=5; break;
    case 13: *a=1; *b=1.5; break;
    }
  }
  return;
}




