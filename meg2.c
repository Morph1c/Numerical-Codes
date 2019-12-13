/*
************************************************
* Questo programma esegue l'algoritmo di 
* Eliminazione di Gauss per poter risolvere
* un sistema lineare di N equazioni e stampa
* a video il vettore x che risolve l'equazione
* matriciale Ax = b
************************************************
* Scritto da: Andrea Toccaceli, 
* 13/12/2019(Ultima modifica)
*/

#include <stdlib.h>
#include <stdio.h>

// Prototipo funzioni
void acquisisciMatriceQuadrata(int n, float A[n][n]);
void acquisiciVettore(int n, float V[n]);
void generaMatriceSistema(int n, float A[n][n], float V[n], float M[n+1][n+2]); // M è una matrice (n+1) x (n+2)
void scambiaRiga(int n, float M[n+1][n+2], int j);
void normalizzaRiga(int i, int n, float M[n+1][n+2], float coefficente_pivot);
void Gauss(int n, float M[n+1][n+2]);
void stampaSoluzione(int n, float T[n+1][n+2], float x_tilde[n]);// T è una matrice triangolare superiore generata dall'eliminazione

// Main
int main(void){

	// Acquisisco il numero di incognite del sistema
	int n;
	printf("Inserisci il numero di incognite: ");
	scanf("%d", &n);
	printf("\n");

	// Acquisisco i coefficenti del sistema in una matrice n x n, 
	// e acquisisco il vettore b
	float A[n][n];
	float B[n];
	acquisisciMatriceQuadrata(n, A);
	acquisiciVettore(n, B);

	// Costruisco una matrice (n+1) x (n+2) per tener conto
	// del numero di colonna e di riga
	float M[n+1][n+2];
	generaMatriceSistema(n, A, B, M);

	// esegue MEG su M
	Gauss(n, M);
	printf("\n");

	// Stampa il vettore x soluzione dell'equazione Ax = b
	float x_tilde[n];
	stampaSoluzione(n, M, x_tilde);

	return 0;

}

// Acquisisce una matrice quadrata
void acquisisciMatriceQuadrata(int n, float A[n][n]){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("Inserisci A[%d][%d] = ", i, j);
			scanf("%f", &A[i][j]);
			printf("\n");
		}
	}
}

// Acquisice un vettore
void acquisiciVettore(int n, float V[n]){
	for(int i = 0; i < n; i++){
		printf("Inserisci V[%d] = ", i);
		scanf("%f", &V[i]);
		printf("\n");
	}
}


// Genera la matrice su cui eseguire MEG
void generaMatriceSistema(int n, float A[n][n], float V[n], float M[n+1][n+2]){
	for(int i = 0; i < (n + 1); i++){
		for(int j = 0; j < (n + 2); j++){// j va da 0 a n + 1
			if(i == 0){// alla prima riga(riga 0) aggiungo la numerazione 1, .., n per numerare le colonne
				M[i][j] = j;
			}
			else if(j == 0){// all'inizio di ogni riga metto il numero di riga
				M[i][j] = i;
			}
			else if(j > 0 && j < (n + 1)){// j compreso tra 1 ed n
				M[i][j] = A[i-1][j-1];
			}
			else{// j == n+1
				M[i][j] = V[i-1];
			}
		}
		printf("\n");
	}
}

/* Algortimo per scambiare righe e trovare la riga con l'elemento pivot più grande
** int j é il numero di riga da cui partire per effetuare la ricerca e scambiare con
** la riga candindato > j
** int j è anche la colonna dove si trova l'elemnento da scambiare
*/
void scambiaRiga(int n, float M[n+1][n+2], int j){
	float max_pivot = 0; 
	int pivot_index = 0;
	// ricerca massimo elemento tra le righe
	for(int m = j + 1; m < n + 1; m++){
		if(M[m][j] > max_pivot){
			max_pivot = M[m][j];
			pivot_index = m;
		}
	}

	// A questo punto so che devo scambiare la riga j con la riga pivot_index
	if(max_pivot > 0 && pivot_index > 0){// se il programma trova elementi da scambiare scambia altrimenti esci
		float temp = 0;
		for(int l = 0; l < n + 2; l++){//poichè devo scambiare anche il numero di riga nella colonna 0
			temp = M[j][l];
			M[j][l] = M[pivot_index][l];
			M[pivot_index][l] = temp;
		}
	}
}

// Prende una riga e la normalizza rendendo 1 il coefficente del termine M[i][i]
void normalizzaRiga(int i, int n, float M[n+1][n+2], float coefficente_pivot){
	for(int k = 1; k < n + 2; k++){
		M[i][k] = M[i][k] / coefficente_pivot;
	}
}


// Esegue l'algoritmo di eliminazione di Gauss
void Gauss(int n, float M[n+1][n+2]){
	float coefficente_pivot;
	float coefficente_primo_elemento;

	for(int i = 1; i < n ; i++){// escludo la prima riga dove ci sono gli indici colonna
		coefficente_pivot = M[i][i];
		if(coefficente_pivot == 0){
			scambiaRiga(n, M, i);
			coefficente_pivot = M[i][i];
		}
		if(coefficente_pivot == 0){
			break;
		}
		normalizzaRiga(i, n, M, coefficente_pivot);
		for(int j = i + 1; j < n + 1; j++){
			coefficente_primo_elemento = M[j][i];
			for(int k = i; k < n + 2; k++){
				M[j][k] = M[j][k] - (M[i][k]*coefficente_primo_elemento);
			}
	    }
	}
}


void stampaSoluzione(int n, float T[n+1][n+2], float x_tilde[n]){
	float somma;
	int stampa_soluzione = 0;
	for(int i = n; i  > 0; i--){
		somma = 0;
		if(T[i][i] == 0){
			printf("Il tuo sistema o ha infinite soluzioni o nessuna\n");
			stampa_soluzione = 1;
			break;
		}
		for(int j = i; j < n; j++){
			somma = somma + (x_tilde[j] * T[i][j + 1]);
		}
		somma = (T[i][n + 1] - somma) / T[i][i];
		x_tilde[i - 1] = somma;
	}

	if(stampa_soluzione == 0){
		for(int i = 0; i < n; i++){
			printf("x_%d = %f\n", i+1, x_tilde[i]);
		}
	}

}
