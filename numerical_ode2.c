/*
* Implementazione metodi ad 1 passo per la soluzione
* di ODE
* metodi impliciti:
* -Eulero in avanti
* -Eulero modificato
* -Heun
* -Runge Kutta
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include </home/andrea/Scrivania/plpc/diretti.c>

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// Prototipo funzioni
int eulero_avanti(int N, double X[N + 1][3], double h, int f);
int eulero_modificato(int N, double X[N + 1][3], double h, int f);
int heun(int N, double X[N + 1][2], double h, int f);
int runge_kutta(int N, double X[N + 1][3], double h, int f);
int eulero_indietro(int N, double X[N + 1][3], double h, int f);
int eulero_cromer(int N, double X[N + 1][3], double h, int f);
int newton(int f, int f2, double t, double yN, double vN, int N, int i, double X[N + 1][3]);
int jacobian(double j_valued[2][2], int mode_1, int mode_2, double t, double yN, double vN, double y_k, double v_k);
double campo_vettoriale_modificato(int mode, int i, double value_n, double t, double y, double v);
int gauss_2d(double sol[2], double M[2][2], double b[2]);
double campo_vettoriale(int f, double t, double y, double v);
double soluzione_esatta(int f, double t);

// VARIABILI GLOBALI

double h;

// MAIN

int main(void){
	FILE *file;
	double y0, t0, v0; // dato iniziale
	double T; // tempo finale
	int N; // numero di iterazioni

	int metodo;
	printf("Scegli il metodo di risoluzione:\n[1] Eulero in Avanti\n[2] Eulero modificato\n[3] Heun\n[4] Runge-Kutta (4° ordine)\n[5] Eulero all'indietro\n[6] Eulero-Cromer\n: ");
	scanf("%d", &metodo);

	int f;
	printf("Campo vettoriale:\nCaso R2\n[1] y'' = -w^2y\n[2] y'' = -w^2y -ay'\n[3] y'' = -sin(y) -y'\n[4] y'' = -w^2 + sin(t)\n[5, 6] Preda-Predatore \ny1' = a*y1 -b*y1*y2\ny2' =  -g*y2 + d*y1*y2\n[7, 8] Van der Pol\ny1'= y2 -y1^3 -y1\ny2' = -y1\n : "); // devo aggiungere tutte le ODE su cui simulare i metodi
	scanf("%d", &f);

	int dati_manuali;
	printf("\nDati insireti manualmente o preimpostati?(manuali = 1/preimpostati = 0): ");
	scanf("%d", &dati_manuali);
	if(dati_manuali == 1){
		printf("Inserisci i dati iniziali (v0, y0, t0): ");
		scanf("%lf", &v0);
		scanf("%lf", &y0);
		scanf("%lf", &t0);
		

		
	}

	else{
		switch(f){
			case 1: y0 = 1; v0 = 2; t0 = 0; T = 5; break;
			case 2: y0 = 1; v0 = 2; t0 = 0; T = 5; break;
			case 3: y0 = 1; v0 = 1; t0 = 0; T = 5; break;
			case 4: y0 = 1; v0 = 2; t0 = 0; T = 5; break;
			case 5: v0 = 80; y0 = 30; t0 = 0; T = 30; break;
			case 6: v0 = 1; y0; 1; t0 = 0; T = 5; break;
		}
	}
	
	printf("Inserisci il limite temporale: ");
	scanf("%lf", &T);
	//printf("Inserisci il passo di discretizzazione: ");
	//scanf("%lf", &h);
	printf("Inserisci numero di iterazioni(N): ");
	scanf("%d", &N);

	// Inizializzo matrice contenente tempi, valori e valuto passo di discretizzazione
	h = ((T + t0) - t0) / N;
	//T = h*N;
	double X[N + 1][3]; // la prima colonna per i tempi e la seconda per i valori di x1 e la terza per i valori di x2
	X[0][0] = t0;
	X[0][1] = y0;
	X[0][2] = v0;

	// Risolvo la ODE con uno dei metodi
	int out;
	switch(metodo){
		case 1: out = eulero_avanti(N, X, h, f); break; // il metodo sovrascrive la matrice X con i valori trovati
		case 2: out = eulero_modificato(N, X, h, f); break;
		case 3: out = eulero_modificato(N, X, h, f); break;
		case 4: out = runge_kutta(N, X, h, f); break;
		case 5: out = eulero_indietro(N, X, h, f); break;
		case 6: out = eulero_cromer(N, X, h, f); break;
	}

	// plotto i risultati e norma infinito errori
	file = fopen("graficode.txt", "wt");
	for(int i = 0; i < N + 1; i++){
		fprintf(file, "%lf %lf %lf\n", X[i][0], X[i][1], X[i][2]);
	}
	fclose(file);

	file = fopen("comandode1.gp", "wt");
	fprintf(file, "plot \'graficode.txt\' using 1:2 with lines title '(t, y)', \'graficode.txt\' using 1:3 with lines title '(t, v)'");
	fclose(file);

	file = fopen("comandode2.gp", "wt");
	fprintf(file, "\nplot \'graficode.txt\' using 2:3 with lines");
	fclose(file);

	system("gnuplot -p comandode1.gp");
	system("gnuplot -p comandode2.gp");

	return 0;
}


// EULERO AVANTI
int eulero_avanti(int N, double X[N + 1][3], double h, int f){
	int f2 = 0; // è il valore che mi dice che se ho preda predatore o van der pol devo scegliere il campo vettoriale diversamente per la 2 equ
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1;
	}
	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		X[i][1] = X[i - 1][1] + h*campo_vettoriale(f2, X[i - 1][0], X[i - 1][1], X[i - 1][2]); // modifico la y o y1
		X[i][2] = X[i - 1][2] + h*campo_vettoriale(f, X[i - 1][0], X[i - 1][1], X[i - 1][2]); // modifico la v o y2
	}

	return 0;
}

// EULERO MODIFICATO
int eulero_modificato(int N, double X[N + 1][3], double h, int f){
	double half_value_v;
	double half_value_y;
	int f2;
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1;
	}
	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		half_value_y = X[i - 1][1] + (h/2)*campo_vettoriale(f2, X[i - 1][0] + 0.5, X[i - 1][1], X[i - 1][2]);
		half_value_v = X[i - 1][2] + (h/2)*campo_vettoriale(f, X[i - 1][0] + 0.5, X[i - 1][1], X[i - 1][2]);//f(numero_campo_vett, tn + 0.5, yn, vn)
		X[i][1] = X[i - 1][1] + h*campo_vettoriale(f2, X[i - 1][0] + 0.5, half_value_y, half_value_v); // aggiorno y
		X[i][2] = X[i - 1][2] + h*campo_vettoriale(f, X[i - 1][0] + 0.5, half_value_y, half_value_v); // aggiorno v
	}

	return 0;
}

//HEUN
int heun(int N, double X[N + 1][2], double h, int f){
	double Y[N+1][3]; // valori con eulero esplicito
	Y[0][0] = X[0][0];
	Y[0][1] = X[0][1];
	Y[0][2] = X[0][2];
	int f2;
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1;
	}
	int nll = eulero_avanti(N, Y, h, f);

	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		X[i][1] = X[i - 1][1] + (h/2)*(campo_vettoriale(f2, X[i - 1][0], X[i - 1][1], X[i - 1][2]) + campo_vettoriale(0, Y[i][0], Y[i][1], Y[i][2]));
		X[i][2] = X[i - 1][2] + (h/2)*(campo_vettoriale(f, X[i - 1][0], X[i - 1][1], X[i - 1][2]) + campo_vettoriale(f, Y[i][0], Y[i][1], Y[i][2]));
	}

	return 0;
}

// RUNGE KUTTA

int runge_kutta(int N, double X[N + 1][3], double h, int f){
	int f2;
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1; // ossia esso è gia un sistema di 2 equ e quindi il campo vettoriale cambia
	}
	double k1, k2, k3, k4;
	double c1, c2, c3, c4;
	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		k1 = campo_vettoriale(f2, X[i - 1][0], X[i - 1][1], X[i - 1][2]);
		c1 = campo_vettoriale(f, X[i - 1][0], X[i - 1][1], X[i - 1][2]);

		k2 = campo_vettoriale(f2, X[i - 1][0] + (h/2), X[i - 1][1] + (h/2)*k1, X[i - 1][2] + (h/2)*c1);
		c2 = campo_vettoriale(f, X[i - 1][0] + (h/2), X[i - 1][1] + (h/2)*k1, X[i - 1][2] + (h/2)*c1);

		k3 = campo_vettoriale(f2, X[i - 1][0] + (h/2), X[i - 1][1] + (h/2)*k2, X[i - 1][2] + (h/2)*c2);
		c3 = campo_vettoriale(f, X[i - 1][0] + (h/2), X[i - 1][1] + (h/2)*k2, X[i - 1][2] + (h/2)*c2);

		k4 = campo_vettoriale(f2, X[i][0], X[i - 1][1] + h*k3, X[i - 1][2] + h*c3);
		c4 = campo_vettoriale(f, X[i][0], X[i - 1][1] + h*k3, X[i - 1][2] + h*c3);

		X[i][1] = X[i - 1][1] + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
		X[i][2] = X[i - 1][2] + (h/6)*(c1 + 2*c2 + 2*c3 + c4);
	}

	return 0;
}

// EULERO INDIETRO
int eulero_cromer(int N, double X[N + 1][3], double h, int f){
	int f2 = 0; // è il valore che mi dice che se ho preda predatore o van der pol devo scegliere il campo vettoriale diversamente per la 2 equ
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1;
	}
	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		//X[i][1] = X[i - 1][1]; // modifico la y o y1
		//X[i][2] = X[i - 1][2]; // modifico la v o y2
		//newton(f, f2, X[i][0], X[i][1], X[i][2], N, i, X);
		X[i][2] = X[i - 1][2] + h*campo_vettoriale(f, X[i - 1][0], X[i - 1][1], X[i - 1][2]);
		X[i][1] = X[i - 1][1] + h*campo_vettoriale(f2, X[i - 1][0], X[i - 1][1], X[i][2]); // modifico la y o y1
	}

	return 0;
}

// EULERO INDIETRO
int eulero_indietro(int N, double X[N + 1][3], double h, int f){
	int f2 = 0; // è il valore che mi dice che se ho preda predatore o van der pol devo scegliere il campo vettoriale diversamente per la 2 equ
	if(f < 5){
		f2 = 0;
	}
	else{
		f2 = f + 1;
	}
	for(int i = 1; i < N + 1; i++){
		X[i][0] = X[0][0] + i*h; // aggiorno tempo
		X[i][1] = X[i - 1][1]; // modifico la y o y1
		X[i][2] = X[i - 1][2]; // modifico la v o y2
		newton(f, f2, X[i][0], X[i][1], X[i][2], N, i, X);
	}

	return 0;
}

int newton(int f, int f2, double t, double yN, double vN, int N, int i, double X[N + 1][3]){
	double F_k[2]; // F(X_k)
	F_k[0] = 1;
	F_k[1] = 1;
	double sol[2];
	double j[2][2];
	double eps = 0.001;
	while(fabs(max(F_k[0], F_k[1])) > eps){
		//printf("\nrisolvo newton per l'iterazione %d\n, F_k[0] = %lf, F_k[1] = %lf\n, sum = %lf", i, F_k[0], F_k[1], F_k[0] + F_k[1]);
		jacobian(j, f2, f, t,  yN, vN, X[i][1], X[i][2]); // dF1/dx, con x = y_n, y = v_n 
		F_k[0] = campo_vettoriale_modificato(f2, 1, yN, t , X[i][1], X[i][2]);
		F_k[1] = campo_vettoriale_modificato(f, 2, vN, t, X[i][1], X[i][2]);
		//printf("\nrisolvo newton per l'iterazione %d\n, F_k[0] = %lf, F_k[1] = %lf\n, sum = %lf", i, F_k[0], F_k[1], F_k[0] + F_k[1]);
		gauss_2d(sol, j, F_k);
		//MEG(2, j, F_k, sol);
		X[i][1] = X[i][1] - sol[0];
		X[i][2] = X[i][2] - sol[1];  
	}

	return 0;
}

int jacobian(double j_valued[2][2], int mode_1, int mode_2, double t, double yN, double vN, double y_k, double v_k){
	j_valued[0][0] = (campo_vettoriale_modificato(mode_1, 1, yN, t, y_k + h, v_k) - campo_vettoriale_modificato(mode_1, 1, yN, t, y_k - h, v_k)) / (2*h);//dF1/dx
	j_valued[0][1] = (campo_vettoriale_modificato(mode_2, 2, vN, t, y_k + h, v_k) - campo_vettoriale_modificato(mode_2, 2, vN, t, y_k - h, v_k)) / (2*h);
	j_valued[1][0] = (campo_vettoriale_modificato(mode_1, 1, yN, t, y_k, v_k + h) - campo_vettoriale_modificato(mode_1, 1, yN, t, y_k, v_k - h)) / (2*h);
	j_valued[1][1] = (campo_vettoriale_modificato(mode_2, 2, vN, t, y_k, v_k + h) - campo_vettoriale_modificato(mode_2, 2, vN, t, y_k, v_k - h)) / (2*h);

	return 0;

}

double campo_vettoriale_modificato(int mode, int i, double value_n, double t, double y, double v){
	double value;
	if(i == 1){
		value = y;
	}
	else{
		value = v;
	}

	return value - value_n - h * campo_vettoriale(mode, t, y, v); // sarebbe F ottenuta dal metodo implicito
}

int gauss_2d(double sol[2], double M[2][2], double b[2]){
	double k = M[1][0] / M[0][0];
	sol[1] = (b[1] - (k * b[0])) / (M[1][1] - (k*M[0][1]));
	sol[0] = (b[0] - M[0][1]*sol[1]) / M[0][0];

	return 0;
}

// CAMPO VETTORIALE

double campo_vettoriale(int f, double t, double y, double v){
	double omega = 2;
	double alpha = 2;
	double alpha2 = 0.25;
	double beta = 0.01;
	double gamma = 1;
	double delta = 0.01;
	double fx;
  	switch(f)
    	{
    	case 0: fx = v; break;
		case 1: fx = -1*(omega*omega)*y; break;
		case 2: fx= -1*(omega*omega)*y - alpha*v; break;
		case 3: fx = -1*sin(y)-1*v; break;
		case 4: fx = -1*(omega*omega)*y + sin(t); break;
		case 5: fx = -1*gamma*v + delta*y*v; break; // y2
		case 6: fx = alpha2*y - beta*y*v; break; // y1
		case 7: fx = -1*y; break; // y2
		case 8: fx = v -1*(y*y*y) -y; break; // y1
    	}
  return fx;
}









