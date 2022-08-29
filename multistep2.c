//
// Metodi multistep
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DIM 10000
#define eps 0.001

// Variabili globali
FILE *file;
double lambda;
double k;
double m;
int T;

// prototipo funzione
void   adbash1(int,double[DIM][2],double,int);
void   adbash2(int,double[DIM][2],double,int);
void   adbash3(int,double[DIM][2],double,int);
void   adbash4(int,double[DIM][2],double,int);
double admou(double,double[DIM][2],double,int,int,int);
double bdf(double,double[DIM][2],double,int,int,int);
double derivata(double,double[DIM][2],double,int,int,int,int);
int    esempio();
double f(double,double,int);
void   grafico(double[DIM][2],int,int);
//void   grafico1(double[DIM][2],int,int);
double newton(double,double[DIM][2],double,int,int,int,int);
double phi(double,double[DIM][2],double,int,int,int,int);
void   rk4(double[DIM][2],double,int,int);
void grafico1(double M[DIM][2], int N, int scelta);



int main(void){

	//definizione variabili
	int funzione,met,N;
	int exit = 0;
	double h, x0;
	double M[DIM][2];
	int metodo,step;

	printf("\n\nMETODI MULTISTEP 1D\n\n");

	while (exit==0) {

		//scegliamo esempio su cui testare il programma
		funzione=esempio();
		
		//prendo dati in input 
	    printf("\nDATI IN INPUT\n\nInserire il tempo iniziale del sistema:\nt0= ");
	    scanf("%lf", &M[0][0]);
        printf("\n\nInserire il dato iniziale del sistema:\ny0= ");
        scanf("%lf", &M[0][1]);
        do{
	     	printf("\n\nInserire il passo di discretizzazione iniziale:\nh= ");
	     	scanf("%lf", &h);
	    } while (h<=0);
	    do{
	    	printf("\n\nInserire il numero dei passi:\nN= ");
			scanf("%d", &N);
	    } while((N<1)||(N>DIM-2));


		//scelgo il metodo da applicare
		do{
			printf("\n\n\n\nSCELTA DEL METODO NUMERICO\n\n1) Adams-Bashforth\n2) Adams-Moulton\n3) BDF\nSelezionare il metodo da applicare: ");
			scanf("%d", &metodo);
		} while ((metodo<1) || (metodo>3));

		do{
			printf("\n\n\nSelezionare il numero di step del metodo (1,2,3,4):\nstep= ");
			scanf("%d", &step);
		} while ((step<1) || (step>4));

		
 
 		// inizializzo i primi ord punti per far partire il metodo
 		for(int i =1; i < step; i++){
 			M[i][0]=M[i-1][0]+h; //questa non c'era ma i primi tempo non li devi mettere?
 			rk4(M, h, i, funzione);
 		}

		switch(metodo)
		{
			case 1:{

				switch(step)
				{
					case 1:{//adams-bashforth 1-step
						for(int n = step ;n <= N; n++){
							M[n][0]=M[n-1][0]+h;
							adbash1(n, M, h, funzione);
							//printf("h: %.2f at %d\n", h, n);
						}
						break;
					}

					case 2:{//adams-bashforth 2-step
						for(int n = step ;n <= N; n++){
							M[n][0]=M[n-1][0]+h;
							adbash2(n, M, h, funzione);
							//printf("h: %.2f at %d\n", h, n);
						}
						break;
					}
					case 3:{//adams-bashforth 3-step
						for(int n = step ;n <= N; n++){
							M[n][0]=M[n-1][0]+h;
							adbash3(n, M, h, funzione);
							//printf("h: %.2f at %d\n", h, n);
						}
						break;
					}
					case 4:{//adams-bashforth 4-step
						for(int n = step ;n <= N; n++){
							M[n][0]=M[n-1][0]+h;
							adbash4(n, M, h, funzione);
							//printf("h: %.2f at %d\n", h, n);
						}
						break;
					}
				}
				break;
			}
		
		case 2:{//adams-moulton
	    		for(int n = step ;n <= N; n++){
	    			M[n][0]=M[n-1][0]+h;
	    			//uso metodo gear per x0 con eulero esplicito
                    x0=M[n-1][1]+h*f(M[n-1][0],M[n-1][1],funzione);
	    			//con newton uso il metodo di adams-moulton
	    			//con il punto di partenza trovato con il metodo gear
	    			M[n][1]=newton(x0,M,h,n,metodo,step,funzione);
	    		}
				break;
			}

		case 3:{//bdf
	    		for(int n = step ;n <= N; n++){
	    			M[n][0]=M[n-1][0]+h;
	    			//uso metodo gear per x0 con eulero esplicito
                    x0=M[n-1][1]+h*f(M[n-1][0],M[n-1][1],funzione);
	    			//con newton uso il metodo di adams-moulton
	    			//con il punto di partenza trovato con il metodo gear
	    			M[n][1]=newton(x0,M,h,n,metodo,step,funzione);
	    		}
				break;

		}

	    }


        //grafico la soluzione trovata
        grafico1(M, N, funzione);
        printf("\n\n\n0) Ripetere il programma\n1) Uscire dal programma\n\nSelezionare una delle due opzioni: ");
        scanf("%d",&exit);
        printf("\n\n\n\n");

    }

return 0;
}




void adbash1(int n, double M[DIM][2], double h, int scelta){
	M[n][1] = M[n-1][1] + h*(f(M[n-1][0], M[n-1][1], scelta));
}

void adbash2(int n, double M[DIM][2], double h, int scelta){
	M[n][1] = M[n-1][1] + (h/2) * (3*f(M[n-1][0], M[n-1][1], scelta) - f(M[n -2][0], M[n- 2][1], scelta));
}

void adbash3(int n, double M[DIM][2], double h, int scelta){
	M[n][1] = M[n-1][1] + (h/12) * (23*f(M[n-1][0], M[n-1][1], scelta) - 16*f(M[n -2][0], M[n- 2][1], scelta) + 5*f(M[n-3][0], M[n-3][1], scelta));
}

void adbash4(int n, double M[DIM][2], double h, int scelta){
	M[n][1] = M[n-1][1] + (h/24) * (55*f(M[n-1][0], M[n-1][1], scelta) - 59*f(M[n -2][0], M[n- 2][1], scelta) + 37*f(M[n-3][0], M[n-3][1], scelta) - 9*f(M[n-4][0], M[n-4][1], scelta));
}


double admou(double x, double M[DIM][2], double h, int n, int step, int scelta){
	double y;

	switch(step){

		case 1:{
			y=M[n-1][1]+h*0.5*(f(M[n-1][0],M[n-1][1],scelta)+f(M[n][0],x,scelta))-x; // regola dei trapezi
			break;
		}
		case 2:{
			y=M[n-1][1]+h*(5*f(M[n][0],x,scelta)+8*f(M[n-1][0],M[n-1][1],scelta)-f(M[n-2][0],M[n-2][1],scelta))/12-x;
			break;
		}
		case 3:{
			y=M[n-1][1]+h*(9*f(M[n][0],x,scelta)+19*f(M[n-1][0],M[n-1][1],scelta)-5*f(M[n-2][0],M[n-2][1],scelta)+f(M[n-3][0],M[n-3][1],scelta))/24-x;
			break;
		}
		case 4:{
            y=M[n-1][1]+h*(251.*f(M[n][0],x,scelta)+646.*f(M[n-1][0],M[n-1][1],scelta)-264.*f(M[n-2][0],M[n-2][1],scelta)+106.*f(M[n-3][0],M[n-3][1],scelta)-19.*f(M[n-4][0],M[n-4][1],scelta))/720.-x;
			break;
		}
	}

	

	return y;
}



double bdf(double x, double M[DIM][2], double h, int n, int step, int scelta){
	double y;

	switch(step){

		case 1:{
			y=M[n-1][1]+h*f(M[n][0],x,scelta)-x;
			break;
		}
		case 2:{
			y=(4*M[n-1][1]-M[n-2][1]+2*h*f(M[n][0],x,scelta))/3-x;
			break;
		}
		case 3:{
			y=(18*M[n-1][1]-9*M[n-2][1]+2*M[n-3][1]+6*h*f(M[n][0],x,scelta))/11-x;
			break;
		}
		case 4:{
            y=(48*M[n-1][1]-36*M[n-2][1]+16*M[n-3][1]-3*M[n-4][1]+12*h*f(M[n][0],x,scelta))/25.-x;		
            break;
		}
	}

	

	return y;
}



double derivata(double x, double M[DIM][2], double h, int n, int met, int step, int scelta){
	double der_x;
	der_x=(phi(x+0.0001,M,h,n,met,step,scelta)-phi(x-0.0001,M,h,n,met,step,scelta))*5000;
return der_x; 
}



int esempio(){
	int scelta_funzione; 
	do {
		printf("\n\n\n\nSCELTA DELL'ESEMPIO\n\n1 : y'=t*e^(3t)-2y\n2 : y'=lambda*y\n3 : y'=y*log(y)\n4 : y'=e^(-(t+y))\n5 : y'=ky(m-y)\nSelezionare l'esempio test: ");
		scanf("%d",&scelta_funzione);
	} while (scelta_funzione < 1 && scelta_funzione > 6);


	if(scelta_funzione == 2){
		printf("\nlambda= ");
		scanf("%lf", &lambda);
	}

	else if(scelta_funzione == 5){
		printf("\nk= ");
		scanf("%lf", &k);
		printf("\nm= ");
		scanf("%lf", &m);
	}

	return scelta_funzione;
}



double f(double t, double y, int scelta){
	double res;
	switch(scelta){
		case 1:{
			res = t*exp(3*t) - 2*y;
		break;

		}

		case 2:{
			res = lambda*y;
		break;
		}

		case 3:{
			res = y*log(y);
			break;
		}
		case 4:{
			res = exp(-1*(t+y));
			break;
		}

		case 6:{
			res = y*(1-y);
			break;
		}
		case 5:{
			res = k*y*(m - y);
			break;
		}
	}

	return res;
}



void grafico(double M[DIM][2], int N, int scelta){

    	file=fopen("dati3.txt", "wt");
    	for (int n=0; n<=N; n++){
    		fprintf(file, "%lf\t%lf\n", M[n][0], M[n][1]);
    	}
    	fclose(file);
    	file=fopen("comandi3.txt", "wt");
        fprintf(file, "set title 'soluzione approssimata'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati3.txt\" using 1:2 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi2.txt");
     	//gnuplot Domenico
     	//system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi3.txt");
     	remove("comandi3.txt");
     	remove("dati3.txt");

}



void grafico1(double M[DIM][2], int N, int scelta){

	double sol_esatta_tn;
    //double err=0;

    if(scelta == 1){
    	file=fopen("dati1.txt", "wt");
    	for (int n=0; n<=N; n++){
    		//sol_esatta_tn=((double) 1/25)*exp(-1*2*M[n][0])*(exp(5*M[n][0])*(5*M[n][0]-1)+1);
    		sol_esatta_tn=(M[1][0]+0.04)*exp(2*(M[0][0]-M[n][0]))+0.2*(M[n][0]-M[0][0])*exp(3*(M[n][0]-M[0][0]))-0.04*exp(3*(M[n][0]-M[0][0]));
    		//if (err<fabs(M[n][1]-sol_esatta_tn)) err=sol_esatta_tn;
    		fprintf(file, "%lf\t%lf\t%lf\n", M[n][0], M[n][1], sol_esatta_tn);
    	}
    	fclose(file);
    	//printf("\n\n%lf", err);
    	file=fopen("comandi1.txt", "wt");
        fprintf(file, "set title 'soluzione approssimata'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati1.txt\" using 1:2 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi1.txt");
     	//gnuplot Domenico
        //system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi1.txt");
     	remove("comandi1.txt");
        file=fopen("comandi1.txt", "wt");
        fprintf(file, "set title 'soluzione esatta'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati1.txt\" using 1:3 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi1.txt");
     	//gnuplot Domenico
        //system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi1.txt");
     	remove("comandi1.txt");
     	remove("dati1.txt");
    }

    else if(scelta == 2){
    	file=fopen("dati2.txt", "wt");
    	for (int n=0; n<=N; n++){
    		sol_esatta_tn=M[0][1]*exp(lambda*(M[n][0]-M[0][0]));
    		//if (err<fabs(M[n][1]-sol_esatta_tn)) err=sol_esatta_tn;
    		fprintf(file, "%lf\t%lf\t%lf\n", M[n][0], M[n][1], sol_esatta_tn);
    	}
    	fclose(file);
    	//printf("\n\n%lf", err);
    	file=fopen("comandi2.txt", "wt");
        fprintf(file, "set title 'soluzione approssimata'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati2.txt\" using 1:2 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi2.txt");
     	//gnuplot Domenico
        //system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi2.txt");
     	remove("comandi2.txt");
        file=fopen("comandi2.txt", "wt");
        fprintf(file, "set title 'soluzione esatta'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati2.txt\" using 1:3 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi2.txt");
     	//gnuplot Domenico
        //system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi2.txt");
     	remove("dati2.txt");
     	remove("comandi2.txt");    
     }

    else{
    	file=fopen("dati3.txt", "wt");
    	for (int n=0; n<=N; n++){
    		fprintf(file, "%lf\t%lf\n", M[n][0], M[n][1]);
    	}
    	fclose(file);
    	file=fopen("comandi3.txt", "wt");
        fprintf(file, "set title 'soluzione approssimata'\nset xlabel 't'\nset ylabel 'y'\nset autoscale\nset grid\nplot \"dati3.txt\" using 1:2 with lines\n");
     	fclose(file);
     	//gnuplot Andrea
     	system("gnuplot -p comandi3.txt");
     	//gnuplot Domenico
     	//system("C:/Users/dluon/OneDrive/Desktop/Programmazione/Grafiche/bin/gnuplot.exe -p comandi3.txt");
     	remove("comandi3.txt");
     	remove("dati3.txt");
    }

}


double newton(double x0, double M[DIM][2], double h, int n, int met, int step, int scelta){
	int cont=0;
	double x1, x2;
	x2=x0;
    
    do{
    	x1=x2;
     	x2=x1-phi(x1,M,h,n,met,step,scelta)/derivata(x1,M,h,n,met,step,scelta);
     	
     	cont++;
     } while(fabs(x1-x2)>eps && cont<DIM);


    return x2;
}

double phi(double x, double M[DIM][2], double h, int n, int met, int step, int scelta){
	double y;
	if (scelta == 2) y=admou(x,M,h,n,step,scelta);
	else  y=bdf(x,M,h,n,step,scelta);
	return y;
}


void rk4(double M[DIM][2], double h, int n, int scelta){
	double k1,k2,k3,k4;

	k1 = f(M[n-1][0], M[n-1][1], scelta);
	k2 = f(M[n-1][0]+h/2, M[n-1][1]+h/2*k1, scelta);
	k3 = f(M[n-1][0]+h/2, M[n-1][1]+h/2*k2, scelta);
	k4 = f(M[n-1][0]+h, M[n-1][1]+h*k3, scelta);

	M[n][1] = M[n-1][1] + h/6*(k1+2*k2+2*k3+k4);
}







