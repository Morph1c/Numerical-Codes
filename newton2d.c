#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 500
#define h 0.01
#define eps 0.000001

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

int newton(int f, int f2, double X[MAX_ITER][2]);
int jacobian(double j_valued[2][2], int f, int f2, double x_k, double y_k);
int gauss_2d(double sol[2], double M[2][2], double b[2]);
double fun(int f, double x, double y);
void pivoting(double M[2][2], double b[2]);

int main(void){
	int option;
	int iter_num;
	double zeros_sys[MAX_ITER][2];
	//zeros_sys[0][0] = 0;
	//zeros_sys[0][1] = 0;

	printf("Scegliere il sistema non lineare di cui calcolare la radice:\n[1, 2] y - 2x^2 = 0, x - y^2 + 1 = 0\n[3, 4] e^(x^2+y^2) - 1 = 0, e^(x^2-y^2) - 1 = 0: ");
	scanf("%d", &option);

	printf("Scegli i punti iniziali x0, y0: ");
	scanf("%lf", &zeros_sys[0][0]);
	scanf("%lf", &zeros_sys[0][1]);

	iter_num = newton(option, option + 1, zeros_sys);
	//if(fabs(fun(1, zeros_sys[iter_num][0], zeros_sys[iter_num][1]) + fun(2, zeros_sys[iter_num][0], zeros_sys[iter_num][1])) > eps || zeros_sys[iter_num][0] != zeros_sys[iter_num][0]){
	//	printf("Attenzione il punto iniziale scelto non è valido\n");
	//}
	//else{
		printf("\n%f %f\n", zeros_sys[iter_num][0], zeros_sys[iter_num][1]);
	//}

	return 0;
}

int newton(int f, int f2, double X[MAX_ITER][2]){
	double F_k[2]; // F(X_k)
	int iter = 1;
	double sol_k[2];
	double X_k1[2];
	X[1][0] = X[0][0]; // sarebbe da dove voglio partire
	X[1][1] = X[0][1];
	X[0][0] = X[0][1] = eps + 1; // artificio per non far chiudere subito il ciclo for
	double j[2][2];
	while(fabs(max(X[iter][0] - X[iter - 1][0], X[iter][1] - X[iter - 1][1])) > eps && iter < MAX_ITER){
		
		jacobian(j, f, f2, X[iter][0], X[iter][1]); // dF1/dx, con x = y_n, y = v_n 
		F_k[0] = fun(f, X[iter][0], X[iter][1]);
		F_k[1] = fun(f2, X[iter][0], X[iter][1]);
		pivoting(j, F_k);
		gauss_2d(sol_k, j, F_k);
		//printf("j a %d: [0][0] %lf [0][1] %lf [1][0] %lf [1][1] %lf and F_k[0] = %lf, F_k[1] = %lf and sol is 0: %lf 1: %lf\n", iter, j[0][0], j[0][1], j[1][0], j[1][1], F_k[0], F_k[1], sol_k[0], sol_k[1]);

		X[iter + 1][0] = X[iter][0] - sol_k[0];
		X[iter + 1][1] = X[iter][1] - sol_k[1];
		printf("X_%d : %lf Y_%d : %lf F_1 = %lf, F_2 = %lf\n", iter, X[iter+1][0], iter, X[iter+1][1], fun(f, X[iter+1][0], X[iter+1][1]), fun(f2, X[iter+1][0], X[iter+1][1])); 

		//printf("max %lf\n", max(X[iter + 1][0] - X[iter][0], X[iter + 1][1] - X[iter][1]));
		iter++;

	}

	return iter;
}

void pivoting(double M[2][2], double b[2]){
	double support = 0;
	if(fabs(M[0][1]) > fabs(M[0][0])){
		support = M[0][0];
		M[0][0] = M[1][0];
		M[1][0] = support;
		support = M[0][1];
		M[0][1] = M[1][1];
		M[1][1] = support;
		support = b[0];
		b[0] = b[1];
		b[1] = support;
	}
}

int jacobian(double j_valued[2][2], int f, int f2, double x_k, double y_k){
	j_valued[0][0] = (fun(f, x_k + h, y_k) - fun(f, x_k - h, y_k)) / (2*h);//dF1/dx
	j_valued[0][1] = (fun(f2, x_k + h, y_k) - fun(f2, x_k - h, y_k)) / (2*h);
	j_valued[1][0] = (fun(f, x_k, y_k + h) - fun(f, x_k, y_k - h)) / (2*h);
	j_valued[1][1] = (fun(f2, x_k, y_k + h) - fun(f2, x_k, y_k - h)) / (2*h);

	return 0;

}

int gauss_2d(double sol[2], double M[2][2], double b[2]){
	double k = M[1][0] / M[0][0];
	sol[1] = (b[1] - (k * b[0])) / (M[1][1] - (k*M[0][1]));
	sol[0] = (b[0] - M[0][1]*sol[1]) / M[0][0];

	return 0;
}


double fun(int f, double x, double y){
	float res;
	switch(f){
		case 1: res = y - 2 * (x * x); break;
		case 2: res = x - (y * y) + 1; break;
		case 3: res = exp(x*x+y*y) - 1; break;
		case 4: res = exp(x*x - y*y) - 1; break;
	}

	return res;
}