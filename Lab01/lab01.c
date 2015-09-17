// To compile: gcc -std=c99 -o lab01 lab01.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X 0
#define Y 1
#define Z 2
#define DIM 3

#define DELTA(X,Y) (X-Y)

struct point {
	double loc[DIM];
};

struct prism {
	struct point center;
	double length[DIM];
};

struct survey {
	struct point lower_left;
	struct point upper_right;
	double ivl[DIM];
};

struct observation {
	struct point p;
	double tensor_matrix[DIM][DIM];
};


double dist(double a, double b, double c) {
	return sqrt(pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0));
}

double tensor(unsigned int d1, unsigned int d2, struct point a, struct prism p) {
	unsigned int index = 0;
	double d[3][2];
	double num[2][2][2];
	double den[2][2][2];
	double val[2][2][2];
	double value = 0.0;

	// Calculate the delta components for all bounds at point a
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < 2; j++) {
			d[i][j] = DELTA(a.loc[i], p.center.loc[i] + p.length[i] * 0.5 * ((j == 0) ? (-1.0) : (1.0)));
		}
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				if (d1 == d2) {
					num[i][j][k] = 1.0;
					den[i][j][k] = dist(d[X][i], d[Y][j], d[Z][k]);
					
					(X == d1) ? (den[i][j][k] *= d[X][i]) : (num[i][j][k] *= d[X][i]);
					(Y == d1) ? (den[i][j][k] *= d[Y][j]) : (num[i][j][k] *= d[Y][j]);
					(Z == d1) ? (den[i][j][k] *= d[Z][k]) : (num[i][j][k] *= d[Z][k]);
					
					val[i][j][k] = atan2(num[i][j][k], den[i][j][k]); 	
				} else {
					index = (~((1 << d1) + (1 << d2))) & (7);
					num[i][j][k] = dist(d[X][i], d[Y][j], d[Z][k]);

					(X == index) ? (num[i][j][k] += d[X][i]) : 1;
					(Y == index) ? (num[i][j][k] += d[Y][j]) : 1;
					(Z == index) ? (num[i][j][k] += d[Z][k]) : 1;

					val[i][j][k] = log(num[i][j][k]);				
				}
				val[i][j][k] *= (((i + j + k) % 2) ? (-1.0) : (1.0));
				value += val[i][j][k];
			}
		}
	}

	return value;
}

int main() {

	char filename[30];
	
	struct survey s;
	struct prism p;
	double density;

	FILE* file;

	printf("Welcome to the rectangular prism gravity gradiometry response calculator!\n");
	printf("Please input the parameter filepath: ");

	scanf("%s", filename);
	file = fopen(filename, "r");

	if (file == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}

	// Get center point
	fscanf(file, "%lf", &(p.center.loc[X]));
	fscanf(file, "%lf", &(p.center.loc[Y]));
	fscanf(file, "%lf", &(p.center.loc[Z]));
	// Get dimensions
	fscanf(file, "%lf", &(p.length[X]));
	fscanf(file, "%lf", &(p.length[Y]));
	fscanf(file, "%lf", &(p.length[Z]));
	// Get Density
	fscanf(file, "%lf", &(density));
	// Get max bounds
	fscanf(file, "%lf", &(s.lower_left.loc[X]));
	fscanf(file, "%lf", &(s.lower_left.loc[Y]));
	s.lower_left.loc[Z] = 0.0;
	// Get min bounds
	fscanf(file, "%lf", &(s.upper_right.loc[X]));
	fscanf(file, "%lf", &(s.upper_right.loc[Y]));
	s.upper_right.loc[Z] = 0.0;
	// Get interval
	fscanf(file, "%lf", &(s.ivl[X]));
	fscanf(file, "%lf", &(s.ivl[Y]));
	s.ivl[Z] = 0.0;
	fclose(file);

	// Show them what's up, ya know
	printf("Here are the values I got...\n");
	printf("Prism Center\n");
	printf("\tX: %lf\n", p.center.loc[X]);
	printf("\tY: %lf\n", p.center.loc[Y]);
	printf("\tZ: %lf\n", p.center.loc[Z]);
	printf("Prism Dimensions\n");
	printf("\tX: %lf\n", p.length[X]);
	printf("\tY: %lf\n", p.length[Y]);
	printf("\tZ: %lf\n", p.length[Z]);
	printf("\tDensity Contrast (g/cc): %lf\n", density);
	printf("Survey Bounds\n");
	printf("\tX Min (Bottom Right): %lf\n", s.lower_left.loc[X]);
	printf("\tY Min (Bottom Left) : %lf\n", s.lower_left.loc[Y]);
	printf("\tX Max (Upper Right) : %lf\n", s.upper_right.loc[X]);
	printf("\tY Max (Upper Left)  : %lf\n", s.upper_right.loc[Y]);
	printf("Survey Spacing\n");
	printf("\tdx : %lf\n", s.ivl[X]);
	printf("\tdy : %lf\n", s.ivl[Y]);
	printf("\n");

	// And here... we... go!

	int n = ((s.upper_right.loc[X] - s.lower_left.loc[X]) / s.ivl[X]) + 1;
	int e = ((s.upper_right.loc[Y] - s.lower_left.loc[Y]) / s.ivl[Y]) + 1;

	printf("%d, %d\n", n, e);

	struct observation **grid = (struct observation**)(malloc(n * sizeof(struct observation*)));
	for (int i = 0; i < n; i++) {
		grid[i] = (struct observation*)(malloc(e * sizeof(struct observation)));
	}

	double north = 0.0;
	double east = 0.0;
	double big_g = 0.667384;
	double rho = density * 1000.0; // From gram/cc to kg/m^3
	struct point pnt;

	for (int i = 0; i < n; i++) {
		north = s.lower_left.loc[X] + (double)(i * s.ivl[X]);
		pnt.loc[X] = north;
		for (int j = 0; j < e; j++) {
			east = s.lower_left.loc[Y] + (double)(j * s.ivl[Y]);
			pnt.loc[Y] = east;
			grid[i][j].p.loc[X] = north;
			grid[i][j].p.loc[Y] = east;
			grid[i][j].tensor_matrix[X][X] = tensor(X, X, pnt, p) * rho * big_g;
			grid[i][j].tensor_matrix[X][Y] = tensor(X, Y, pnt, p) * rho * big_g;
			grid[i][j].tensor_matrix[X][Z] = tensor(X, Z, pnt, p) * rho * big_g;
			grid[i][j].tensor_matrix[Y][X] = (-1.0) * grid[i][j].tensor_matrix[X][Y];
			grid[i][j].tensor_matrix[Y][Y] = tensor(Y, Y, pnt, p) * rho * big_g;
			grid[i][j].tensor_matrix[Y][Z] = tensor(Y, Z, pnt, p) * rho * big_g;
			grid[i][j].tensor_matrix[Z][X] = (-1.0) * grid[i][j].tensor_matrix[X][Z];
			grid[i][j].tensor_matrix[Z][Y] = (-1.0) * grid[i][j].tensor_matrix[Y][Z];
			grid[i][j].tensor_matrix[Z][Z] = tensor(Z, Z, pnt, p) * rho * big_g;
		}
	}

	FILE* txx = fopen("T_xx_m.csv", "w");
	FILE* txy = fopen("T_xy_m.csv", "w");
	FILE* txz = fopen("T_xz_m.csv", "w");
	FILE* tyy = fopen("T_yy_m.csv", "w");
	FILE* tyz = fopen("T_yz_m.csv", "w");
	FILE* tzz = fopen("T_zz_m.csv", "w");

	if (txx == NULL || txy == NULL || txz == NULL || tyy == NULL || tyz == NULL || tzz == NULL) {
		printf("Error opening output files!\n");
		exit(1);
	}

	// Matrix Form
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < e; j++) {
			fprintf(txx, "%lf", grid[i][j].tensor_matrix[X][X]);
			fprintf(txx, ((j+1 == e) ? "\n" : ","));
			fprintf(txy, "%lf", grid[i][j].tensor_matrix[X][Y]);
			fprintf(txy, ((j+1 == e) ? "\n" : ","));
			fprintf(txz, "%lf", grid[i][j].tensor_matrix[X][Z]);
			fprintf(txz, ((j+1 == e) ? "\n" : ","));
			fprintf(tyy, "%lf", grid[i][j].tensor_matrix[Y][Y]);
			fprintf(tyy, ((j+1 == e) ? "\n" : ","));
			fprintf(tyz, "%lf", grid[i][j].tensor_matrix[Y][Z]);
			fprintf(tyz, ((j+1 == e) ? "\n" : ","));
			fprintf(tzz, "%lf", grid[i][j].tensor_matrix[Z][Z]);
			fprintf(tzz, ((j+1 == e) ? "\n" : ","));
		}
	}

	fclose(txx);
	fclose(txy);
	fclose(txz);
	fclose(tyy);
	fclose(tyz);
	fclose(tzz);

	for (int i = 0; i < n; i++) {
		free(grid[i]);
	}
	free(grid);

	printf("Done\n");
	return 0;
}