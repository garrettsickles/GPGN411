// To compile: gcc -std=c99 -o lab02 lab02.c -lm

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
	double B_a[DIM];
	double total;
};


double dist(double a, double b, double c) {
	return sqrt(pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0));
}

double tensor(unsigned int d1, unsigned int d2, struct point a, struct prism p) {
	unsigned int index = 0;
	double swap = 1.0;
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
					// index = (~((1 << d1) + (1 << d2))) & (7);
					num[i][j][k] = dist(d[X][i], d[Y][j], d[Z][k]);

					if ((d1==X && d2==Y) || (d1==X && d2==Y)) num[i][j][k] += d[Z][k];
					if ((d1==X && d2==Z) || (d1==Z && d2==X)) num[i][j][k] += d[Y][j];
					if ((d1==Y && d2==Z) || (d1==Z && d2==Y)) num[i][j][k] += d[X][i];
					// (X == index) ? (num[i][j][k] += d[X][i]) : 1;
					// (Y == index) ? (num[i][j][k] += d[Y][j]) : 1;
					// (Z == index) ? (num[i][j][k] += d[Z][k]) : 1;

					// Only use terms that are nonzero
					if (num[i][j][k] != 0.0) val[i][j][k] = log(num[i][j][k]);
					else {
						val[i][j][k] = 0.0;
						swap = -1.0;
					}
				}
				val[i][j][k] *= (((i + j + k) % 2) ? (-1.0) : (1.0));
				value += val[i][j][k];
			}
		}
	}

	return value*swap;
}

int main() {

	char filename[30];
	
	struct survey s;
	struct prism p;

	double kappa, inc, dec, strength;

	FILE* file;

	printf("Welcome to the rectangular prism tensor response calculator!\n");
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
	//Get Magnetic Stuff
	fscanf(file, "%lf", &(kappa));
	fscanf(file, "%lf", &(strength));
	fscanf(file, "%lf", &(inc));
	fscanf(file, "%lf", &(dec));
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
	printf("\tSusceptibility: %lf\n", kappa);
	printf("\tField Strength: %lf\n", strength);
	printf("\tInclination: %lf\n", inc);
	printf("\tDeclination: %lf\n", dec);
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

	struct observation **grid = (struct observation**)(malloc(n * sizeof(struct observation*)));
	for (int i = 0; i < n; i++) {
		grid[i] = (struct observation*)(malloc(e * sizeof(struct observation)));
	}

	double north = 0.0;
	double east = 0.0;
	struct point pnt;

	double B_n[DIM], B_h[DIM];

	B_h[X] = cos(inc * 3.14159 / 180.0) * cos(dec * 3.14159 / 180.0);
	B_h[Y] = cos(inc * 3.14159 / 180.0) * sin(dec * 3.14159 / 180.0);
	B_h[Z] = sin(inc * 3.14159 / 180.0);

	B_n[X] = B_h[X] * strength;
	B_n[Y] = B_h[Y] * strength;
	B_n[Z] = B_h[Z] * strength;

	for (int i = 0; i < n; i++) {
		north = s.lower_left.loc[X] + (double)(i * s.ivl[X]);
		pnt.loc[X] = north;
		for (int j = 0; j < e; j++) {
			east = s.lower_left.loc[Y] + (double)(j * s.ivl[Y]);
			pnt.loc[Y] = east;
			grid[i][j].p.loc[X] = north;
			grid[i][j].p.loc[Y] = east;
			grid[i][j].tensor_matrix[X][X] = tensor(X, X, pnt, p);
			grid[i][j].tensor_matrix[X][Y] = tensor(X, Y, pnt, p);
			grid[i][j].tensor_matrix[X][Z] = tensor(X, Z, pnt, p);
			grid[i][j].tensor_matrix[Y][X] = grid[i][j].tensor_matrix[X][Y];
			grid[i][j].tensor_matrix[Y][Y] = tensor(Y, Y, pnt, p);
			grid[i][j].tensor_matrix[Y][Z] = tensor(Y, Z, pnt, p);
			grid[i][j].tensor_matrix[Z][X] = grid[i][j].tensor_matrix[X][Z];
			grid[i][j].tensor_matrix[Z][Y] = grid[i][j].tensor_matrix[Y][Z];
			grid[i][j].tensor_matrix[Z][Z] = tensor(Z, Z, pnt, p);

			grid[i][j].B_a[X] = B_n[X] * grid[i][j].tensor_matrix[X][X];
			grid[i][j].B_a[X] += B_n[Y] * grid[i][j].tensor_matrix[X][Y];
			grid[i][j].B_a[X] += B_n[Z] * grid[i][j].tensor_matrix[X][Z];
			grid[i][j].B_a[X] *= (kappa / (4 * 3.14159));

			grid[i][j].B_a[Y] = B_n[X] * grid[i][j].tensor_matrix[Y][X];
			grid[i][j].B_a[Y] += B_n[Y] * grid[i][j].tensor_matrix[Y][Y];
			grid[i][j].B_a[Y] += B_n[Z] * grid[i][j].tensor_matrix[Y][Z];
			grid[i][j].B_a[Y] *= (kappa / (4 * 3.14159));

			grid[i][j].B_a[Z] = B_n[X] * grid[i][j].tensor_matrix[Z][X];
			grid[i][j].B_a[Z] += B_n[Y] * grid[i][j].tensor_matrix[Z][Y];
			grid[i][j].B_a[Z] += B_n[Z] * grid[i][j].tensor_matrix[Z][Z];
			grid[i][j].B_a[Z] *= (kappa / (4 * 3.14159));

			grid[i][j].total = (B_h[X] * grid[i][j].B_a[X]) + (B_h[Y] * grid[i][j].B_a[Y]) + (B_h[Z] * grid[i][j].B_a[Z]);
		}
	}

	FILE* bx = fopen("B_x_m.csv", "w");
	FILE* by = fopen("B_y_m.csv", "w");
	FILE* bz = fopen("B_z_m.csv", "w");
	FILE* tt = fopen("Total_m.csv", "w");

	if (bx == NULL || by == NULL || bz == NULL || tt == NULL) {
		printf("Error opening output files!\n");
		exit(1);
	}

	// Matrix Form
	for (int i = n-1; i >= 0; i--) {
		for (int j = e-1; j >= 0; j--) {
			fprintf(bx, "%lf", grid[i][j].B_a[X]);
			fprintf(bx, ((j == 0) ? "\n" : ","));
			fprintf(by, "%lf", grid[i][j].B_a[Y]);
			fprintf(by, ((j == 0) ? "\n" : ","));
			fprintf(bz, "%lf", grid[i][j].B_a[Z]);
			fprintf(bz, ((j == 0) ? "\n" : ","));
			fprintf(tt, "%lf", grid[i][j].total);
			fprintf(tt, ((j == 0) ? "\n" : ","));
		}
	}

	fclose(bx);
	fclose(by);
	fclose(bz);
	fclose(tt);

	bx = fopen("B_x_p.csv", "w");
	by = fopen("B_y_p.csv", "w");
	bz = fopen("B_z_p.csv", "w");
	tt = fopen("Total_p.csv", "w");

	// Point Form
	for (int i = n-1; i >= 0; i--) {
		for (int j = e-1; j >= 0; j--) {
			fprintf(bx, "%lf,%lf,", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(bx, "%lf\n", grid[i][j].B_a[X]);

			fprintf(by, "%lf,%lf,", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(by, "%lf\n", grid[i][j].B_a[Y]);

			fprintf(bz, "%lf,%lf,", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(bz, "%lf\n", grid[i][j].B_a[Z]);

			fprintf(tt, "%lf,%lf,", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(tt, "%lf\n", grid[i][j].total);
		}
	}

	bx = fopen("B_x_p.dat", "w");
	by = fopen("B_y_p.dat", "w");
	bz = fopen("B_z_p.dat", "w");
	tt = fopen("Total_p.dat", "w");

	fprintf(bx, "%d\t%d\n", e, n);
	fprintf(by, "%d\t%d\n", e, n);
	fprintf(bz, "%d\t%d\n", e, n);
	fprintf(tt, "%d\t%d\n", e, n);

	// Point Form
	for (int j = 0; j < e; j++) {
		for (int i = 0; i < n; i++) {
			fprintf(bx, "%lf\t%lf\t", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(bx, "%lf\n", grid[i][j].B_a[X]);

			fprintf(by, "%lf\t%lf\t", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(by, "%lf\n", grid[i][j].B_a[Y]);

			fprintf(bz, "%lf\t%lf\t", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(bz, "%lf\n", grid[i][j].B_a[Z]);

			fprintf(tt, "%lf\t%lf\t", grid[i][j].p.loc[Y], grid[i][j].p.loc[X]);
			fprintf(tt, "%lf\n", grid[i][j].total);
		}
	}

	FILE* x_obs = fopen("x_obs.csv", "w");
	FILE* y_obs = fopen("y_obs.csv", "w");

	for (int i = 0; i < n; i++) fprintf(x_obs, "%lf\n", grid[i][0].p.loc[X]);
	for (int j = 0; j < e; j++) fprintf(y_obs, "%lf\n", grid[0][j].p.loc[Y]);

	fclose(x_obs);
	fclose(y_obs);

	fclose(bx);
	fclose(by);
	fclose(bz);
	fclose(tt);

	for (int i = 0; i < n; i++) free(grid[i]);
	free(grid);

	printf("Done\n");
	return 0;
}
