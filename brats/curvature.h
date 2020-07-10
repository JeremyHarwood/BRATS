/* Find the polynomial regression of a data set.

   What needs to be passed:

   int n - The number of data points (e.g. maps)
   double x[] - An array containing the x data
   double y[] - An array containing the y data
   int poly - The degree of polynomial regression to be calculated
   double coeff[] - An empty array to hold the return values. This array must be atleast the polynomial order +1 in size.

   The data should already be converted to log values if being used for curvature purposes.

   If you have any question please email jeremy.harwood@physics.org and I will try and answer them as soon as I can
*/


int curvature(int n, float *x, float *y, int poly, float *coeff, int printpoly) {

  // Declaring variables

  int a, b, c, d, i, count;
  double ratio, leadone;
  double lhsmatrix[(poly+1)][(poly+1)];
  double rhsmatrix[(poly+1)];

  // Check the order is sensible

  if (poly < 0) {
    fprintf(stderr,"\nError: The polynomial order must be positive!\nReturning to main program...\n\n");
    return 100;
  }

  // Zero out the matrix

for (a = 0; a < (poly+1); a++) {
    for (b = 0; b < (poly+1); b++) {
      lhsmatrix[a][b] = 0;
    }
    rhsmatrix[a] = 0;
 }

  // Populate the matrix

  for (a = 0; a < (poly+1); a++) {
    for (c = 0; c < n; c++) {
      rhsmatrix[a] += (y[c] * pow(x[c], a));
	}
    for  (b = 0; b < (poly+1); b++) {
      if ((a == 0) && (b == 0)) {
	lhsmatrix[a][b] = n;
      }
      else {
	for (d = 0; d < n; d++) {
	  lhsmatrix[a][b] += pow(x[d], (b+a));
	}
      }
    }
  }

  // Do the gaussian elimination

  // Make the triangular matrix

  count = 0;

  for (b = 0; b < (poly+1); b++) {
    count++;

    for (a = count; a < (poly+1); a++) {
      ratio = (lhsmatrix[a][b] / lhsmatrix[b][b]);
      lhsmatrix[a][b] = 0;
      rhsmatrix[a] -= (ratio*rhsmatrix[b]);

      for (i = b+1; i < (poly+1); i++) {
	lhsmatrix[a][i] -= (ratio*lhsmatrix[b][i]);
      }
    }
  }

  count = 0;

  for (a = 0; a < (poly+1); a++) {
    for (b = count; b < (poly+1); b++) {
      if (b == count) {
	leadone = lhsmatrix[a][count];
      }
      lhsmatrix[a][b] /= leadone;
    }

    rhsmatrix[a] /= leadone;
    count++;
  }

  // Work out the coefficients
  
  count = 0;

  coeff[poly] = rhsmatrix[poly];

  for (a = (poly-1); a >= 0; a--) {
    count++;
    coeff[a] = rhsmatrix[a];

    for (b = poly; b > (poly-count); b--) {
      coeff[a] -= (coeff[b]*lhsmatrix[a][b]);
    }
  }

  // Output the results
  if (printpoly == 1) {

    if (poly == 2) {
      printf("Quadratic best fit: y = ");
    }

    if (poly == 3) {
      printf("Cubic best fit: y = ");
    }

    if (poly > 3) {
      printf("%dth order polynomial best fit: y = ", poly);
    }
  
  }

 for (a = (poly); a >= 0; a--) {

   if ((coeff[a] < 1e-32) && (coeff[a] > -1e-32)) {
     coeff[a] = 0.0;
   }

   if (printpoly == 1) {

     if (((coeff[a] < 1000.0) && (coeff[a] > 1e-3)) || ((coeff[a] > -1000.0) && (coeff[a] < -1e-3))) {
     
       if (a == 0) {
	 printf(" + %.3f", coeff[a]);
       }
       else if (a == 1) {
	 printf(" + %.3f x", coeff[a]);
       }
       else if (a == (poly)){
	 printf("%.3f x^%d", coeff[a], (a));
       }
       else {
	 printf(" + %.3f x^%d", coeff[a], (a));
       }
     }

     else {

       if (a == 0) {
	 printf(" + %.3e", coeff[a]);
       }
       else if (a == 1) {
	 printf(" + %.3e x", coeff[a]);
       }
       else if (a == (poly)){
	 printf("%.3e x^%d", coeff[a], (a));
       }
       else {
	 printf(" + %.3e x^%d", coeff[a], (a));
       }
     }
   }
 }
  
 if (printpoly == 1) {
   printf("\n"); 
 }


 return 0;

}
