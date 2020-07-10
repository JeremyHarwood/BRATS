/*
  Draw the beam on a map
*/


int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa) {

  // PGLOT prototypes
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgsci(int ci);

  bmaj /= 2;
  bmin /= 2;

  // Declare variables
  int a, numpoints;  

  float *xpoints, *ypoints;    /* Outline ellipse in arrays */
  float degtorad = 0.0174532925;

  double cosrot, sinrot, cosangle, sinangle, polarradius, denominator, xrot, yrot;
  double x  = 0.0;
  double y  = 0.0;

  // Allocate the memory
  xpoints = (float *)calloc(362, sizeof(float));
  ypoints = (float *)calloc(362, sizeof(float));

  cosrot  = cos(bpa * degtorad);
  sinrot  = sin(bpa * degtorad);

  numpoints = 0;

  // Loop through 360 degrees
  for (a = 0; a <= 360; a++) {

    cosangle = cos(a * degtorad);
    sinangle = sin(a * degtorad);

    denominator = (bmin*cosangle * bmin*cosangle + bmaj*sinangle * bmaj*sinangle);

    // Get the point in polar co-ordinates
    if (denominator == 0.0) {
      polarradius = 0.0;
    }
    else {
      polarradius = sqrt( bmin*bmaj * bmin*bmaj / denominator );
    }

    x = polarradius * cosangle;
    y = polarradius * sinangle;

    // Rotate the point and put it where requested on the map.
    xrot = x * cosrot - y * sinrot  + xpos;
    yrot = x * sinrot + y * cosrot  + ypos;
    xpoints[a] = (float) xrot;
    ypoints[a] = (float) yrot;
    numpoints++;
  }

  // Set plot colour
  cpgsci(1);

  // Plot the ellipse on the map
  cpgline(numpoints, xpoints, ypoints);
   
  // Free up the allocated memory
  free(xpoints);
  free(ypoints);

  return 0;
}


