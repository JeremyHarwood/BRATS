/* 
   Define the global settings for making smooth maps
*/

// Declare the brightness, contrast, and ramp levels (these stay constant)
const float bright=0.5;
const float contra=1.0;
//const float rl[]={-0.5, 0.0, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 1.0};
const float rl[]={0.0, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 1.0};

// Allocate the lrgb memory
//float *rl = (float *)calloc(11, sizeof(float));
float *rr = (float *)calloc(10, sizeof(float));
float *rg = (float *)calloc(10, sizeof(float));
float *rb = (float *)calloc(10, sizeof(float));

//float r1, r2, r3, r4, r5, r6, r7, r8, r9, g1, g2, g3, g4, g5, g6, g7, g8, g9, b1, b2, b3, b4, b5, b6, b7, b8, b9;

float *tr = (float *)calloc(9, sizeof(float));
float *tg = (float *)calloc(9, sizeof(float));
float *tb = (float *)calloc(9, sizeof(float));


// Select the colour scheme
if (pltcolour == 0) { // Red to blue HSV aka rainbow (Default)
  // Red
  tr[0]=1.0;
  tr[1]=1.0;
  tr[2]=1.0;
  tr[3]=0.5;
  tr[4]=0.0;
  tr[5]=0.0;
  tr[6]=0.0;
  tr[7]=0.0;
  tr[8]=0.0;
  // Green
  tg[0]=0.0;
  tg[1]=0.5;
  tg[2]=1.0;
  tg[3]=1.0;
  tg[4]=1.0;
  tg[5]=1.0;
  tg[6]=1.0;
  tg[7]=0.5;
  tg[8]=0.0;
  // Blue
  tb[0]=0.0;
  tb[1]=0.0;
  tb[2]=0.0;
  tb[3]=0.0;
  tb[4]=0.0;
  tb[5]=0.5;
  tb[6]=1.0;
  tb[7]=1.0;
  tb[8]=1.0;
 }
 else if (pltcolour == 1) { // Blue to red RGB (Heat)
  // Red
  tr[0]=0.0;
  tr[1]=0.12;
  tr[2]=0.25;
  tr[3]=0.37;
  tr[4]=0.50;
  tr[5]=0.62;
  tr[6]=0.75;
  tr[7]=0.87;
  tr[8]=1.0;
  // Green
  tg[0]=0.0;
  tg[1]=0.0;
  tg[2]=0.0;
  tg[3]=0.0;
  tg[4]=0.0;
  tg[5]=0.0;
  tg[6]=0.0;
  tg[7]=0.0;
  tg[8]=0.0;
  // Blue
  tb[0]=1.0;
  tb[1]=0.87;
  tb[2]=0.75;
  tb[3]=0.62;
  tb[4]=0.50;
  tb[5]=0.37;
  tb[6]=0.25;
  tb[7]=0.12;
  tb[8]=0.0;
 }
 else if (pltcolour == 2) { // Yellow to Purple HSV (Viridis)
  // Red
  tr[0]=0.87;
  tr[1]=0.48;
  tr[2]=0.12;
  tr[3]=0.37;
  tr[4]=0.0;
  tr[5]=0.0;
  tr[6]=0.0;
  tr[7]=0.0;
  tr[8]=0.20;
  // Green
  tg[0]=0.88;
  tg[1]=0.84;
  tg[2]=0.79;
  tg[3]=0.75;
  tg[4]=0.71;
  tg[5]=0.60;
  tg[6]=0.29;
  tg[7]=0.03;
  tg[8]=0.0;
  // Blue
  tb[0]=0.0;
  tb[1]=0.0;
  tb[2]=0.0;
  tb[3]=0.20;
  tb[4]=0.48;
  tb[5]=0.67;
  tb[6]=0.62;
  tb[7]=0.58;
  tb[8]=0.54;
 }
 else if (pltcolour == 3) { // Red 
  // Red
  tr[0]=1.0;
  tr[1]=1.0;
  tr[2]=1.0;
  tr[3]=1.0;
  tr[4]=1.0;
  tr[5]=1.0;
  tr[6]=1.0;
  tr[7]=1.0;
  tr[8]=1.0;
  // Green
  tg[0]=0.0;
  tg[1]=0.10;
  tg[2]=0.20;
  tg[3]=0.30;
  tg[4]=0.40;
  tg[5]=0.50;
  tg[6]=0.60;
  tg[7]=0.70;
  tg[8]=0.80;
  // Blue
  tb[0]=0.0;
  tb[1]=0.10;
  tb[2]=0.20;
  tb[3]=0.30;
  tb[4]=0.40;
  tb[5]=0.50;
  tb[6]=0.60;
  tb[7]=0.70;
  tb[8]=0.80;
 }
else if (pltcolour == 4) { // Green 
  // Red
  tr[0]=0.0;
  tr[1]=0.10;
  tr[2]=0.20;
  tr[3]=0.30;
  tr[4]=0.40;
  tr[5]=0.50;
  tr[6]=0.60;
  tr[7]=0.70;
  tr[8]=0.80;
  // Green
  tg[0]=1.0;
  tg[1]=1.0;
  tg[2]=1.0;
  tg[3]=1.0;
  tg[4]=1.0;
  tg[5]=1.0;
  tg[6]=1.0;
  tg[7]=1.0;
  tg[8]=1.0;
  // Blue
  tb[0]=0.0;
  tb[1]=0.10;
  tb[2]=0.20;
  tb[3]=0.30;
  tb[4]=0.40;
  tb[5]=0.50;
  tb[6]=0.60;
  tb[7]=0.70;
  tb[8]=0.80;
 }
else if (pltcolour == 5) { // Blue 
  // Red
  tr[0]=0.0;
  tr[1]=0.10;
  tr[2]=0.20;
  tr[3]=0.30;
  tr[4]=0.40;
  tr[5]=0.50;
  tr[6]=0.60;
  tr[7]=0.70;
  tr[8]=0.80;
  // Green
  tg[0]=0.0;
  tg[1]=0.10;
  tg[2]=0.20;
  tg[3]=0.30;
  tg[4]=0.40;
  tg[5]=0.50;
  tg[6]=0.60;
  tg[7]=0.70;
  tg[8]=0.80;
  // Blue
  tb[0]=1.0;
  tb[1]=1.0;
  tb[2]=1.0;
  tb[3]=1.0;
  tb[4]=1.0;
  tb[5]=1.0;
  tb[6]=1.0;
  tb[7]=1.0;
  tb[8]=1.0;
 }
 else if (pltcolour == 6) { // Grey scale 
  // Red
  tr[0]=0.0;
  tr[1]=0.12;
  tr[2]=0.25;
  tr[3]=0.37;
  tr[4]=0.50;
  tr[5]=0.62;
  tr[6]=0.75;
  tr[7]=0.87;
  tr[8]=1.0;
  // Green
  tg[0]=0.0;
  tg[1]=0.12;
  tg[2]=0.25;
  tg[3]=0.37;
  tg[4]=0.50;
  tg[5]=0.62;
  tg[6]=0.75;
  tg[7]=0.87;
  tg[8]=1.0;
  // Blue
  tb[0]=0.0;
  tb[1]=0.12;
  tb[2]=0.25;
  tb[3]=0.37;
  tb[4]=0.50;
  tb[5]=0.62;
  tb[6]=0.75;
  tb[7]=0.87;
  tb[8]=1.0;
 }
else { // Default to rainbow
  // Red
  tr[0]=1.0;
  tr[1]=1.0;
  tr[2]=1.0;
  tr[3]=0.5;
  tr[4]=0.0;
  tr[5]=0.0;
  tr[6]=0.0;
  tr[7]=0.0;
  tr[8]=0.0;
  // Green
  tg[0]=0.0;
  tg[1]=0.5;
  tg[2]=1.0;
  tg[3]=1.0;
  tg[4]=1.0;
  tg[5]=1.0;
  tg[6]=1.0;
  tg[7]=0.5;
  tg[8]=0.0;
  // Blue
  tb[0]=0.0;
  tb[1]=0.0;
  tb[2]=0.0;
  tb[3]=0.0;
  tb[4]=0.0;
  tb[5]=0.5;
  tb[6]=1.0;
  tb[7]=1.0;
  tb[8]=1.0;
 }

// If the colours are inverted, flip the array values
if (invertcolours == 1) {

  float tmp_col1, tmp_col2;
  int i;

  for (i=0; i<4; i++) {
    tmp_col1 = tr[i];
    tmp_col2 = tr[8-i];
    tr[i] = tmp_col2;
    tr[8-i] = tmp_col1;

    tmp_col1 = tg[i];
    tmp_col2 = tg[8-i];
    tg[i] = tmp_col2;
    tg[8-i] = tmp_col1;

    tmp_col1 = tb[i];
    tmp_col2 = tb[8-i];
    tb[i] = tmp_col2;
    tb[8-i] = tmp_col1;
  }
  
 }


// Set the rbg values based on the selected colour scheme
// Red
rr[0]=tr[0];
rr[1]=tr[1];
rr[2]=tr[2];
rr[3]=tr[3];
rr[4]=tr[4];
rr[5]=tr[5];
rr[6]=tr[6];
rr[7]=tr[7];
rr[8]=tr[8];
rr[9]=1.0;
// Green
rg[0]=tg[0];
rg[1]=tg[1];
rg[2]=tg[2];
rg[3]=tg[3];
rg[4]=tg[4];
rg[5]=tg[5];
rg[6]=tg[6];
rg[7]=tg[7];
rg[8]=tg[8];
rg[9]=1.0;
// Blue
rb[0]=tb[0];
rb[1]=tb[1];
rb[2]=tb[2];
rb[3]=tb[3];
rb[4]=tb[4];
rb[5]=tb[5];
rb[6]=tb[6];
rb[7]=tb[7];
rb[8]=tb[8];
rb[9]=1.0;


