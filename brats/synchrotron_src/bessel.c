/* Numerical Recipes Bessel functions */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define NUSE1 5
#define NUSE2 5
#define EPS 1.0e-10
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define PI 3.141592653589793

static void nrerror(char *s) {
  fprintf(stderr,"NR error: %s\n",s);
  exit(1);
}

float chebev(float a, float b, float c[], int m, float x) {

  float d=0.0,dd=0.0,sv,y,y2;
  int j;
  if ((x-a)*(x-b) > 0.0) nrerror("x not in range in routine chebev");
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>=1;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi) {
  float xx;
  static float c1[] = { -1.142022680371168e0,6.5165112670737e-3, 3.087090173086e-4,-3.4706269649e-6,6.9437664e-9, 3.67795e-11,-1.356e-13};
  static float c2[] = { 1.843740587300905e0,-7.68528408447867e-2, 1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8, 2.423096e-10,-1.702e-13,-1.49e-15};
  xx=8.0*x*x-1.0;
  *gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
  *gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
  *gampl= *gam2-x*(*gam1);
  *gammi= *gam2+x*(*gam1);
}


void bessjy(float x, float xnu, float *rj, float *ry, float *rjp, float *ryp) {

  /* Returns the Bessel functions rj = J_nu , ry = Y_nu and their
     derivatives rjp = J0 , ryp = Y0 , for positive x and for xnu =
     nu, where nu > 0. The relative accuracy is within one or two
     signi cant digits of EPS, except near a zero of one of the
     functions, where EPS controls its absolute accuracy. FPMIN is a
     number close to the machine's smallest oating-point number. All
     internal arithmetic is in double precision. */

  int i,isign,l,nl;
  double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2, fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl, rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1, temp,w,x2,xi,xi2,xmu,xmu2;
  
  if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessjy");
  nl=(x < XMIN ? (int)(xnu+0.5) : IMAX(0,(int)(xnu-x+1.5)));
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  w=xi2/PI;
  isign=1;
  h=xnu*xi;
  if (h < FPMIN) h=FPMIN;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=1;i<=MAXIT;i++) { b += xi2;
  d=b-d;
  if (fabs(d) < FPMIN) d=FPMIN;
  c=b-1.0/c;
  if (fabs(c) < FPMIN) c=FPMIN;
  d=1.0/d;
  del=c*d;
  h=del*h;
  if (d < 0.0) isign = -isign;
  if (fabs(del-1.0) < EPS) break;
  }
  if (i > MAXIT) nrerror("x too large in bessjy; try asymptotic expansion");
  rjl=isign*FPMIN;
  rjpl=h*rjl;
  rjl1=rjl;
  rjp1=rjpl;
  fact=xnu*xi;
  for (l=nl;l>=1;l--) { rjtemp=fact*rjl+rjpl;
  fact -= xi;
  rjpl=fact*rjtemp-rjl;
  rjl=rjtemp;
  }
  if (rjl == 0.0) rjl=EPS;
  f=rjpl/rjl;
  if (x < XMIN) {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,&gam1,&gam2,&gampl,&gammi);
    ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
    e=exp(e);
    p=e/(gampl*PI);
    q=1.0/(e*PI*gammi);
    pimu2=0.5*pimu;
    fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
    r=PI*pimu2*fact3*fact3;
    c=1.0;
    d = -x2*x2;
    sum=ff+r*q;
    sum1=p;
    for (i=1;i<=MAXIT;i++) { ff=(i*ff+p+q)/(i*i-xmu2);
    c *= (d/i);
    p /= (i-xmu);
    q /= (i+xmu);
    del=c*(ff+r*q);
    sum += del;
    del1=c*p-i*del;
    sum1 += del1;
    if (fabs(del) < (1.0+fabs(sum))*EPS) break;
    }
    if (i > MAXIT) nrerror("bessy series failed to converge");
    rymu = -sum;
    ry1 = -sum1*xi2;
    rymup=xmu*xi*rymu-ry1;
    rjmu=w/(rymup-f*rymu);
  } else {
    a=0.25-xmu2;
    p = -0.5*xi;
    q=1.0;
    br=2.0*x;
    bi=2.0;
    fact=a*xi/(p*p+q*q);
    cr=br+q*fact;
    ci=bi+p*fact;
    den=br*br+bi*bi;
    dr=br/den;
    di = -bi/den;
    dlr=cr*dr-ci*di;
    dli=cr*di+ci*dr;
    temp=p*dlr-q*dli;
    q=p*dli+q*dlr;
    p=temp;
    for (i=2;i<=MAXIT;i++) {
      a += 2*(i-1);
      bi += 2.0;
      dr=a*dr+br;
      di=a*di+bi;
      if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
      fact=a/(cr*cr+ci*ci);
      cr=br+cr*fact;
      ci=bi-ci*fact;
      if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
      den=dr*dr+di*di;
      dr /= den;
      di /= -den;
      dlr=cr*dr-ci*di;
      dli=cr*di+ci*dr;
      temp=p*dlr-q*dli;
      q=p*dli+q*dlr;
      p=temp;
      if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
    }
    if (i > MAXIT) nrerror("cf2 failed in bessjy");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    rjmu=SIGN(rjmu,rjl);
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }
  fact=rjmu/rjl;
  *rj=rjl1*fact;
  *rjp=rjp1*fact;
  for (i=1;i<=nl;i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;
    ry1=rytemp;
  }
  *ry=rymu;
  *ryp=xnu*xi*rymu-ry1;
}

void bessik(float x, float xnu, float *ri, float *rk, float *rip, float *rkp) {
  
  int i,l,nl;
  double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2, gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl, ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;

  if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessik");
  nl=(int)(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  if (h < FPMIN) h=FPMIN;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=1;i<=MAXIT;i++) {
    b += xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > MAXIT) nrerror("x too large in bessik; try asymptotic expansion");
  ril=FPMIN;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl;l>=1;l--) {
    ritemp=fact*ril+ripl;
    fact -= xi;
    ripl=fact*ritemp+ril;
    ril=ritemp;
  }
  f=ripl/ril;
  if (x < XMIN) {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,&gam1,&gam2,&gampl,&gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for (i=1;i<=MAXIT;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*ff;
      sum += del;
      del1=c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*EPS) break;
    }
    if (i > MAXIT) nrerror("bessk series failed to converge");
    rkmu=sum;
    rk1=sum1*xi2;
  } else {
    b=2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=2;i<=MAXIT;i++) {
      a -= 2*(i-1);
      c = -a*c/i;
      qnew=(q1-b*q2)/a;
      q1=q2;
      q2=qnew;
      q += c*qnew;
      b += 2.0;
      d=1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels=q*delh;
      s += dels;
      if (fabs(dels/s) < EPS) break;
    }
    if (i > MAXIT) nrerror("bessik: failure to converge in cf2");
    h=a1*h;
    rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu-rk1;
  rimu=xi/(f*rkmu-rkmup);
  *ri=(rimu*ril1)/ril;
  *rip=(rimu*rip1)/ril;
  for (i=1;i<=nl;i++) {
    rktemp=(xmu+i)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  *rk=rkmu;
  *rkp=xnu*xi*rkmu-rk1;
}
