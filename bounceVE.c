/**
 * @file bounce.c
 * @author Vatsal Sanjay (vatsalsanjay@gmail.com)
 * vatsalsanjay.com
 * Physics of Fluids
 * @date 2023-05-30
 * 
*/

// 1 is drop
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
// #include "log-conform-ViscoElastic_v6.h" // VE part
#include "tension.h"

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// air-water
#define Rho21 (1e-3)
// Calculations!
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define MINlevel 2                                            // maximum level
#define tsnap (0.01)

// scalar Gpd[], lambdav[]; // VE part

int main(int argc, char const *argv[]) {

  // if (argc < 8){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
  //   return 1;
  // }

  MAXlevel = 10;
  We = 20.0; // We is 1 for 0.22 m/s <1250*0.22^2*0.001/0.06>
  
  tmax = 0.50;

  Ohd = 1e-2; // <\mu/sqrt(1250*0.060*0.001)>
  Ohs = 1e-4; //\mu_r * Ohd
  
  Ldomain = 4.0; // size of domain. must keep Ldomain \gg 1

  // Ec = 0.0; De = 0.0; // VE part
  
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);

  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;

  // // polymers
  // Gp = Gpd; // VE part
  // lambda = lambdav; // VE part

  run();
}

// event properties (i++) { // VE part
//   foreach () {
//     Gpd[] = (f[] < 1e-6) ? 0.0: Ec/We; //clamp(f[],0.,1.)*Ec; //(f[] > 1.-1e-6) ? Ec: 0.0;
//     lambdav[] =  (f[] < 1e-6) ? 0.0: sqrt(We)*De; //clamp(f[],0.,1.)*De; //(f[] > 1.-1e-6) ? De: 0.0;
//   }
// }

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y));
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  // adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
  //   (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
  //   MAXlevel, MAXlevel-4); // VE part

  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MAXlevel-4);
  //unrefine(x>150.0);

}

// Outputs
// static
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0., Vcm = 0., vol = 0.;
  foreach (reduction(+:ke), reduction(+:Vcm)){
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    Vcm += 2*pi*y*((f[])*(u.x[]))*sq(Delta);
    vol += 2*pi*y*(f[])*sq(Delta);
  }
  Vcm /= vol;

  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke p\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);
      fprintf (fp, "i dt t ke Vcm\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
  }

}