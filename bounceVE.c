/**
 * @file bounceVE_WeissenbergSweep.c
 * @brief Simulation of a viscoelastic drop using Basilisk C.
 *  This code simulates the dynamics of a viscoelastic drop in an axisymmetric domain. 
 * The simulation parameters include various dimensionless numbers such as the Elasto-capillary number (Ec), Ohnesorge number (Oh), 
 * and Deborah number (De). The code also handles two-phase flow with surface tension and viscoelastic effects.
 * @author Vatsal Sanjay (vatsalsanjay@gmail.com)
 * vatsalsanjay.com
 * Physics of Fluids
 * @date 2024-06-20
 * 
*/

// Include necessary Basilisk C headers
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "log-conform-ViscoElastic_v6.h"
#include "tension.h"

// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// Density ratio between two phases (air-water)
#define Rho21 (1e-3)

// Distance and radius calculations
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// Boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

// Global variables for simulation parameters
int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define MINlevel 2  // Minimum refinement level
#define tsnap (0.01)  // Time interval for snapshots

// Scalars for viscoelastic properties
scalar Gpd[], lambdav[];

int main(int argc, char const *argv[]) {

  // if (argc < 8){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
  //   return 1;
  // }

  MAXlevel = 10;
  We = 100.0; // We is 1 for 0.22 m/s <1250*0.22^2*0.001/0.06>
  tmax = 1.00;
  Ohd = 1e-2;  // Ohnesorge number for drop
  Ohs = 1e-4;  // Ohnesorge number for surrounding fluid
  Ldomain = 4.0;  // Size of the computational domain
  Ec = 0.0; // Elasto-capillary number
  De = 0.0; // Deborah number

  // for Newtonian fluids, to be on the safe side, we should keep both El and Wi as 0–although either being 0 gives Newtonian response. 
  
  // Print simulation parameters to standard error
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);

  // Initialize the computational domain
  L0 = Ldomain;
  X0 = 0.;
  Y0 = 0.;
  init_grid(1 << 4);

  // Create directory for intermediate files
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // Set fluid properties
  rho1 = 1.0;
  mu1 = Ohd / sqrt(We);
  rho2 = Rho21;
  mu2 = Ohs / sqrt(We);
  f.sigma = 1.0 / We;

  // Set viscoelastic properties
  Gp = Gpd;
  lambda = lambdav;

  // Start the simulation
  run();
}

// Event to set viscoelastic properties
event properties (i++) {
  foreach () {
    Gpd[] = (f[] < 1e-6) ? 0.0: Ec/We;
    lambdav[] =  (f[] < 1e-6) ? 0.0: sqrt(We)*De;
  }
}

// Event to initialize the simulation
event init(t = 0) {
  if (!restore(file = "dump")) {
    refine((R2Drop(x, y) < 1.05) && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y));
    foreach() {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
  }
}

// Scalar for curvature calculation
scalar KAPPA[];

// Event for adaptive mesh refinement
event adapt(i++) {
  curvature(f, KAPPA);
  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
                (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
                MAXlevel, MAXlevel - 4);
}

// Event to write output files at regular intervals
event writingFiles(t = 0, t += tsnap; t <= tmax) {
  dump(file = "dump");
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}

// Event to log simulation data
event logWriting(i++) {
  double ke = 0., Vcm = 0., vol = 0.;
  foreach(reduction(+:ke), reduction(+:Vcm)) {
    ke += 2 * pi * y * (0.5 * rho(f[]) * (sq(u.x[]) + sq(u.y[]))) * sq(Delta);
    Vcm += 2 * pi * y * ((f[]) * (u.x[])) * sq(Delta);
    vol += 2 * pi * y * (f[]) * sq(Delta);
  }
  Vcm /= vol;

  static FILE *fp;

  if (pid() == 0) {
    if (i == 0) {
      fprintf(ferr, "i dt t ke p\n");
      fp = fopen("log", "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", 
              MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);
      fprintf(fp, "i dt t ke Vcm\n");
      fprintf(fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    } else {
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
  }
}
