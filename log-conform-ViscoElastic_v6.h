/** Title: log-conform-elastic_v5.h
# Author: Vatsal Sanjay & Ayush Dixit
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Aug 18, 2023
*/

// The code is same as http://basilisk.fr/src/log-conform.h but written for purely elastic limit (lambda \to \infty)

// In this code, conform_p, conform_qq are in fact the Conformation tensor. 

(const) scalar Gp = unity;
(const) scalar lambda = unity;

#include "bcg.h"

symmetric tensor conform_p[], tau_p[];
#if AXI
scalar conform_qq[], tau_qq[];
#endif

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  foreach() {
    foreach_dimension(){
      tau_p.x.x[] = 0.;
      conform_p.x.x[] = 1.;
    }
    tau_p.x.y[] = 0.;
    conform_p.x.y[] = 0.;
#if AXI
    tau_qq[] = 0;
    conform_qq[] = 1.;
#endif
  }

  for (scalar s in {tau_p}) {
    s.v.x.i = -1; // just a scalar, not the component of a vector
    foreach_dimension(){
      if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
      }
    }
  }

  for (scalar s in {conform_p}) {
    s.v.x.i = -1; // just a scalar, not the component of a vector
    foreach_dimension(){
      if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
      }
    }
  }

#if AXI
  scalar s1 = tau_p.x.y;
  s1[bottom] = dirichlet (0.);  
#endif 

#if AXI
  scalar s2 = conform_p.x.y;
  s2[bottom] = dirichlet (0.);  
#endif 
}

/**
## Numerical Scheme

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  /**
  The eigenvalues are saved in vector $\Lambda$ computed from the
  trace and the determinant of the symmetric conformation tensor
  $\mathbf{A}$. */

  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace of the tensor
  double D = A->x.x*A->y.y - sq(A->x.y); // Determinant

  /**
  The eigenvectors, $\mathbf{v}_i$ are saved by columns in tensor
  $\mathbf{R} = (\mathbf{v}_1|\mathbf{v}_2)$. */

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}

/**
The stress tensor depends on previous instants and has to be
integrated in time. In the log-conformation scheme the advection of
the stress tensor is circumvented, instead the conformation tensor,
$\mathbf{A}$ (or more precisely the related variable $\Psi$) is
advanced in time.

In what follows we will adopt a scheme similar to that of [Hao \& Pan
(2007)](#hao2007). We use a split scheme, solving successively

a) the upper convective term:
$$
\partial_t \Psi = 2 \mathbf{B} + (\Omega \cdot \Psi -\Psi \cdot \Omega)
$$
b) the advection term:
$$
\partial_t \Psi + \nabla \cdot (\Psi \mathbf{u}) = 0
$$
c) the model term (but set in terms of the conformation 
tensor $\mathbf{A}$). In an Oldroyd-B viscoelastic fluid, the model is
$$ 
\partial_t \mathbf{A} = -\frac{\mathbf{f}_r (\mathbf{A})}{\lambda}
$$

The implementation below assumes that the values of $\Psi$ and
$\conform_p$ are never needed simultaneously. This means that $\conform_p$ can
be used to store (temporarily) the values of $\Psi$ (i.e. $\Psi$ is
just an alias for $\conform_p$). */

event tracer_advection(i++)
{
    tensor Psi = conform_p;
#if AXI
    scalar Psiqq = conform_qq;
#endif

    /**
    ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

    foreach() {
      /**
        We assume that the stress tensor $\mathbf{\tau}_p$ depends on the
        conformation tensor $\mathbf{A}$ as follows
        $$
        \mathbf{\tau}_p = G_p (\mathbf{A}) =
        G_p (\mathbf{A} - I)
        $$
      */

      double fa = (f[] < 1e-6 ? 0.0: 1.0);

      pseudo_t A;
      A.x.y = fa*conform_p.x.y[];

      foreach_dimension()
        A.x.x = (fa != 0 ? fa*conform_p.x.x[]: 1.);
      /**
       In the axisymmetric case, $\Psi_{\theta \theta} = \log A_{\theta
       \theta}$. Therefore $\Psi_{\theta \theta} = \log [ ( 1 + \text{fa}
       \tau_{p_{\theta \theta}})]$. */

#if AXI
      double Aqq = (fa != 0 ? fa*conform_qq[]: 1.); 
      Psiqq[] = log (Aqq); 
#endif

      /**
      The conformation tensor is diagonalized through the
      eigenvector tensor $\mathbf{R}$ and the eigenvalues diagonal
      tensor, $\Lambda$. */

      pseudo_v Lambda;
      pseudo_t R;
      diagonalization_2D (&Lambda, &R, &A);
      
      /**
      $\Psi = \log \mathbf{A}$ is easily obtained after diagonalization, 
      $\Psi = R \cdot \log(\Lambda) \cdot R^T$. */
      
      Psi.x.y[] = R.x.x*R.y.x*log(Lambda.x) + R.y.y*R.x.y*log(Lambda.y);
      foreach_dimension()
      	Psi.x.x[] = sq(R.x.x)*log(Lambda.x) + sq(R.x.y)*log(Lambda.y);

        /**
        We now compute the upper convective term $2 \mathbf{B} +
        (\Omega \cdot \Psi -\Psi \cdot \Omega)$.

        The diagonalization will be applied to the velocity gradient
        $(\nabla u)^T$ to obtain the antisymmetric tensor $\Omega$ and
        the traceless, symmetric tensor, $\mathbf{B}$. If the conformation
        tensor is $\mathbf{I}$, $\Omega = 0$ and $\mathbf{B}= \mathbf{D}$.  */

      pseudo_t B;
      double OM = 0.;
      if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
	B.x.y = (u.y[1,0] - u.y[-1,0] +
		 u.x[0,1] - u.x[0,-1])/(4.*Delta); 
	foreach_dimension() 
	  B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
      }
      else {
	pseudo_t M;
	foreach_dimension() {
	  M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
		   sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
		   R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + 
				u.y[1] - u.y[-1]))/(2.*Delta);
	  M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
		   R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
		   R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
		   R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
	}
	double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
	OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;
	
	B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
	foreach_dimension()
	  B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);	
      }

        /**
        We now advance $\Psi$ in time, adding the upper convective
        contribution. */

      double s = - Psi.x.y[];
      Psi.x.y[] += dt*(2.*B.x.y + OM*(Psi.y.y[] - Psi.x.x[]));
      foreach_dimension() {
	s *= -1;
	Psi.x.x[] += dt*2.*(B.x.x + s*OM);
      }

      /**
      In the axisymmetric case, the governing equation for $\Psi_{\theta
      \theta}$ only involves that component, 
      $$ 
      \Psi_{\theta \theta}|_t - 2 L_{\theta \theta} = 
      \frac{\mathbf{f}_r(e^{-\Psi_{\theta \theta}})}{\lambda} 
      $$
      with $L_{\theta \theta} = u_y/y$. Therefore step (a) for
      $\Psi_{\theta \theta}$ is */

#if AXI
      Psiqq[] += dt*2.*u.y[]/max(y, 1e-20);
#endif

}

  /**
  ### Advection of $\Psi$
  
  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

#if AXI
  advection ({Psi.x.x, Psi.x.y, Psi.y.y, Psiqq}, uf, dt);
#else
  advection ({Psi.x.x, Psi.x.y, Psi.y.y}, uf, dt);
#endif

    /**
    ### Convert back to \conform_p */

    foreach() {
      /**
      It is time to undo the log-conformation, again by
      diagonalization, to recover the conformation tensor $\mathbf{A}$
      and to perform step (c).*/

      pseudo_t A = {{Psi.x.x[], Psi.x.y[]}, {Psi.y.x[], Psi.y.y[]}}, R;
      pseudo_v Lambda;
      diagonalization_2D (&Lambda, &R, &A);
      Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);
      
      A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
      foreach_dimension()
        A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

      /**
      We perform now step (c) by integrating 
      $\mathbf{A}_t = -\mathbf{f}_r (\mathbf{A})/\lambda$ to obtain
      $\mathbf{A}^{n+1}$. This step is analytic,
      $$
      \int_{t^n}^{t^{n+1}}\frac{d \mathbf{A}}{\mathbf{I}- \mathbf{A}} = 
      \frac{\Delta t}{\lambda}
      $$
      */

     double intFactor = lambda[] != 0. ? exp(-dt/lambda[]): 0.;
     
#if AXI
      Aqq = (1. - intFactor) + intFactor*exp(Psiqq[]);
#endif

      A.x.y *= intFactor;
      foreach_dimension()
        A.x.x = (1. - intFactor) + A.x.x*intFactor;

      /**
        Then the Conformation tensor $\mathcal{A}_p^{n+1}$ is restored from
        $\mathbf{A}^{n+1}$.  */
      
      double fa = (f[] < 1e-6 ? 0.0: 1.0);
      
      conform_p.x.y[] = fa*A.x.y;
      tau_p.x.y[] = Gp[]*A.x.y;
#if AXI
      conform_qq[] = fa != 0.0 ? fa*(Aqq): 1.0;
      tau_qq[] = Gp[]*(Aqq - 1.);
#endif

      foreach_dimension(){
        conform_p.x.x[] = fa != 0.0 ? fa*(A.x.x): 1.0;
        tau_p.x.x[] = Gp[]*(A.x.x - 1.);
      }

  }
}

/**
### Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{\tau}_p$ is defined at cell centers
while the corresponding force (acceleration) will be defined at cell
faces. Two terms contribute to each component of the momentum
equation. For example the $x$-component in Cartesian coordinates has
the following terms: $\partial_x \mathbf{\tau}_{p_{xx}} + \partial_y
\mathbf{\tau}_{p_{xy}}$. The first term is easy to compute since it can be
calculated directly from center values of cells sharing the face. The
other one is harder. It will be computed from vertex values. The
vertex values are obtained by averaging centered values.  Note that as
a result of the vertex averaging cells `[]` and `[-1,0]` are not
involved in the computation of shear. */

event acceleration (i++)
{
  face vector av = a;
  foreach_face()
    if (fm.x[] > 1e-20) {
      double shear = (tau_p.x.y[0,1]*cm[0,1] + tau_p.x.y[-1,1]*cm[-1,1] -
		      tau_p.x.y[0,-1]*cm[0,-1] - tau_p.x.y[-1,-1]*cm[-1,-1])/4.;
      av.x[] += (shear + cm[]*tau_p.x.x[] - cm[-1]*tau_p.x.x[-1])*
	alpha.x[]/(sq(fm.x[])*Delta);
    }
#if AXI
  foreach_face(y)
    if (y > 0.)
      av.y[] -= (tau_qq[] + tau_qq[0,-1])*alpha.y[]/sq(y)/2.;
#endif
}

