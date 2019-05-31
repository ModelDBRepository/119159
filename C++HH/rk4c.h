void     
runge_kutta (double *x, double t, double tau,  double *xout){
double F1[NDIM], F2[NDIM], F3[NDIM], F4[NDIM], xtemp[NDIM];
double t_half, half_tau, t_full;
int  i;

  half_tau = 0.5 * tau;
  model (x, t, F1);
  t_half = t + half_tau;
  for (i = 0;i < NDIM; i++)
    xtemp[i] = x[i] + half_tau*F1[i];
  model (xtemp, t_half, F2);
  for (i = 0;i < NDIM; i++)
    xtemp[i] = x[i] + half_tau*F2[i];
  model (xtemp, t_half, F3);
  t_full = t + tau;
  for (i = 0;i < NDIM; i++)
    xtemp[i] = x[i] + tau * F3[i];
  model (xtemp, t_full, F4);
  for (i = 0;i < NDIM; i++)
    xout[i] = x[i] + tau/6.0 *(F1[i] + F4[i] + 2.0*(F2[i] + F3[i]));
}
