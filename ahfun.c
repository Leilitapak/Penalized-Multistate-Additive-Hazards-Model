#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP get_coef(SEXP cft, SEXP del, SEXP z, SEXP nsub) {
	R_len_t	n = INTEGER(nsub)[0], n0 = nrows(z), p = ncols(z), i, i0, j;
	int		*r_del = INTEGER(del), *ind;
	double	*r_cft = REAL(cft), *r_z = REAL(z), *r_b, *r_x, *r_v, *work, s, t;
	SEXP	ans;
	
	ind = (int *) R_alloc(n, sizeof(int));
	work = (double *) R_alloc(n, sizeof(double));
	for (i = 0; i < n; i++) {
		ind[i] = i; work[i] = r_cft[i];
	}
	R_qsort_I(work, ind, 1, n);
	
	PROTECT(ans = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, p)); r_b = REAL(VECTOR_ELT(ans, 0));
	SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP, n, p)); r_x = REAL(VECTOR_ELT(ans, 1));
	SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, p)); r_v = REAL(VECTOR_ELT(ans, 2));
	for (j = 0; j < p; j++) {
		s = 0.0;
		for (i = n - 1; i >= 0; i--) {
			s += r_z[ind[i] + j*n0]; work[i] = s/(n - i);
		}
		s = 0.0; t = 0.0; r_v[j] = 0.0;
		for (i = 0; i < n; i++) {
			i0 = ind[i];
			if (r_del[i0]) s += r_z[i0 + j*n0] - work[i];
			t += work[i]*(i > 0 ? r_cft[i0] - r_cft[ind[i - 1]] : r_cft[i0]);
			r_x[i0 + j*n] = (r_z[i0 + j*n0]*r_cft[i0] - t)/n;
			r_v[j] += r_x[i0 + j*n]*r_z[i0 + j*n0];
		}
		r_b[j] = s/n;
	}
	UNPROTECT(1);
	return ans;
}

SEXP get_pred_err(SEXP b, SEXP x, SEXP z, SEXP sol) {
	R_len_t	n = nrows(x), n0 = nrows(z), p = ncols(z),
			n_lam = INTEGER(getAttrib(sol, R_DimSymbol))[1],
			n_a = INTEGER(getAttrib(sol, R_DimSymbol))[2], i, j, k, m;
	double	*r_b = REAL(b), *r_x = REAL(x), *r_z = REAL(z), *r_sol = REAL(sol), *r_ans,
			s, s1, s2, t;
	SEXP	ans;
	
	PROTECT(ans = allocMatrix(REALSXP, n_lam, n_a)); r_ans = REAL(ans);
	for (j = 0; j < n_lam; j++) for (k = 0; k < n_a; k++) {
		s = 0.0;
		for (m = 0; m < n; m++) {
			s1 = 0.0; s2 = 0.0;
			for (i = 0; i < p; i++) {
				s1 += r_x[m + i*n]*r_sol[i + j*p + k*p*n_lam];
				s2 += r_z[m + i*n0]*r_sol[i + j*p + k*p*n_lam];
			}
			s += s1*s2;
		}
		t = 0.0;
		for (i = 0; i < p; i++) t += r_b[i]*r_sol[i + j*p + k*p*n_lam];
		r_ans[j + k*n_lam] = 0.5*s - t;
	}
	UNPROTECT(1);
	return ans;
}

SEXP coord_descent(SEXP b, SEXP x, SEXP z, SEXP v, SEXP method, SEXP lam_seq, SEXP a_seq,
				   SEXP upto) {
	R_len_t n = nrows(x), n0 = nrows(z), p = ncols(z), n_lam = length(lam_seq),
			n_a = length(a_seq), i, j, k, m;
	int		pen, iter, max_iter = 50, s, brk;
	double	*r_b = REAL(b), *r_x = REAL(x), *r_z = REAL(z), *r_v = REAL(v),
			*r_lam_seq = REAL(lam_seq), *r_a_seq = REAL(a_seq), max_s = REAL(upto)[0], lam, a,
			a1=NA_REAL, a2=NA_REAL, b_hat, t0, a_t0, t1, t2, c0, c1, c2, q, r, rt_q, alp,
 			old, dif, err, tol = 1e-6, *r_sol, *bet, *z_bet, *z_bet0;
	SEXP	ans;
	
	PROTECT(ans = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ans, 0, alloc3DArray(REALSXP, p, n_lam, n_a));
	r_sol = REAL(VECTOR_ELT(ans, 0));
	SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, 1));
	
	bet = (double *) R_alloc(p, sizeof(double));
	z_bet = (double *) R_alloc(n, sizeof(double));
	z_bet0 = (double *) R_alloc(n, sizeof(double));
	
	for (i = 0; i < p; i++) bet[i] = 0.0;
 	for (i = 0; i < n; i++) z_bet[i] = 0.0;
	brk = 0;
	for (j = 0; j < n_lam && !brk; j++) {
		lam = r_lam_seq[j];
		if (j > 0) {
			for (i = 0; i < p; i++) bet[i] = r_sol[i + (j - 1)*p];
 			for (i = 0; i < n; i++) z_bet[i] = z_bet0[i];
		}
		for (k = 0; k < n_a; k++) {
			a = r_a_seq[k];
			pen = INTEGER(method)[0];
			switch(pen) {
				case 1:		/* SCAD */
					if (R_FINITE(a)) {
						a1 = lam*a/(a - 1.0); a2 = (a - 1.0)/(a - 2.0);
					} else
						pen = 0;
					break;
				case 2:		/* MCP */
					if (R_FINITE(a)) a1 = a/(a - 1.0); else pen = 0;
					break;
				case 3:		/* SICA */
					if (R_FINITE(a)) {
						a2 = a*a; a1 = lam*(a2 + a);
					} else
						pen = 0;
					break;
				case 4:		/* Enet */
					a1 = lam*a; a2 = 1.0/(1.0 + lam*(1.0 - a)); break;
			}
			for (iter = 0; iter < max_iter; iter++) {
				err = 0.0;
				for (i = 0; i < p; i++) {
					old = bet[i]; b_hat = 0.0;
					for (m = 0; m < n; m++) b_hat += r_x[m + i*n]*z_bet[m];
					t0 = (r_b[i] - b_hat)/r_v[i] + bet[i]; a_t0 = fabs(t0);
					switch(pen) {
						case 0:		/* Lasso */
							bet[i] = copysign(fdim(a_t0, lam), t0); break;
						case 1:		/* SCAD */
							if (a_t0 > a*lam)
								bet[i] = t0;
							else if (a_t0 > 2.0*lam)
								bet[i] = copysign((a_t0 - a1)*a2, t0);
							else
								bet[i] = copysign(fdim(a_t0, lam), t0);
							break;
						case 2:		/* MCP */
							if (a_t0 > a*lam)
								bet[i] = t0;
							else
								bet[i] = copysign(fdim(a_t0, lam)*a1, t0);
							break;
						case 3:		/* SICA */
							c2 = 2.0*a - a_t0; c1 = a2 - 2.0*a*a_t0; c0 = a1 - a2*a_t0;
							q = 1.0/9*c2*c2 - 1.0/3*c1;
							r = 1.0/27*c2*c2*c2 - 1.0/6*c1*c2 + 0.5*c0;
							if (q*q*q <= r*r)
								bet[i] = 0.0;
							else {
								rt_q = sqrt(q);
								alp = 1.0/3*acos(r/(rt_q*rt_q*rt_q));
								t1 = -2.0*rt_q*cos(alp - 2.0/3*M_PI) - 1.0/3*c2;
								t2 = -2.0*rt_q*cos(alp + 2.0/3*M_PI) - 1.0/3*c2;
								if (t1 > 0.0)
									if (0.5*t2 + lam*(a + 1.0)/(a + t2) < a_t0)
										bet[i] = copysign(t2, t0);
									else
										bet[i] = 0.0;
								else
									bet[i] = copysign(fmax(t2, 0.0), t0);
							}
							break;
						case 4:		/* Enet */
							bet[i] = copysign(fdim(a_t0, a1), t0)*a2; break;
					}
					if (bet[i] != old) {
						dif = bet[i] - old;
						for (m = 0; m < n; m++) z_bet[m] += r_z[m + i*n0]*dif;
						err = fmax(err, fabs(dif));
					}
				}
				if (err < tol) break;
			}
			for (i = 0; i < p; i++)	r_sol[i + j*p + k*p*n_lam] = bet[i];
			if (k == 0) {
				s = 0;
				for (i = 0; i < p; i++) if (bet[i] != 0.0) s++;
				if (R_FINITE(max_s) && s > max_s)
					brk = 1;
				else
					for (i = 0; i < n; i++) z_bet0[i] = z_bet[i];
			}
		}
	}
	for (i = 0; i < p; i++) for (m = j; m < n_lam; m++) for (k = 0; k < n_a; k++)
		r_sol[i + m*p + k*p*n_lam] = NA_REAL;
	
	INTEGER(VECTOR_ELT(ans, 1))[0] = j;
	UNPROTECT(1);
	return ans;
}

SEXP cross_valid(SEXP cft, SEXP del, SEXP z, SEXP method, SEXP lam_seq, SEXP a_seq, SEXP ind) {
	R_len_t	n = nrows(z), p = ncols(z), n_lam = length(lam_seq), n_a = length(a_seq),
			n_fold = length(ind), max_n1, max_n2, i, j, j0, j1, j2, k, m;
	int		*r_del = INTEGER(del), *r_del1, *r_del2, *r_i, *r_j;
	double	*r_cft = REAL(cft), *r_z = REAL(z), *r_cv, *r_se, *r_cft1, *r_cft2, *r_z1, *r_z2, cv0;
	SEXP	ans, upto, cft1, cft2, del1, del2, z1, z2, n1, n2, c1, c2, sol, pe;
	
	PROTECT(ans = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP, n_lam, n_a));
	r_cv = REAL(VECTOR_ELT(ans, 0));
	SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP, n_lam, n_a));
	r_se = REAL(VECTOR_ELT(ans, 1));
	SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, 1));
	r_i = INTEGER(VECTOR_ELT(ans, 2));
	SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
	r_j = INTEGER(VECTOR_ELT(ans, 3));
	PROTECT(upto = allocVector(REALSXP, 1)); REAL(upto)[0] = R_PosInf;
	PROTECT(n1 = allocVector(INTSXP, 1)); PROTECT(n2 = allocVector(INTSXP, 1));
	
	max_n1 = 0; max_n2 = 0;
	for (i = 0; i < n_fold; i++) {
		max_n1 = imax2(max_n1, n - length(VECTOR_ELT(ind, i)));
		max_n2 = imax2(max_n2, length(VECTOR_ELT(ind, i)));
	}
	PROTECT(cft1 = allocVector(REALSXP, max_n1)); r_cft1 = REAL(cft1);
	PROTECT(cft2 = allocVector(REALSXP, max_n2)); r_cft2 = REAL(cft2);
	PROTECT(del1 = allocVector(INTSXP, max_n1)); r_del1 = INTEGER(del1);
	PROTECT(del2 = allocVector(INTSXP, max_n2)); r_del2 = INTEGER(del2);
	PROTECT(z1 = allocMatrix(REALSXP, max_n1, p)); r_z1 = REAL(z1);
	PROTECT(z2 = allocMatrix(REALSXP, max_n2, p)); r_z2 = REAL(z2);
	
	for (i = 0; i < n_lam; i++) for (j = 0; j < n_a; j++) {
		r_cv[i + j*n_lam] = 0.0; r_se[i + j*n_lam] = 0.0;
	}
	for (k = 0; k < n_fold; k++) {
		j1 = 0; j2 = 0;
		for (i = 0; i < n_fold; i++)
			if (i != k)
				for (j = 0; j < length(VECTOR_ELT(ind, i)); j++) {
					j0 = INTEGER(VECTOR_ELT(ind, i))[j] - 1;
					r_cft1[j1] = r_cft[j0]; r_del1[j1] = r_del[j0];
					for (m = 0; m < p; m++) r_z1[j1 + m*max_n1] = r_z[j0 + m*n];
					j1++;
				}
			else
				for (j = 0; j < length(VECTOR_ELT(ind, i)); j++) {
					j0 = INTEGER(VECTOR_ELT(ind, i))[j] - 1;
					r_cft2[j2] = r_cft[j0]; r_del2[j2] = r_del[j0];
					for (m = 0; m < p; m++) r_z2[j2 + m*max_n2] = r_z[j0 + m*n];
					j2++;
				}
		INTEGER(n1)[0] = n - length(VECTOR_ELT(ind, k));
		INTEGER(n2)[0] = length(VECTOR_ELT(ind, k));
		PROTECT(c1 = get_coef(cft1, del1, z1, n1));
		PROTECT(sol = VECTOR_ELT(coord_descent(VECTOR_ELT(c1, 0), VECTOR_ELT(c1, 1), z1,
								 VECTOR_ELT(c1, 2), method, lam_seq, a_seq, upto), 0));
		PROTECT(c2 = get_coef(cft2, del2, z2, n2));
		PROTECT(pe = get_pred_err(VECTOR_ELT(c2, 0), VECTOR_ELT(c2, 1), z2, sol));
		for (i = 0; i < n_lam; i++) for (j = 0; j < n_a; j++) {
			r_cv[i + j*n_lam] += REAL(pe)[i + j*n_lam];
			r_se[i + j*n_lam] += REAL(pe)[i + j*n_lam]*REAL(pe)[i + j*n_lam];
		}
		UNPROTECT(4);
	}
	
	cv0 = R_PosInf;
	for (i = 0; i < n_lam; i++) for (j = 0; j < n_a; j++) {
		r_cv[i + j*n_lam] /= n_fold;
		r_se[i + j*n_lam] =
			sqrt((r_se[i + j*n_lam] - n_fold*r_cv[i + j*n_lam]*r_cv[i + j*n_lam])/(n - 1));
		if (r_cv[i + j*n_lam] < cv0){
			cv0 = r_cv[i + j*n_lam]; r_i[0] = i; r_j[0] = j;
		}
	}
	r_i[0]++; r_j[0]++;
	UNPROTECT(10);
	return ans;
}
