/****************************************************************

file elkies_allg.c

Elkies package for the three cubes problem,
Version 1.0

Copyright (C) 2007 A.-S. Elsenhans and J. Jahnel

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Address of the authors

	A.-S. Elsenhans/ J. Jahnel
	Math. Institut der Universitaet
	Bunsenstrasze 3--5
	D-37073 Goettingen
        Germany

WWW     http://www.uni-math.gwdg.de/jahnel

*****************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<gmp.h>
#include<math.h>
#include"festkomma.h"

double        seeking_wide  =     1.0E14;
#define        UNDERCARRIAGE     1.0E11

/* Tripel mit |x**3 + y**3 - z**3| < EDITION_BARRIER werden ausgegeben. */
#define        EDITION_BARRIER  1000


/* The global variables half step width, half tile width and
   Tile offsets are calculated in calculated tile format.
   This function will execute all NEW_TILE runs of the outer loop
   called.
   The FACTOR controls the size of the tile.

   Note, hereby:
   Half step size is half the extent of a tile in the x direction.
   Half tile width is half the extent of a tile in the y direction.
   It must apply:
      half tile width = 1.001 * y '' (x) * half step size ** 2/4
   With
      y (x): = (1 - x ** 3) ** (1/3).
   The content of the tile is thus 1.001 * y '' (x) * half step size ** 3.

   This value is set to 1,001 * (FACTOR / SEARCH) ** 3.
   From FACTOR, SEARCH and the current second derivative are thus
   half step width and half tile width calculable. */
#define        NEW_TILE        1000000
#define        FAKTOR            5.5
double  half_step, half_tile_width;
mpx_t                   x_0, tiles_offset, Ax;
mpx_t                 y_diff, y_inv, y_inv_diff;
char                      ausg[1000], file[100];
long                                      row;


/* 
Count for the total number of tiles processed. */
long                                    tile;

/* These mpf variables should themselves be local variables, mainly in the
   Functions calculate_y_value and calculate_three_linearf respectively
   calculate_three_linearf_mpf.
   But with global variables, the code runs faster.
 */
mpf_t                                tmp1, tmp2;


/* Only a prototype. */
void init (long v[3][3], mpx_t y_0_type, mpx_t x_0_type,
  double half_step, mpx_t tiles_offset,
  double half_tile_width, double seeking_wide);


void out () {
 FILE *fp;

 fp = fopen (file, "a");
 fprintf (fp, "%6ld \t ", row);
 row++;
 fprintf (fp, "%s", ausg);
 fclose (fp);
}


/* pr: = m1 * m2. Product of two 3x3 matrices with integer entries.
   We calculate in long and pay attention to overflow over \ pm2 ** 63.
   pr = m1 and pr = m2 are allowed. */
inline long matrix_prod (long pr[3][3],
                         long m1[3][3], long m2[3][3]) {
 long               i, j;
 long         prod[3][3];
 long  prodh, argh, argl;
 long                err;

 err = 0;
 for (i = 0; i < 3; i++) {
  for (j = 0; j < 3; j++) {
   smul_ppmm (prodh, prod[i][j], m1[i][0], m2[0][j]);
   smul_ppmm (argh, argl, m1[i][1], m2[1][j]);
   addf_ssaaaa (prodh, prod[i][j], prodh, prod[i][j], argh, argl);
   smul_ppmm (argh, argl, m1[i][2], m2[2][j]);
   addf_ssaaaa (prodh, prod[i][j], prodh, prod[i][j], argh, argl);
   if (((prod[i][j] >= 0) && (prodh != 0))
     || ((prod[i][j] < 0) && (prodh != -1))) {
    err = 1; /* Overflow! */
   }
  }
 }
 /* Copying. Needed in the event that pr to the address of m1 or m2
    should be written. */
 for (i = 0; i < 3; i++) {
  for (j = 0; j < 3; j++) {
   pr[i][j] = prod[i][j];
  }
 }
 return (err);
}


/* Copying. Needed in the event that pr to the address of m1 or m2
    should be written.*/
void calculate_tile (mpx_t increment) {
 double  x, zweite_abl, additional_summand;
 long                              fl;

 /* Output only every 100 * Tile_NEW tile. */
 fl = tile % (100 * NEW_TILE);

 x = mpx_get_d (x_0);


 zweite_abl = 2*x / pow (1 - x*x*x, 5.0/3.0);
 half_step = pow (zweite_abl, -1.0/3.0) * FAKTOR / seeking_wide;

 additional_summand = EDITION_BARRIER
                              / (UNDERCARRIAGE * UNDERCARRIAGE * UNDERCARRIAGE);
 half_tile_width =  1.001 * zweite_abl *
                                half_step*half_step / 4
                      + additional_summand;
                     /* additional_summand is designed to ensure that solutions work with you
                        Height> UNDERCARRIAGE guaranteed to be found. */
 if (fl == 0) {
  sprintf (ausg, "half_step = %.9e, half_tile_width = %.9e.\n",
                               half_step, half_tile_width);
  out ();
 }

 mpx_set_d (increment, 2 * half_step);
 mpx_set_d (tiles_offset,
                  zweite_abl * half_step*half_step / 4);
}


/***************************************************************************
 *
 * Eigentlicher Code
 *
 **************************************************************************/


/* Called once per NEW_TILE tile.
    Computes 1 / y_0 (- 1) consuming and exact. */
void y_inv_init (mpx_t y_0) {
 /* Initialization of y_inv. */
 mpf_set_ui (tmp1, 1);
 mpf_set_mpx (tmp2, y_0);
 mpf_div (tmp1, tmp1, tmp2);
 mpf_sub_ui (tmp1, tmp1, 1);
 mpx_set_mpf (y_inv, tmp1);
 /* gmp_sprintf (ausg, "y_inv = %.*Ff.\n", 40, tmp1); */
}


/* Called once per NEW_TILE tile.
    Initializes y_diff as y '(x_0) * increment.
    Computes y_inv_diff as y_diff / y_0 ** 2. */
void y_diff_init (mpx_t y_0, mpx_t increment) {
 double  diff, y;

/* Initialization of y_diff. * /
 mpx_mul (y_diff, Ax, increment);
/* In case y_0 <x_0 is in fact Ax> 1. Correct here after. * /
 if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
  mpx_add (y_diff, y_diff, increment);
 /* mpf_set_mpx (tmp1, y_diff);
 gmp_sprintf (ausg, "y_diff = %.*Ff.\n\n", 40, tmp1); */

/* Initialization of y_inv_diff. * /
 y = mpx_get_d (y_0);
 diff = mpx_get_d (y_diff) / (y*y);
 mpx_set_d (y_inv_diff, diff);
 /* mpf_set_mpx (tmp1, y_inv_diff);
 gmp_sprintf (ausg, "y_inv_diff = %.*Ff.\n", 40, tmp1); */
}


/* y_0 := (1 - x_0**3)**(1/3). */
inline void berechne_y_wert (mpx_t y_0) {
 double        nenner, diff;
 mpx_t  tmpx1, tmpx2, tmpx3;

 mpx_add (y_0, y_0, y_diff);
 /* y_0 - = increment * A. Linear improvement of the starting value.
     Minus, because we walk backwards through the interval.
     y_diff is calculated once per 1000000 tile.
     A is always negative, ie (-increment) * A> 0.
     Note, y_diff> 0 after initialization.. */

 mpx_mul (tmpx1, x_0, x_0);
 mpx_mul (tmpx2, tmpx1, x_0);
 tmpx3[0] = ~tmpx2[0]; tmpx3[1] = ~tmpx2[1];
 /* y3 = 1 - x_0*x_0*x_0; */
 /* We make a mistake of exactly 2 ** (- 128). */

 mpx_mul (tmpx1, y_0, y_0);
 nenner = 3.0 * mpx_get_d (tmpx1);

 mpx_mul (tmpx2, y_0, tmpx1);
 mpx_sub (tmpx1, tmpx2, tmpx3);
 diff = mpx_get_d (tmpx1);
 if (diff > 0.5) {
  sprintf (ausg, "Underflow.\n"); out ();
  exit (0);
 }
 diff /= nenner;

 mpx_set_d (tmpx1, diff);
 mpx_sub (y_0, y_0, tmpx1);
 /* A Newton-Iteration. */
 /* y_0 = (2*y_0 + y3/(y_0*y_0)) / 3
        = (2*y_0*y_0*y_0 + y3) / (3*y_0*y_0)
        = y_0 + (y3 - y_0*y_0*y_0) / (3*y_0*y_0);
   Difference can be calculated in double. * /
  /* A Newton iteration. * /

 /* mpf_set_mpx (tmp1, y_0);
 gmp_sprintf (ausg, "y_0 = %.*Ff.\n", 40, tmp1); out (); */
}


/* Write y '(x_0) into the global variable Ax. * /
inline void berechne_y_strich (mpx_t y_0) {
 mpx_t  tmpx1, tmpx2;

/* First calculate 1 / y_0. * /

  /* Linear improvement of the starting value. * /
 mpx_sub (y_inv, y_inv, y_inv_diff);

/* Now Newton iteration inv_new = inv + (1 - y_0 * inv) * inv.
  Since y_inv approaches the value 1 / y - 1, we have to
       y_inv = y_inv + (1 - y_0 * (1 + y_inv)) * (1 + y_inv)
             = y_inv + (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) * y_inv
  expected. * /
 mpx_mul (tmpx1, y_0, y_inv);
 mpx_add (tmpx2, tmpx1, y_0); /* tmpx2 ist y_0 + y_0 * y_inv. */
 tmpx2[0] = ~tmpx2[0]; tmpx2[1] = ~tmpx2[1];
 /* tmpx2 ist jetzt 1 - y_0 - y_0 * y_inv. Fehler von 2**(-128) ignoriert. */
 mpx_mul (tmpx1, tmpx2, y_inv);
 mpx_add (tmpx1, tmpx1, tmpx2);
 /* tmpx1 ist (1 - y_0 - y_0 * y_inv) + (1 - y_0 - y_0 * y_inv) * y_inv. */
 mpx_add (y_inv, y_inv, tmpx1);
 /* mpf_set_mpx (tmp1, y_inv);
 gmp_sprintf (ausg, "y_inv = %.*Ff.\n", 40, tmp1); */

/* Calculate now
         y'(x_0) = Ax := x_0**2 / y_0**2
                 = x_0**2 * (1 + y_inv)**2
                 = x_0**2 + x_0**2 * (y_inv**2 + 2*y_inv). */
 mpx_mul (tmpx1, y_inv, y_inv);
 mpx_add (tmpx1, tmpx1, y_inv);
 mpx_add (tmpx1, tmpx1, y_inv); /* tmpx1 ist jetzt y_inv**2 + 2*y_inv. */
 mpx_mul (tmpx2, x_0, x_0); /* tmpx2 ist x_0**2. */
 mpx_mul (Ax, tmpx2, tmpx1); /* Ax ist jetzt x_0**2 * (y_inv**2 + 2*y_inv). */
 mpx_add (Ax, Ax, tmpx2);
 /* mpf_set_mpx (tmp1, Ax);
 gmp_sprintf (ausg, "A = %.*Ff.\n", 40, tmp1); */

 /* Ax = - x_0*x_0 / (y_0*y_0); */
/* Derivation of
        y(x) := (1 - x**3)**(1/3)
after the sentence about implicit functions. * /
 /* A = (y_1 - y_0) / (x_1 - x_0); */
 /* A = - x_0*x_0 / pow (1 - x_0*x_0*x_0, 2.0/3.0); */
}


inline void calculate_three_linearf
  (double l[3][3], long v[3][3], mpx_t y_0, double d, double N) {
 long                         i;
 mpx_t  tmpx1, tmpx2, tmpx3, Bx;


 /* mpf_set_mpx (tmp1, x_0);
 gmp_sprintf (ausg, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */
 mpx_mul (tmpx1, Ax, x_0);
 /* Im Falle y_0 < x_0 ist in Wirklichkeit Ax > 1. Korrigieren hier nach. */
 if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
  mpx_add (tmpx1, tmpx1, x_0);

 /* mpf_set_mpx (tmp1, tiles_offset);
 gmp_sprintf (ausg, "tiles_offset = %.*Ff.\n", 40, tmp1); out (); */
 mpx_sub (tmpx2, y_0, tiles_offset);
 mpx_add (Bx, tmpx1, tmpx2); /* Ax ist in Wirklichkeit negativ. */
 /* Bx = y_0 - tiles_offset - Ax * x_0; */
 /* Es gilt 1 <= Bx <= 1.6. Speichern in Wirklichkeit Bx - 1. */
 /* mpf_set_mpx (tmp1, Ax); mpf_set_mpx (tmp2, Bx);
 gmp_sprintf (ausg, "A = %.*Ff.\nB = %.*Ff.\n", 40, tmp1, 40, tmp2); out (); */

 for (i = 0; i < 3; i++) {
  l[0][i] = v[i][2] / N;
  /* sprintf (ausg, "l[0][%ld] = %f.\n", i, l[0][i]); out (); */
  /* No problems with accuracy here. * /

  mpx_mul_si (tmpx1, x_0, -v[i][2]);
  mpxb_add_si (tmpx2, tmpx1, v[i][0]);
  /* tmpx = v[i][0] - v[i][2] * x_0; */
  l[1][i] = mpxg_get_d_simple (tmpx2);
  l[1][i] /= (half_step * N);
  /* sprintf (ausg, "l[1][%ld] = %f.\n", i, l[1][i]); out (); */
/* It is increment * N \ approx 100.
      So we need tmp2 at> 2 (7?) Decimal places.
      For this double should be enough.
      tmp2 is on 128 bits, so> = 23 decimal places exactly. * /

  mpx_mul_si (tmpx1, Ax, v[i][0]);
/* In case y_0 <x_0 is in fact Ax> 1. Correct here after. * /
  if ((y_0[1] < x_0[1]) || ((y_0[1] == x_0[1]) && (y_0[0] <= x_0[0])))
   mpxb_add_si (tmpx1, tmpx1, v[i][0]);

  mpx_mul_si (tmpx2, Bx, -v[i][2]);
  mpx_add (tmpx3, tmpx1, tmpx2); /* Ax ist in Wirklichkeit negativ. */
  mpxb_add_si (tmpx1, tmpx3, v[i][1] - v[i][2]); /* Bx \in [1..2]. */
  /* tmpx = v[i][1] - Ax * v[i][0] - Bx * v[i][2]; */
  l[2][i] = mpxg_get_d (tmpx1);
  l[2][i] /= (d*N);
  /* sprintf (ausg, "l[2][%ld] = %f.\n", i, l[2][i]); out (); */
  /* l[2][i] most sensitive to too little precision.
      It's d * N \ approx 1.0e-15.
      So we need tmp2 at> 15 (20?) Decimal places behind the comma.
      So Ax and Bx to 35 decimal places, corresponding to 128 bits. * /
 }
}


/* Makros fuer lll. */
#define scal_prod(prod, l, vec1, vec2)                            \
do {                                                              \
 long          t1, t2;                                            \
 double  l1[3], l2[3];                                            \
                                                                  \
 for (t1 = 0; t1 < 3; t1++)                                       \
  l1[t1] = l2[t1] = 0;                                            \
 for (t1 = 0; t1 < 3; t1++) {                                     \
  for (t2 = 0; t2 < 3; t2++) {                                    \
   l1[t1] += l[t1][t2] * vec1[t2];                                \
   l2[t1] += l[t1][t2] * vec2[t2];                                \
  }                                                               \
 }                                                                \
 prod = 0;                                                        \
 for (t1 = 0; t1 < 3; t1++)                                       \
  prod += l1[t1]*l2[t1];                                          \
} while (0);



#define gram()                                                    \
 do {                                                             \
  long  t1;                                                       \
                                                                  \
  /* sprintf (ausg, "Gram!\n"); out (); */                        \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[0][t1] = vec[0][t1];      /* double = long */         \
  scal_prod (B[0], l, vec_gram[0], vec_gram[0]);                  \
  scal_prod (mu[1][0], l, vec[1], vec_gram[0]); mu[1][0] /= B[0]; \
                                                                  \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[1][t1] = vec[1][t1] - mu[1][0] * vec_gram[0][t1];     \
  scal_prod (B[1], l, vec_gram[1], vec_gram[1]);                  \
  scal_prod (mu[2][0], l, vec[2], vec_gram[0]); mu[2][0] /= B[0]; \
  scal_prod (mu[2][1], l, vec[2], vec_gram[1]); mu[2][1] /= B[1]; \
                                                                  \
  for (t1 = 0; t1 < 3; t1++)                                      \
   vec_gram[2][t1] = vec[2][t1] - mu[2][0] * vec_gram[0][t1]      \
                                - mu[2][1] * vec_gram[1][t1];     \
  scal_prod (B[2], l, vec_gram[2], vec_gram[2]);                  \
} while (0)


#define red_k_k1()                                                \
 do {                                                             \
  /* sprintf (ausg, "Mache RED (%ld, %ld).\n", k, k-1); out (); */\
  tm = floor (mu[k][k-1] + 0.5);                                  \
  /* if (fabs (tm) > (1LU << 63) - 1) {                           \
   sprintf (ausg, "q is too long.\n"); out ();                    \
   exit (0);                                                      \
  } */                                                            \
  q = lround (tm);                                                \
  if (q != 0) {                                                   \
   for (i = 0; i < 3; i++)                                        \
    vec[k][i] -= q * vec[k-1][i];                                 \
   mu[k][k-1] -= tm;                                              \
   if (k == 2)                                                    \
    mu[2][0] -= tm * mu[1][0];                                    \
  }                                                               \
} while (0)


#define red_2_0()                                                 \
 do {                                                             \
  /* sprintf (ausg, "Mache RED (2, 0).\n"); out (); */            \
  tm = floor (mu[2][0] + 0.5);                                    \
  q = lround (tm);                                                \
  /* if (fabs (tm) > (1LU << 63) - 1) {                           \
   sprintf (ausg, "q is too long.\n"); out ();                    \
   exit (0);                                                      \
  } */                                                            \
  if (q != 0)                                                     \
   for (i = 0; i < 3; i++)                                        \
    vec[2][i] -= q * vec[0][i];                                   \
} while (0)


#define lll_swap()                                                \
do {                                                              \
 /* sprintf (ausg, "Mache SWAP (%ld).\n", k); out (); */          \
 for (i = 0; i < 3; i++) {                                        \
  tmp = vec[k][i];                                                \
  vec[k][i] = vec[k-1][i];                                        \
  vec[k-1][i] = tmp;                                              \
 }                                                                \
                                                                  \
 muc = mu[k][k-1];                                                \
 Bc = B[k] + muc*muc * B[k-1];                                    \
 mu[k][k-1] = muc * B[k-1] / Bc;                                  \
 B[k] = B[k-1]*B[k] / Bc;                                         \
 B[k-1] = Bc;                                                     \
                                                                  \
 if (k == 1) {                                                    \
  t = mu[2][1];                                                   \
  mu[2][1] = mu[2][0] - muc*t;                                    \
  mu[2][0] = t + mu[1][0]*mu[2][1];                               \
 }                                                                \
 else { /* k == 2*/                                               \
  tm = mu[2][0];                                                  \
  mu[2][0] = mu[1][0];                                            \
  mu[1][0] = tm;                                                  \
 }                                                                \
 k = 1;                                                           \
} while (0)


/* LLL -- double-Version.
The function expects the three linear forms as given in l.
    It returns in v a set of vectors short for l.

    The process of the algorithm is (except for an error in gram) by H. Cohen
    accepted.

    We apply it in such a way that it already applies to a particular one
    Base is calculated. The return matrix vec is then one
    Transformation matrix to a new base of short vectors.

    ATTENTION: This function only works if l is not too big
    consists. We walk along the curve and hope that the vectors, the
    were short at last tile, now are not very long. * /
inline void lll (double l[3][3], long vec[3][3]) {
 long                    i, k;
 double  B[3], vec_gram[3][3];
 double              mu[3][3];
 long                  q, tmp;
 double        muc, Bc, t, tm;

/* Standard basis. The first candidate for the transformation matrix. */
 vec[0][0] = 1;  vec[0][1] = 0;  vec[0][2] = 0;
 vec[1][0] = 0;  vec[1][1] = 1;  vec[1][2] = 0;
 vec[2][0] = 0;  vec[2][1] = 0;  vec[2][2] = 1;

 /* LLL. */
 k = 1;
 gram ();

 while (k < 3) {
  /* sprintf (ausg, "%ld\n", k); out (); */
  red_k_k1 ();

  /* Test LLL condition */
  if (B[k] < ((0.99 - mu[k][k-1]*mu[k][k-1]) * B[k-1])) {
  /* if (B[k] < 0.5 * B[k-1]) { */
   /* sprintf (ausg, "Swap mit k = %ld \n", k); out (); */
   lll_swap ();
   /* out (); */
  }
  else {
   if (k == 2)
    red_2_0 ();
   k++;
  }
 }
}


#define lf_neu(lf, l, e)                                          \
do {                                                              \
 /* long  t1, t2;                                                 \
 for (t1 = 0; t1 < 3; t1++)                                       \
  for (t2 = 0; t2 < 3; t2++)                                      \
   lf[t1][t2]                                                     \
   = l[t1][0]*e[t2][0] + l[t1][1]*e[t2][1] + l[t1][2]*e[t2][2]; */\
                                                                  \
 lf[0][0] = l[0][0]*e[0][0] + l[0][1]*e[0][1] + l[0][2]*e[0][2];  \
 lf[0][1] = l[0][0]*e[1][0] + l[0][1]*e[1][1] + l[0][2]*e[1][2];  \
 lf[0][2] = l[0][0]*e[2][0] + l[0][1]*e[2][1] + l[0][2]*e[2][2];  \
                                                                  \
 lf[1][0] = l[1][0]*e[0][0] + l[1][1]*e[0][1] + l[1][2]*e[0][2];  \
 lf[1][1] = l[1][0]*e[1][0] + l[1][1]*e[1][1] + l[1][2]*e[1][2];  \
 lf[1][2] = l[1][0]*e[2][0] + l[1][1]*e[2][1] + l[1][2]*e[2][2];  \
                                                                  \
 lf[2][0] = l[2][0]*e[0][0] + l[2][1]*e[0][1] + l[2][2]*e[0][2];  \
 lf[2][1] = l[2][0]*e[1][0] + l[2][1]*e[1][1] + l[2][2]*e[1][2];  \
 lf[2][2] = l[2][0]*e[2][0] + l[2][1]*e[2][1] + l[2][2]*e[2][2];  \
} while (0);


/******************************************************************************
 *
 * The code for finding grid points in the pyramid.
 *
 *****************************************************************************/

/* res = m * v, 3x3 matrix times column vector. */
#define MMUL(res, m, v)                                           \
do {                                                              \
 /* long  k;                                                      \
                                                                  \
 for (k = 0; k < 3; k++)                                          \
  res[k] = m[k][0]*v[0] + m[k][1]*v[1] + m[k][2]*v[2]; */         \
                                                                  \
 res[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];             \
 res[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];             \
 res[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];             \
} while (0)


/* Berechne die vier Ecken der Pyramide. */
inline void pyr_ecken
 (double pt1[3], double pt2[3], double pt3[3], double pt4[3], double lf[3][3]) {
 long                             i, j;
 double  det_inv, adj[3][3], inv[3][3];

 /* Die Ecken sind gegeben durch (l1, l2, l3) (pti) = vi.
    Wichtig: v1, ..., v4 sind in der Reihenfolge des Umlaufs gelistet. */
 double      v1[3] = {1.0,  1.0,  1.0};
 double      v2[3] = {1.0, -1.0,  1.0};
 double      v3[3] = {1.0, -1.0, -1.0};
 /* double      v4[3] = {1.0,  1.0, -1.0}; */

 /* Inverse Matrix bis auf Skalierung. */
 adj[0][0] = lf[1][1]*lf[2][2] - lf[1][2]*lf[2][1];
 adj[0][1] = - (lf[0][1]*lf[2][2] - lf[0][2]*lf[2][1]);
 adj[0][2] = lf[0][1]*lf[1][2] - lf[0][2]*lf[1][1];

 adj[1][0] = - (lf[1][0]*lf[2][2] - lf[1][2]*lf[2][0]);
 adj[1][1] = lf[0][0]*lf[2][2] - lf[0][2]*lf[2][0];
 adj[1][2] = - (lf[0][0]*lf[1][2] - lf[0][2]*lf[1][0]);

 adj[2][0] = lf[1][0]*lf[2][1] - lf[1][1]*lf[2][0];
 adj[2][1] = - (lf[0][0]*lf[2][1] - lf[0][1]*lf[2][0]);
 adj[2][2] = lf[0][0]*lf[1][1] - lf[0][1]*lf[1][0];

 /* Determinante */
 det_inv = 1.0 / (lf[0][0]*adj[0][0] + lf[0][1]*adj[1][0] + lf[0][2]*adj[2][0]);
 /* sprintf (ausg, "det_inv = %.15f\n", det_inv); out (); */
 /* Inverse */
 for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
   inv[i][j] = adj[i][j] * det_inv;

 MMUL (pt1, inv, v1);
 MMUL (pt2, inv, v2);
 MMUL (pt3, inv, v3);
 /* MMUL (pt4, inv, v4); */
 pt4[0] = pt1[0] + pt3[0] - pt2[0];
 pt4[1] = pt1[1] + pt3[1] - pt2[1];
 pt4[2] = pt1[2] + pt3[2] - pt2[2];

 /* sprintf (ausg, "Ecken:\n"); out ();
 sprintf (ausg, "pt1 = (%f, %f, %f)\n", pt1[0], pt1[1], pt1[2]); out ();
 sprintf (ausg, "pt2 = (%f, %f, %f)\n", pt2[0], pt2[1], pt2[2]); out ();
 sprintf (ausg, "pt3 = (%f, %f, %f)\n", pt3[0], pt3[1], pt3[2]); out ();
 sprintf (ausg, "pt4 = (%f, %f, %f)\n", pt4[0], pt4[1], pt4[2]); out (); */
}


/* z runs in the outermost loop.
    So need maximum and minimum of the z-coordinates of the 5 corners
    the pyramid. */
inline void z_schranken (long *z_anf, long *z_end,
                         double *p1, double *p2, double *p3, double *p4) {
 double  tmp1, tmp2;

 tmp1 = MIN (MIN (MIN (p1[2], p2[2]), MIN (p3[2], p4[2])), 0.0);
 tmp2 = MAX (MAX (MAX (p1[2], p2[2]), MAX (p3[2], p4[2])), 0.0);

 *z_anf = lround (tmp1 + 0.5); /* Round up. */
 *z_end = lround (tmp2 - 0.5); /* round off. */
 /* Whole numbers may be rounded incorrectly.
  But this is not a bug, because we are anyway only the points inside the pyramid
  need. */

 /* sprintf (ausg, "z-Schranken: [%ld, %ld]\n", *z_anf, *z_end); out (); */
}


/* An edge represented as array of length 6.
   y = kante[2]*z + kante[3] und
   x = kante[4]*z + kante[5]
   for kante[0] <= z <= kante[1]. */
inline void kante (double *kante, double *anf, double *end) {
 double  dz_inv;

 dz_inv = 1 / (anf[2] - end[2]);

 if (anf[2] < end[2]) {
  kante[0] = anf[2];
  kante[1] = end[2];
 }
 else {
  kante[0] = end[2];
  kante[1] = anf[2];
 }
 /* If this is a division by 0, then edge will be [0] = edge [1]
  (and hopefully this is not an integer). */
 kante[2] = (anf[1] - end[1]) * dz_inv; /* dy/dz */
 kante[3] = anf[1] - kante[2]*anf[2];

 kante[4] = (anf[0] - end[0]) * dz_inv; /* dx/dz */
 kante[5] = anf[0] - kante[4]*anf[2];

 /* sprintf (ausg, "Kante: y = %f*z + %f, x = %f*z + %f fuer %f <= z <= %f\n",
               kante[2], kante[3], kante[4], kante[5], kante[0], kante[1]);
 out (); */
}


/* Calculate the eight edges of the pyramid. */
inline void pyr_kanten (double kant[8][6],
                 double *p1, double *p2, double *p3, double *p4) {
 double  s[3] = {0.0, 0.0, 0.0};

 kante (kant[0], s, p1);
 kante (kant[1], s, p2);
 kante (kant[2], s, p3);
 kante (kant[3], s, p4);
 kante (kant[4], p1, p2);
 kante (kant[5], p2, p3);
 kante (kant[6], p3, p4);
 kante (kant[7], p4, p1);
}


/* x and y run in the two inner loops.
    So we cut the pyramid with the z = ... - plane and have one
    Polygon. We need the maximum and the minimum of the x and y coordinates
    the corners of this polygon.
    Corners are created by cutting the plane with the edges. */
inline void x_und_y_schranken
                  (long *x_anf, long *x_end, long *y_anf, long *y_end,
                   double kant[8][6], long z) {
 int                          i;
 double                 *kan, k;
 double  xanf, xend, yanf, yend;

 xanf = LONG_MAX;
 xend = LONG_MIN;
 yanf = LONG_MAX;
 yend = LONG_MIN;

 /* Loop through the eight edges. */
 for (i = 0; i < 8; i++) {
  kan = kant[i]; /* Aktuelle Kante. */
  /* Do we meet this edge at all? */
  if ((kan[0] <= z) && (z <= kan[1])) {
   k = kan[2]*z + kan[3];
   yanf = MIN (yanf, k);
   yend = MAX (yend, k);
   k = kan[4]*z + kan[5];
   xanf = MIN (xanf, k);
   xend = MAX (xend, k);
  }
 }
 *y_anf = lround (yanf + 0.5); /* Aufrunden. */
 *y_end = lround (yend - 0.5); /* Abrunden. */
 *x_anf = lround (xanf + 0.5); /* Aufrunden. */
 *x_end = lround (xend - 0.5); /* Abrunden. */

 /* sprintf (ausg, "x_y-Schranken fuer z = %ld: [%ld, %ld] x [%ld, %ld]\n",
                                         z, *x_anf, *x_end, *y_anf, *y_end);
 out (); */
}


/* Issue of the solution. */
void post_proc (long vec_out0, long vec_out1, long vec_out2, long f0) {
 /* double  quot, x0, x1; */

 /* quot = ((double) vec_out0) / ((double) vec_out2); */
 /* sprintf (ausg, "%.18f\n", quot); out (); */

 /* x0 = mpx_get_d (x_0);
 x0 -= half_step;
 x1 = x0 + 2 * half_step; */
 /* sprintf (ausg, "Intervall [%.18f,%.18f]\n\n", x0, x1); out (); */

 /* if ((x0 < quot) && (x1 > quot)) */ {
  sprintf (ausg, "(%ld, %ld, %ld) Loesung fuer k = %ld.\n",
                                          vec_out0, vec_out1, vec_out2, f0);
  out ();
  sprintf (ausg, "\n"); out ();
  /* sprintf (ausg, "Solution found.\n"); out (); */
  /* exit (0); */
 }
}


/* What is v2 ** 3 - v0 ** 3 - v1 ** 3 at the current values of z, y and x?
    All invoices modulo 2 ** 64.
    Output of the triples (v0, v1, v2) where v0 ** 3 + v1 ** 3 - v2 ** 3 really
    is equal to \ pm3. */
inline ulong a_functional_value (long x, long v00, long v01, long v02,
                                long vec_h0, long vec_h1, long vec_h2) {
 long  vec_out0, vec_out1, vec_out2;
 ulong                           f0;

 /* sprintf (ausg, "x = %ld.\n", x); out (); */

 vec_out0 = x * v00 + vec_h0; /* Would be here bill in ulong logical? */
 vec_out1 = x * v01 + vec_h1; /* Have no overflow anyway. */
 vec_out2 = x * v02 + vec_h2;

 /* Bill modulo 2 ** 64. That's why Cast after unsigned. */
 f0 = ((ulong) vec_out2) * ((ulong) vec_out2) * ((ulong) vec_out2) -
      ((ulong) vec_out0) * ((ulong) vec_out0) * ((ulong) vec_out0) -
      ((ulong) vec_out1) * ((ulong) vec_out1) * ((ulong) vec_out1);

 /* Calculate in small_vect only the case z> = 0.
  That's why we need plus and minus at this point.
  Output of course as signed long ints. */
 if ((f0 < EDITION_BARRIER) || (-f0 < EDITION_BARRIER))
  if (f0 != 0)
   post_proc (vec_out0, vec_out1, vec_out2, f0); /* output */

 return (f0);
}


#define XY_SCHLEIFE {                                                       \
  for (y = y_anf; y <= y_end; y++) {                                        \
   /* sprintf (ausg, "Calculate y = %ld.\n", y); out (); */                    \
   vec_h0 = y * v[1][0] + z * v[2][0];                                      \
   vec_h1 = y * v[1][1] + z * v[2][1];                                      \
   vec_h2 = y * v[1][2] + z * v[2][2];                                      \
                                                                            \
   x = x_anf;                                                               \
   /* This function call causes an output if f0 or -f0 \
    is below the EDITION_BARRIER. */                                 \
   f0 = a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
   f1 = a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
   f2 = a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);       \
   anz++;                                                                   \
                                                                            \
   x++;                                                                     \
                                                                            \
   /* difference scheme.*/                                                  \
   d1 = f1 - f0;    d = f2 - f1;                                            \
   dd =  d - d1;                                                            \
   f = f2;                                                                  \
   dd += ddd; /* dd one step ahead. */                            \
   for (; x <= x_end; x++) {                                                \
    d  += dd;                                                               \
    dd += ddd;                                                              \
    f  += d;                                                                \
    /* if (f != a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2)) {\
     sprintf (ausg, "%lu %lu\n", f,                                         \
                 a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2));\
     out ();
     sprintf (ausg, "Bug im Differenzenschema!\n"); out ();                 \
     exit (0);                                                              \
    } */                                                                    \
    anz++;                                                                  \
    if ((f < EDITION_BARRIER) || (-f < EDITION_BARRIER))                  \
     /* Serves only the output. Calculate the function value superfluously
      again. Is faster than calling an extra function. */     \
     a_functional_value (x, v00, v01, v02, vec_h0, vec_h1, vec_h2);          \
   }                                                                        \
  }                                                                         \
} while (0);


/* Search for whole vectors with
         0 <l1 <1, | l2 | <l1 and | l3 | <l1.
    All bills take place with respect to the reduced base v.
    Should (vec [0], vec [1], vec [2]) consist of 3 divisible numbers,
    Let's not call post_proc.
    (Because of vec [0] ** 3 + vec [1] ** 3 + vec [2] ** 3 = 0 (mod 3) are either all
    three components divisible by 3 or none at all.) */
inline long kleine_vect (double lf[3][3], long v[3][3]) {
 long                       x, y, z;
 long                  x_anf, x_end;
 long                  y_anf, y_end;
 long                  z_anf, z_end;
 long        vec_h0, vec_h1, vec_h2;
 long                           anz;
 long                 v00, v01, v02;
 double  p1[3], p2[3], p3[3], p4[3];
 double                  kant[8][6];
 ulong                f, f0, f1, f2;
 ulong               d, dd, ddd, d1;

 /* sprintf (ausg, "lll-Basis:\n[%ld %ld %ld],\n[%ld %ld %ld],\n[%ld %ld %ld]\n",
           v[0][0], v[0][1], v[0][2],
           v[1][0], v[1][1], v[1][2],
           v[2][0], v[2][1], v[2][2]);
 out ();
 sprintf (ausg, "Linear shapes on the first vector: %f %f %f\n",
                                               lf[0][0], lf[1][0], lf[2][0]);
 out (); */

 /* Calculate the corners of the pyramid. The top is the zero point. */
 pyr_ecken (p1, p2, p3, p4, lf);
 /* Calculate the barriers for z. */
 z_schranken (&z_anf, &z_end, p1, p2, p3, p4);

 /* Calculate the eight edges of the pyramid. */
 pyr_kanten (kant, p1, p2, p3, p4);

 /* Triple loop through the pyramid. */
 anz = 0;
 v00 = v[0][0]; v01 = v[0][1]; v02 = v[0][2];
 /* Difference scheme, first part.
  ddd is constant regardless of z, over the entire tile. */
 ddd =  ((ulong) 6) *
       (((ulong) v02) * ((ulong) v02) * ((ulong) v02)
       -((ulong) v00) * ((ulong) v00) * ((ulong) v00)
       -((ulong) v01) * ((ulong) v01) * ((ulong) v01));

 for (z = z_anf; z <= z_end; z++) {
  /* Compute the bounds for x and y. */
  x_und_y_schranken (&x_anf, &x_end, &y_anf, &y_end, kant, z);
  XY_SCHLEIFE;
 }
 return (anz);
}


/* Init. The outer loop. */
void calculate_interval (mpx_t x_0_anf, mpx_t x_0_ende) {
 mpx_t         increment;
 mpx_t                  y_0;
 double  lf[3][3], ln[3][3];
 long      e[3][3], v[3][3];
 long     counter, anz, err;

 x_0[0] = x_0_ende[0]; x_0[1] = x_0_ende[1];
 /* Run backwards from x_0_end to x_0_anf. */

 calculate_tile (increment);

 /* First loop pass. Here we calculate LLL with a lot of precision.
     In addition, y_0 is initialized.*/
 init (v, y_0, x_0,
           half_step, tiles_offset, half_tile_width, seeking_wide);
 y_inv_init (y_0);
 y_diff[0] = 0; y_diff[1] = 0; y_inv_diff[0] = 0; y_inv_diff[1] = 0;
 /* y_diff = y_inv_diff = 0. */
 sprintf (ausg, "Init fertig.\n"); out ();

 /* Loop. Here is enough for almost everything the accuracy of double. */
 counter = NEW_TILE; tile = -NEW_TILE; err = 0;
 /* while (x_0 >= x_0_anf) ... . */
 while ((x_0[1] > x_0_anf[1]) ||
       ((x_0[1] == x_0_anf[1]) && (x_0[0] >= x_0_anf[0]))) {
  /* mpf_set_mpx (tmp1, x_0);
  gmp_sprintf (ausg, "\nx_0 = %.*Ff.\n", 40, tmp1); out (); */

  berechne_y_wert (y_0);
  /* The first time y_diff = 0 => y_0 is calculated correctly. */
  /* mpf_set_mpx (tmp1, y_0);
  gmp_sprintf (ausg, "y_0 = %.*Ff.\n\n", 40, tmp1); out (); */
  berechne_y_strich (y_0);
  /* The first time y_inv_diff = 0 => x'(x_0) is calculated correctly. */
  /* mpf_set_mpx (tmp1, Ax);
  gmp_sprintf (ausg, "A = %.*Ff.\n\n", 40, tmp1); out (); */

  if (counter == NEW_TILE) {
   counter = 0; tile += NEW_TILE;
   calculate_tile (increment);
   y_inv_init (y_0); y_diff_init (y_0, increment);
  }

  /* Last v was not correct due to overflow. */
  if (err > 0) {
   init (v, y_0, x_0,
           half_step, tiles_offset, half_tile_width, seeking_wide);
   mpf_set_mpx (tmp1, x_0);
   gmp_sprintf (ausg, "Restart at x_0 = %.*Ff.\n", 25, tmp1); out ();
  }

  calculate_three_linearf (lf, v, y_0, half_tile_width, seeking_wide);

  lll (lf, e);
  err = matrix_prod (v, e, v); /* v is definitely modulo 2**64. */
 
  lf_neu (ln, lf, e);
  anz = kleine_vect (ln, v);
  /* sprintf (ausg, "Find %ld grid points.\n", anz); out (); */

  /* tile go from x_0 to right and left. */
  mpx_sub (x_0, x_0, increment); counter++;
 }

 sprintf (ausg, "Insgesamt %ld tile behandelt.\n", tile + counter);
 out ();
}


/* We want with something like
       elkies_allg 0.4 0.000001
    start.

    Where 0.4 is the initial value.
    0.000001 is the interval length that the processor should create. */
int main (int argc, char *argv[]) {
 mpx_t   diff, x_0_anf, x_0_ende;

 // mpf_set_default_prec: Set the default precision to be at least prec bits. 
 // All subsequent calls to mpf_init will use this precision, but previously 
 // initialized variables are unaffected.
 mpf_set_default_prec (128);
 mpf_init (tmp1); mpf_init (tmp2);

 mpf_set_str (tmp1, argv[1], 10); mpx_set_mpf (x_0_anf, tmp1);
 mpf_set_str (tmp2, argv[2], 10); mpx_set_mpf (diff, tmp2);

 row = 1;
 mpx_add (x_0_ende, x_0_anf, diff);

 /* Name of the output file */
 mpf_set_mpx (tmp1, x_0_anf); mpf_set_mpx (tmp2, x_0_ende);
 gmp_sprintf (file, "liste_allg_%.*Ff_%.*Ff.txt", 8, tmp1, 8, tmp2);

 /* First edition */
 gmp_sprintf (ausg, "Start  from %.*Ff bis %.*Ff.\n",
                                                      15, tmp1, 15, tmp2);
 out ();

 calculate_interval (x_0_anf, x_0_ende);

 mpf_clear (tmp1); mpf_clear (tmp2);
 return 0;
}

