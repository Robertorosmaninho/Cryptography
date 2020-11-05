#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_SIZE 32

//MDC Extended using Euler
void mdc_estendido(mpz_t g, mpz_t x, mpz_t y, const mpz_t a, const mpz_t b) {
  mpz_t q, r, a1, b1;
  mpz_inits(q, r, a1, b1, NULL);
  mpz_t qx1, qy1;
  mpz_inits(qx1, qy1, NULL);

  mpz_set(g, b1);
  mpz_set_ui(x, 0);
  mpz_set_ui(y, 1);

  mpz_set(a1, a);
  mpz_set(b1, b);

  mpz_t x0, y0, x1, y1;
  mpz_inits(x0, y0, x1, y1, NULL);

  mpz_set_ui(x0, 1);
  mpz_set_ui(y0, 0);
  mpz_set_ui(x1, 0);
  mpz_set_ui(y1, 1);

  mpz_tdiv_qr(q, r, a1, b1);
  //x2 = x0 - q*x1
  //y2 = y0 - q*y1
  //x0 = x1
  //x1 = x2
  while (mpz_cmp_ui(r, 0) != 0) {

    mpz_set(a1, b1);
    mpz_set(b1, r);

    mpz_mul(qx1, q, x1);
    mpz_sub(x, x0, qx1);

    mpz_mul(qy1, q, y1);
    mpz_sub(y, y0, qy1);

    mpz_set(g, r);

    mpz_set(x0, x1);
    mpz_set(x1, x);

    mpz_set(y0, y1);
    mpz_set(y1, y);

    mpz_tdiv_qr(q, r, a1, b1);
  }

  mpz_clears(q, r, x0, x1, y0, y1, qx1, qy1, NULL);
}

//Calculating the inverse mod
int inverso_modular(mpz_t r, const mpz_t a, const mpz_t n) {
  mpz_t a1, n1, g, x, y;
  mpz_inits(a1, n1, g, x, y, NULL);

  mpz_set(a1, a);
  mpz_set(n1, n);

  mdc_estendido(g, x, y, a1, n1);

  if (!mpz_cmp_ui(g, 1)) {
    mpz_set(r, x);
    mpz_clears(g, x, y, NULL);
    return 1;
  }

  mpz_clears(g, x, y, NULL);
  return 0;
}

void exp_binaria(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t n) {

  mpz_t base, exp;
  mpz_inits(base, exp, NULL);
  mpz_set(base, b);
  mpz_set(exp, e);
  mpz_set_ui(r, 1);

  mpz_tdiv_r(base, base, n);

  while (mpz_cmp_ui(exp, 0) > 0) {

    //Teste de paridade
    if (mpz_odd_p(exp) != 0) {
      mpz_mul(r, r, base);
      mpz_tdiv_r(r, r, n);
    }

    //Chão de exp/2->quadrado da base->mod

    mpz_fdiv_q_ui(exp, exp, 2);
    mpz_mul(base, base, base);
    mpz_tdiv_r(base, base, n);
  }

  mpz_clears(base, exp, NULL);
}

int talvez_primo(const mpz_t a, const mpz_t n, const mpz_t n1,
                 unsigned int t, const mpz_t q) {

  //Testa se o numero é composto e retorna 0 case seja ou 1 caso seja primo
  //ou pseudoprimo forte para a base a.

  mpz_t r;
  mpz_init(r);

  unsigned int i = 0;

  exp_binaria(r, a, q, n);

  while (i < t) {
    if ((i == 0 && mpz_cmp_ui(r, 1) == 0) || (i >= 0 && mpz_cmp(r, n1) == 0)) {
      mpz_clear(r);
      return 1;
    } else {
      i++;
      mpz_mul(r, r, r);
      mpz_tdiv_r(r, r, n);

      if (i >= t) {
        mpz_clear(r);
        return 0;
      }
    }
  }
  return 0;//Se o teste falhar completamente.
}

//Função que gera número aleatorio
void numero_aleatorio(mpz_t r, const mpz_t n1, gmp_randstate_t rnd) {
  mp_bitcnt_t num_bits = mpz_sizeinbase(n1, 2);
  do {
    mpz_urandomb(r, rnd, num_bits);
  } while (!(mpz_cmp_ui(r, 2) >= 0 && mpz_cmp(r, n1) <= 0));
}

// Miller-Rabin's Test
int provavelmente_primo(const mpz_t n, unsigned int iter,
                        gmp_randstate_t rnd) {
  int its_prime;

  mpz_t n1, a;
  mpz_inits(n1, a, NULL);

  mpz_set(n1, n);
  mpz_sub_ui(n1, n1, 1);

  //Creating q and t to use in maybe_prime
  unsigned int t = 0;
  mpz_t aux, q;
  mpz_inits(aux, q, NULL);
  mpz_set(aux, n1);

  while (mpz_even_p(aux) != 0) {
    t++;
    mpz_tdiv_q_ui(aux, aux, 2);
    mpz_set(q, aux);
  }

  for (int i = 0; i < iter; i++) {
    numero_aleatorio(a, n, rnd);
    its_prime = talvez_primo(a, n, n1, t, q);

    if (!its_prime) {
      mpz_clears(n1, a, aux, q, NULL);
      return 0;//If isn't prime
    }
  }

  mpz_clears(n1, a, aux, q, NULL);
  return 1;//If n is probably prime
}


//Creating random primes
void primo_aleatorio(mpz_t r, unsigned int b, gmp_randstate_t rnd) {

  int prob_prime = 0;

  while (!prob_prime) {
    mpz_urandomb(r, rnd, b);
    if (mpz_cmp_ui(r, 2) >= 0) {
      prob_prime = provavelmente_primo(r, 20, rnd);
    }
  }
}

//Generating Keys
void gera_chaves(mpz_t n, mpz_t e, mpz_t d, gmp_randstate_t rnd) {

  mpz_t p, q, phi_n, p1, q1;
  mpz_inits(p, q, phi_n, p1, q1, NULL);

  //Generating p and q
  primo_aleatorio(p, 2048, rnd);//mudar para 2048
  primo_aleatorio(q, 2048, rnd);

  //Generating n
  mpz_mul(n, p, q);

  //Generating phi_n
  mpz_sub_ui(p1, p, 1);
  mpz_sub_ui(q1, q, 1);
  mpz_mul(phi_n, p1, q1);

  //Generating e and d, then  testing if n, e and d it's a valid key
  mpz_set_ui(e, 65536);
  mpz_set_ui(d, 1);
  mpz_t g, aux;
  mpz_inits(g, aux, NULL);
  do {
    mpz_add_ui(e, e, 1);
    if (inverso_modular(d, e, phi_n)) {
      break;
    }
    mdc_estendido(g, d, aux, e, phi_n);
  } while (mpz_cmp_ui(g, 1) != 0);

  //	inverso_modular(d,e,phi_n);

  while (mpz_cmp_ui(d, 0) <= 0) {
    mpz_add(d, d, phi_n);
  }

  mpz_clears(p, q, phi_n, p1, q1, NULL);
}

// Transforms Strings to numbers in ASCII
void codifica(mpz_t r, const char *str) {
  mpz_t aux, mult;
  mpz_inits(aux, mult, NULL);

  mpz_set_ui(r, 0);

  for (int i = 0; i < strlen(str); i++) {
    mpz_ui_pow_ui(aux, 256, i);
    mpz_mul_ui(mult, aux, (int) str[i]);
    mpz_add(r, r, mult);
  }
  mpz_clears(aux, mult, NULL);
}

char *decodifica(const mpz_t n) {
  int i = 0;
  char *str = (char *) malloc(500);

  mpz_t aux, mult, q, r, cod;
  mpz_inits(aux, mult, q, r, cod, NULL);

  mpz_set(cod, n);
  do {
    mpz_tdiv_qr_ui(q, r, cod, 256);
    mpz_set(cod, q);

    str[i++] = (char) mpz_get_ui(r);
  } while (mpz_cmp_ui(q, 0) != 0);

  mpz_clears(aux, mult, q, r, cod, NULL);
  return str;
}

void criptografa(mpz_t C, const mpz_t M, const mpz_t n, const mpz_t e) {
  exp_binaria(C, M, e, n);
}

void descriptografa(mpz_t M, const mpz_t C, const mpz_t n, const mpz_t d) {
  exp_binaria(M, C, d, n);
}
