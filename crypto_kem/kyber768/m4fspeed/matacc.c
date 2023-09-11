#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "ntt.h"
#include "matacc.h"


/*************************************************
* Name:        matacc_cache32
*
* Description: Multiplies a row of A or A^T, generated on-the-fly,
*              with a vector of polynomials and accumulates into the result.
*              Using asymmetric multiplication and better accumulation.
*
* Arguments:   - poly *r:                    pointer to output polynomial to accumulate in
*              - const polyvec *b:           pointer to input vector of polynomials to multiply with
*              - polyvec *b_prime:           pointer to output vector of polynomials to store b multiplied by zetas
*              - unsigned char i:            byte to indicate the index < KYBER_K of the row of A or A^T
*              - const unsigned char *seed:  pointer to the public seed used to generate A
*              - int transposed:             boolean indicatin whether A or A^T is generated
**************************************************/
void matacc_cache32(poly* r, const polyvec *b, polyvec *b_prime, unsigned char i, const unsigned char *seed, int transposed) {
  unsigned char buf[XOF_BLOCKBYTES+2];
  xof_state state;
  int16_t c[4];
  int32_t r_tmp[KYBER_N]; // stores intermediate accumulated values to save reductions
  int j = 0;

  // 16-32

  if (transposed)
    xof_absorb(&state, seed, i, j);
  else
    xof_absorb(&state, seed, j, i);

  xof_squeezeblocks(buf, 1, &state);

  // 첫번째 누적의 경우, 16-bit * 16-bit를 수행한 후, 32-bit 결과값을 그대로 누적함
  matacc_asm_cache_16_32(r_tmp, b->vec[j].coeffs, c, buf, zetas, &state, b_prime->vec[j].coeffs);

  // 32-32 KYBER_K - 2 times
  for(j=1;j<KYBER_K - 1;j++) {

    if (transposed)
      xof_absorb(&state, seed, i, j);
    else
      xof_absorb(&state, seed, j, i);

    xof_squeezeblocks(buf, 1, &state);
    
    // 두번째부터 마지막을 제외한 다항식 곱셈에 대해 16-bit * 16-bit 곱셈을 수행한 후, 기존 값이 누적되어있는 32-bit 배열에 누적하여 저장
    matacc_asm_cache_32_32(r_tmp, b->vec[j].coeffs, c, buf, zetas, &state, b_prime->vec[j].coeffs);
  }

  // 32-16

  if (transposed)
    xof_absorb(&state, seed, i, j);
  else
    xof_absorb(&state, seed, j, i);

  xof_squeezeblocks(buf, 1, &state);

  // 마지막 연산의 경우, 16-bit*16-bit 곱셈을 수행한 후, 32-bit 배열에 누적한 후 최종적으로 감산을 수행하여 16-bit 배열에 저장함
  matacc_asm_cache_32_16(r->coeffs, b->vec[j].coeffs, c, buf, zetas, &state, b_prime->vec[j].coeffs, r_tmp);
}

/*************************************************
* Name:        matacc_opt32
*
* Description: Multiplies a row of A or A^T, generated on-the-fly,
*              with a vector of polynomials and accumulates into the result.
*              Using asymmetric multiplication and better accumulation.
*
* Arguments:   - poly *r:                    pointer to output polynomial to accumulate in
*              - const polyvec *b:           pointer to input vector of polynomials to multiply with
*              - const polyvec *b_prime:     pointer to input vector of polynomials to store b multiplied by zetas
*              - unsigned char i:            byte to indicate the index < KYBER_K of the row of A or A^T
*              - const unsigned char *seed:  pointer to the public seed used to generate A
*              - int transposed:             boolean indicatin whether A or A^T is generated
**************************************************/
void matacc_opt32(poly* r, const polyvec *b, const polyvec *b_prime, unsigned char i, const unsigned char *seed, int transposed) {
  unsigned char buf[XOF_BLOCKBYTES+2];
  xof_state state;
  int16_t c[4];
  int32_t r_tmp[KYBER_N]; // stores intermediate accumulated values to save reductions
  int j = 0;

  // 16-32

  if (transposed)
    xof_absorb(&state, seed, i, j);
  else
    xof_absorb(&state, seed, j, i);

  xof_squeezeblocks(buf, 1, &state);

  matacc_asm_opt_16_32(r_tmp, b->vec[j].coeffs, c, buf, &state, b_prime->vec[j].coeffs);
  
  // 32-32 KYBER_K - 2 times
  for(j=1;j<KYBER_K - 1;j++) {

    if (transposed)
      xof_absorb(&state, seed, i, j);
    else
      xof_absorb(&state, seed, j, i);

    xof_squeezeblocks(buf, 1, &state);

    matacc_asm_opt_32_32(r_tmp, b->vec[j].coeffs, c, buf, &state, b_prime->vec[j].coeffs);
  }

  // 32-16

  if (transposed)
    xof_absorb(&state, seed, i, j);
  else
    xof_absorb(&state, seed, j, i);

  xof_squeezeblocks(buf, 1, &state);

  matacc_asm_opt_32_16(r->coeffs, b->vec[j].coeffs, c, buf, &state, b_prime->vec[j].coeffs, r_tmp);
}
