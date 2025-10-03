#include <fips202.h>
#include <tutil.h>
#include <mp.h>
#include <encoded_sizes.h>
#include <ec_params.h>
#include <verification.h>
#include <assert.h>

void
public_key_init(public_key_t *pk)
{
    ec_curve_init(&pk->curve);
}

void
hash_to_challenge(scalar_mtx_2x2_t *chall,
                  const public_key_t *pk,
                  const unsigned char *message,
                  int message_length,
                  const ec_curve_t *commitment,
                  int torsion_exponent)
{
    unsigned char buf[2 * FP2_ENCODED_BYTES];

    size_t hash_bytes = (4 * torsion_exponent + 7) / 8;
    size_t limbs = (hash_bytes + sizeof(digit_t) - 1) / sizeof(digit_t);
    size_t bits = (4 * torsion_exponent) % RADIX;
    digit_t mask = ((digit_t)-1) >> ((RADIX - bits) % RADIX);
    scalar_t scalar;
    memset(scalar, 0, NWORDS_ORDER * sizeof(digit_t));
    {
        fp2_t j1, j2;
        ec_j_inv(&j1, &pk->curve);
        ec_j_inv(&j2, commitment);
        fp2_encode(buf, &j1);
        fp2_encode(buf + FP2_ENCODED_BYTES, &j2);
    }
    int ok = 0;
    uint32_t counter = 0;
    while (!ok) {
        counter = counter + 1;
        // hash C||counter
        shake256incctx ctx;
#ifdef TARGET_BIG_ENDIAN
        mask = BSWAP_DIGIT(mask);
#endif

        shake256_inc_init(&ctx);
        shake256_inc_absorb(&ctx, buf, 2 * FP2_ENCODED_BYTES);
        shake256_inc_absorb(&ctx, (uint8_t *)&counter, 4);
        shake256_inc_absorb(&ctx, message, message_length);
        shake256_inc_finalize(&ctx);
        shake256_inc_squeeze((void *)(scalar), hash_bytes, &ctx);
        (scalar)[limbs - 1] &= mask;
#ifdef TARGET_BIG_ENDIAN
        mask = BSWAP_DIGIT(mask);
        for (int i = 0; i < NWORDS_ORDER; i++)
            (*scalar)[i] = BSWAP_DIGIT((*scalar)[i]);
#endif

        {
            // parse output as 4 ints in [0,torsion]
            scalar_t mask;
            memset(mask, 0, NWORDS_ORDER * sizeof(digit_t));
            mask[0] = 1;
            mp_neg(mask, NWORDS_ORDER);
            // mp_shiftr only works for shifts up to RADIX - 1
            for (int shift = sizeof(scalar_t) * 8 - torsion_exponent; shift > 0; shift -= RADIX - 1)
                mp_shiftr(mask, shift >= RADIX ? RADIX - 1 : shift, NWORDS_ORDER);

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < NWORDS_ORDER; k++) {
                        (*chall)[i][j][k] = mask[k] & scalar[k];
                    }
                    // mp_shiftr only works for shifts up to RADIX - 1
                    for (int shift = torsion_exponent; shift > 0; shift -= RADIX - 1)
                        mp_shiftr(scalar, shift >= RADIX ? RADIX - 1 : shift, NWORDS_ORDER);
                }
            }
        }

        // check parity of det(M_chall) which must be odd
        ok = (mp_is_odd((*chall)[0][0], NWORDS_ORDER) & mp_is_odd((*chall)[1][1], NWORDS_ORDER)) !=
             (mp_is_odd((*chall)[0][1], NWORDS_ORDER) & mp_is_odd((*chall)[1][0], NWORDS_ORDER));
    }
}
