#include <verification.h>
#include <signature.h>
#include <encoded_sizes.h>
#include <mp.h>

int
verify_test_hash_to_challenge()
{
    int res = 0;
    scalar_mtx_2x2_t chall;
    unsigned char m[7] = "Hello!";

    commitment_curve_t commitment;
    public_key_t pk;
    secret_key_t sk;
    secret_commitment_t sk_com;
    ibz_mat_2x2_t chall_ibz;
    ibz_t det;

    ibz_mat_2x2_init(&chall_ibz);
    ibz_init(&det);
    commitment_curve_init(&commitment);
    public_key_init(&pk);
    secret_key_init(&sk);
    secret_commitment_init(&sk_com);

    protocols_keygen(&pk, &sk);
    protocols_commitment_gen(&commitment, &sk_com);

    hash_to_challenge(&chall, &pk, m, 7, &commitment, TORSION_challenge_torsion);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_copy_digit_array(&(chall_ibz[i][j]), chall[i][j]);
        }
    }

    ibz_mat_2x2_det(&det, &chall_ibz);
    ibz_gcd(&det, &det, &ibz_const_two);

    res = res || !ibz_is_one(&det);

    if (res) {
        printf("Verification test hash_to_challenge failed\n");
    }
    commitment_curve_finalize(&commitment);
    public_key_finalize(&pk);
    secret_key_finalize(&sk);
    secret_commitment_finalize(&sk_com);
    ibz_mat_2x2_finalize(&chall_ibz);
    ibz_finalize(&det);
    return (res);
}

int
verify_test_scalar_mtx_2x2_inv_mod_2exp_up_to_scalar()
{
    int res = 0;
    scalar_mtx_2x2_t a, inv;
    ibz_mat_2x2_t x, inv_x, prod;
    ibz_t mod;
    int exp = 10;
    ibz_mat_2x2_init(&x);
    ibz_mat_2x2_init(&inv_x);
    ibz_mat_2x2_init(&prod);
    ibz_init(&mod);
    ibz_pow(&mod, &ibz_const_two, exp);

    // create matrix
    ibz_mat_2x2_set(&x, 1, 5, 3, 0);
    ibz_mat_2x2_to_digit_array(&a, &x);
    mp_neg(a[1][0], sizeof(scalar_t) / sizeof(digit_t));
    ibz_neg(&(x[1][0]), &(x[1][0]));
    // apply inv
    scalar_mtx_2x2_inv_mod_2exp_up_to_scalar(&inv, &a, exp, sizeof(scalar_t) / sizeof(digit_t));

    // translate (back) to ibz
    ibz_mat_2x2_from_digit_array(&inv_x, &inv);
    // compare to multiple of id...
    ibz_mat_2x2_mul_mod(&prod, &x, &inv_x, &mod);
    res = res || !ibz_is_zero(&(prod[0][1]));
    res = res || !ibz_is_zero(&(prod[1][0]));
    res = res || !(ibz_cmp(&(prod[1][1]), &(prod[0][0])) == 0);
    res = res || !ibz_is_odd(&(prod[1][1]));
    ibz_mat_2x2_mul_mod(&prod, &inv_x, &x, &mod);
    res = res || !ibz_is_zero(&(prod[0][1]));
    res = res || !ibz_is_zero(&(prod[1][0]));
    res = res || !(ibz_cmp(&(prod[1][1]), &(prod[0][0])) == 0);
    res = res || !ibz_is_odd(&(prod[1][1]));

    if (res) {
        printf("Verification test scalar_mtx_2x2_inv_mod_2exp_up_to_scalar failed\n");
    }
    ibz_mat_2x2_finalize(&x);
    ibz_mat_2x2_finalize(&inv_x);
    ibz_mat_2x2_finalize(&prod);
    ibz_finalize(&mod);
    return (res);
}

int
verify_test_scalar_mtx_2x2_mul_mod_2exp()
{
    int res = 0;
    scalar_mtx_2x2_t a, b, inv;
    ibz_mat_2x2_t x, y, r, prod;
    ibz_t mod;
    int exp = 10;
    ibz_mat_2x2_init(&x);
    ibz_mat_2x2_init(&y);
    ibz_mat_2x2_init(&r);
    ibz_mat_2x2_init(&prod);
    ibz_init(&mod);
    ibz_pow(&mod, &ibz_const_two, exp);

    ibz_mat_2x2_set(&x, 1, 5, 3, 0);
    ibz_mat_2x2_to_digit_array(&a, &x);
    mp_neg(a[1][0], sizeof(scalar_t) / sizeof(digit_t));
    ibz_neg(&(x[1][0]), &(x[1][0]));

    ibz_mat_2x2_set(&y, 2, 5, 7, 0);
    ibz_mat_2x2_to_digit_array(&b, &y);
    mp_neg(b[1][0], sizeof(scalar_t) / sizeof(digit_t));
    ibz_neg(&(y[1][0]), &(y[1][0]));

    scalar_mtx_2x2_mul_mod_2exp(&inv, &a, &b, exp, sizeof(scalar_t) / sizeof(digit_t));
    // translate (back) to ibz

    ibz_mat_2x2_from_digit_array(&r, &inv);
    // compare to multiple of id...
    ibz_mat_2x2_mul_mod(&prod, &x, &y, &mod);
    res = res || !(ibz_cmp(&(prod[0][0]), &(r[0][0])) == 0);
    res = res || !(ibz_cmp(&(prod[1][0]), &(r[1][0])) == 0);
    res = res || !(ibz_cmp(&(prod[0][1]), &(r[0][1])) == 0);
    res = res || !(ibz_cmp(&(prod[1][1]), &(r[1][1])) == 0);

    if (res) {
        printf("Verification test scalar_mtx_2x2_mul_mod_2exp failed\n");
    }
    ibz_mat_2x2_finalize(&x);
    ibz_mat_2x2_finalize(&y);
    ibz_mat_2x2_finalize(&r);
    ibz_mat_2x2_finalize(&prod);
    ibz_finalize(&mod);
    return (res);
}

int
main()
{
    printf("Running verification tests:\n\n");
    int res = 0;
    res = res | verify_test_hash_to_challenge();
    res = res | verify_test_scalar_mtx_2x2_inv_mod_2exp_up_to_scalar();
    res = res | verify_test_scalar_mtx_2x2_mul_mod_2exp();
    if (res) {
        printf("\nSome tests failed!\n");
    }
    return (res);
}