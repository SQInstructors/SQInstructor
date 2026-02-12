#include <verification.h>
#include <mp.h>
#include <hd.h>
#include <biextension.h>
#include <encoded_sizes.h>
#include <assert.h>

// Check that the basis change matrix elements are canonical
// representatives modulo 2^(SQIsign_response_length + 2).
static int
check_canonical_basis_change_matrix(const signature_t *sig)
{
    // This works as long as all values in sig->mat_B_aux_can_to_B_aux are
    // positive integers.
    // GIACOMO: adapted to new sig struct
    int ret = 1;
    scalar_t aux;

    memset(aux, 0, NWORDS_ORDER * sizeof(digit_t));
    aux[0] = 0x1;
    multiple_mp_shiftl(aux, sig->resp_length, NWORDS_ORDER);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (mp_compare(aux, sig->mat_B_aux_can_to_B_aux[i][j], NWORDS_ORDER) <= 0) {
                ret = 0;
            }
        }
    }

    return ret;
}

// same as matrix_application_even_basis() in id2iso.c, with some modifications:
// - THE MATRIX IS TRANSPOSED!!!!
// - this version works with a matrix of scalars (not ibz_t).
// - reduction modulo 2^f of matrix elements is removed here, because it is
//   assumed that the elements are already canonical representatives modulo
//   2^f; this is ensured by calling check_canonical_basis_change_matrix() at
//   the beginning of protocols_verify().
static int
matrix_scalar_application_even_basis(ec_basis_t *bas, const ec_curve_t *E, scalar_mtx_2x2_t *mat, int f)
{
    scalar_t scalar0, scalar1;
    memset(scalar0, 0, NWORDS_ORDER * sizeof(digit_t));
    memset(scalar1, 0, NWORDS_ORDER * sizeof(digit_t));

    ec_basis_t tmp_bas;
    copy_basis(&tmp_bas, bas);

    // For a matrix [[a, b], [c, d]] we compute:
    //
    // first basis element R = [a]P + [b]Q
    if (!ec_biscalar_mul(&bas->P, (*mat)[0][0], (*mat)[0][1], f, &tmp_bas, E))
        return 0;
    // second basis element S = [c]P + [d]Q
    if (!ec_biscalar_mul(&bas->Q, (*mat)[1][0], (*mat)[1][1], f, &tmp_bas, E))
        return 0;
    // Their difference R - S = [a - c]P + [b - d]Q
    mp_sub(scalar0, (*mat)[0][0], (*mat)[1][0], NWORDS_ORDER);
    mp_mod_2exp(scalar0, f, NWORDS_ORDER);
    mp_sub(scalar1, (*mat)[0][1], (*mat)[1][1], NWORDS_ORDER);
    mp_mod_2exp(scalar1, f, NWORDS_ORDER);
    return ec_biscalar_mul(&bas->PmQ, scalar0, scalar1, f, &tmp_bas, E);
}

/*
// Compute the bases for the challenge and auxillary curve from
// the canonical bases. Challenge basis is reconstructed from the
// compressed scalars within the challenge.
static int
challenge_and_aux_basis_verify(ec_basis_t *B_chall_can,
                               ec_basis_t *B_aux_can,
                               ec_curve_t *E_chall,
                               ec_curve_t *E_aux,
                               signature_t *sig,
                               const int pow_dim2_deg_resp)
{

    // recovering the canonical basis as TORSION_EVEN_POWER for consistency with signing
    if (!ec_curve_to_basis_2f_from_hint(B_chall_can, E_chall, TORSION_EVEN_POWER, sig->hint_chall))
        return 0;

    // setting to the right order
    ec_dbl_iter_basis(B_chall_can,
                      TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion - sig->two_resp_length,
                      B_chall_can,
                      E_chall);

    if (!ec_curve_to_basis_2f_from_hint(B_aux_can, E_aux, TORSION_EVEN_POWER, sig->hint_aux))
        return 0;

    // setting to the right order
    ec_dbl_iter_basis(B_aux_can, TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion, B_aux_can, E_aux);

#ifndef NDEBUG
    if (!test_basis_order_twof(B_chall_can, E_chall, HD_extra_torsion + pow_dim2_deg_resp + sig->two_resp_length))
        debug_print("canonical basis has wrong order, expect something to fail");
#endif

    // applying the change matrix on the basis of E_chall
    return matrix_scalar_application_even_basis(B_chall_can,
                                                E_chall,
                                                &sig->mat_Bchall_can_to_B_chall,
                                                pow_dim2_deg_resp + HD_extra_torsion + sig->two_resp_length);
}

// When two_resp_length is non-zero, we must compute a small 2^n-isogeny
// updating E_chall as the codomain as well as push the basis on E_chall
// through this isogeny
static int
two_response_isogeny_verify(ec_curve_t *E_chall, ec_basis_t *B_chall_can, const signature_t *sig, int pow_dim2_deg_resp)
{
    ec_point_t ker, points[3];

    // choosing the right point for the small two_isogenies
    if (mp_is_even(sig->mat_Bchall_can_to_B_chall[0][0], NWORDS_ORDER) &&
        mp_is_even(sig->mat_Bchall_can_to_B_chall[1][0], NWORDS_ORDER)) {
        copy_point(&ker, &B_chall_can->Q);
    } else {
        copy_point(&ker, &B_chall_can->P);
    }

    copy_point(&points[0], &B_chall_can->P);
    copy_point(&points[1], &B_chall_can->Q);
    copy_point(&points[2], &B_chall_can->PmQ);

    ec_dbl_iter(&ker, pow_dim2_deg_resp + HD_extra_torsion, &ker, E_chall);

#ifndef NDEBUG
    if (!test_point_order_twof(&ker, E_chall, sig->two_resp_length))
        debug_print("kernel does not have order 2^(two_resp_length");
#endif

    if (ec_eval_small_chain(E_chall, &ker, sig->two_resp_length, points, 3, false)) {
        return 0;
    }

#ifndef NDEBUG
    if (!test_point_order_twof(&points[0], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[0] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
    if (!test_point_order_twof(&points[1], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[1] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
    if (!test_point_order_twof(&points[2], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[2] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
#endif

    copy_point(&B_chall_can->P, &points[0]);
    copy_point(&B_chall_can->Q, &points[1]);
    copy_point(&B_chall_can->PmQ, &points[2]);
    return 1;
}
    */

// The commitment curve can be recovered from the codomain of the 2D
// isogeny built from the bases computed during verification.
static int
compute_commitment_curve_verify(ec_curve_t *E_com,
                                ec_basis_t *E_com_basis_image,
                                const ec_basis_t *E_pk_basis_to_be_pushed,
                                const ec_basis_t *B_pk,
                                const ec_basis_t *B_aux,
                                const ec_curve_t *E_pk,
                                const ec_curve_t *E_aux,
                                int pow_dim2_deg_resp)

{
#ifndef NDEBUG
    // Check all the points are the correct order
    if (!test_basis_order_twof(B_pk, E_pk, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("B_pk does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");

    if (!test_basis_order_twof(B_aux, E_aux, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("B_aux does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
#endif

    // now compute the dim2 isogeny from E_pk x E_aux -> E_com x E_aux'
    // of kernel B_pk x B_aux

    // first we set-up the kernel
    theta_couple_curve_t EpkxEaux;
    copy_curve(&EpkxEaux.E1, E_pk);
    copy_curve(&EpkxEaux.E2, E_aux);

    theta_kernel_couple_points_t dim_two_ker;
    copy_bases_to_kernel(&dim_two_ker, B_pk, B_aux);
    theta_couple_point_t B_pk_theta[3];
    copy_point(&B_pk_theta[0].P1, &E_pk_basis_to_be_pushed->P);
    copy_point(&B_pk_theta[1].P1, &E_pk_basis_to_be_pushed->Q);
    copy_point(&B_pk_theta[2].P1, &E_pk_basis_to_be_pushed->PmQ);
    ec_point_init(&B_pk_theta[0].P2);
    ec_point_init(&B_pk_theta[1].P2);
    ec_point_init(&B_pk_theta[2].P2);

    // computing the isogeny
    theta_couple_curve_t codomain;
    int codomain_splits;
    ec_curve_init(&codomain.E1);
    ec_curve_init(&codomain.E2);
    // handling the special case where we don't need to perform any dim2 computation
    if (pow_dim2_deg_resp == 0) {
        codomain_splits = 1;
        copy_curve(&codomain.E1, &EpkxEaux.E1);
        copy_curve(&codomain.E2, &EpkxEaux.E2);
        // We still need to check that E_chall is supersingular
        // This assumes that HD_extra_torsion == 2
        if (!ec_is_basis_four_torsion(B_pk, E_pk)) {
            return 0;
        }
    } else {
        codomain_splits = theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp, &EpkxEaux, &dim_two_ker, true, &codomain, B_pk_theta, 3);
    }

    // computing the commitment curve
    // its always the first one because of our (2^n,2^n)-isogeny formulae
    copy_curve(E_com, &codomain.E1);

    fp2_copy(&E_com_basis_image->P.x, &B_pk_theta[0].P1.x);
    fp2_copy(&E_com_basis_image->P.z, &B_pk_theta[0].P1.z);
    fp2_copy(&E_com_basis_image->Q.x, &B_pk_theta[1].P1.x);
    fp2_copy(&E_com_basis_image->Q.z, &B_pk_theta[1].P1.z);
    fp2_copy(&E_com_basis_image->PmQ.x, &B_pk_theta[2].P1.x);
    fp2_copy(&E_com_basis_image->PmQ.z, &B_pk_theta[2].P1.z);

    return codomain_splits;
}
/*
// SQIsign verification
int
protocols_verify(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l)
{
    int verify;

    if (!check_canonical_basis_change_matrix(sig))
        return 0;

    // Computation of the length of the dim 2 2^n isogeny
    int pow_dim2_deg_resp = SQIsign_response_length - (int)sig->two_resp_length - (int)sig->backtracking;

    // basic sanity test: checking that the response is not too long
    if (pow_dim2_deg_resp < 0)
        return 0;
    // The dim 2 isogeny embeds a dim 1 isogeny of odd degree, so it can
    // never be of length 2.
    if (pow_dim2_deg_resp == 1)
        return 0;

    // check the public curve is valid
    if (!ec_curve_verify_A(&(pk->curve).A))
        return 0;

    // Set auxiliary curve from the A-coefficient within the signature
    ec_curve_t E_aux;
    if (!ec_curve_init_from_A(&E_aux, &sig->E_aux_A))
        return 0; // invalid curve

    // checking that we are given A-coefficients and no precomputation
    assert(fp2_is_one(&pk->curve.C) == 0xFFFFFFFF && !pk->curve.is_A24_computed_and_normalized);

    // computation of the challenge
    ec_curve_t E_chall;
    if (!compute_challenge_verify(&E_chall, sig, &pk->curve, pk->hint_pk)) {
        return 0;
    }

    // Computation of the canonical bases for the challenge and aux curve
    ec_basis_t B_chall_can, B_aux_can;

    if (!challenge_and_aux_basis_verify(&B_chall_can, &B_aux_can, &E_chall, &E_aux, sig, pow_dim2_deg_resp)) {
        return 0;
    }

    // When two_resp_length != 0 we need to compute a second, short 2^r-isogeny
    if (sig->two_resp_length > 0) {
        if (!two_response_isogeny_verify(&E_chall, &B_chall_can, sig, pow_dim2_deg_resp)) {
            return 0;
        }
    }

    // We can recover the commitment curve with a 2D isogeny
    // The supplied signature did not compute an isogeny between eliptic products
    // and so definitely is an invalid signature.
    ec_curve_t E_com;
    if (!compute_commitment_curve_verify(&E_com, &B_chall_can, &B_aux_can, &E_chall, &E_aux, pow_dim2_deg_resp))
        return 0;

    scalar_t chk_chall;

    // recomputing the challenge vector
    hash_to_challenge(&chk_chall, pk, &E_com, m, l);

    // performing the final check
    verify = mp_compare(sig->chall_coeff, chk_chall, NWORDS_ORDER) == 0;

    return verify;
}
*/

void
scalar_mtx_2x2_mod_2exp(scalar_mtx_2x2_t *x, int exp, int nwords)
{
    mp_mod_2exp((*x)[0][0], exp, nwords);
    mp_mod_2exp((*x)[0][1], exp, nwords);
    mp_mod_2exp((*x)[1][0], exp, nwords);
    mp_mod_2exp((*x)[1][1], exp, nwords);
}

// return 1 if invertible, 0 else
int
scalar_mtx_2x2_inv_mod_2exp_up_to_scalar(scalar_mtx_2x2_t *inv, const scalar_mtx_2x2_t *x, int exp, int nwords)
{
    scalar_t det, prod;
    memset(&det, 0, nwords * sizeof(digit_t));
    memset(&prod, 0, nwords * sizeof(digit_t));
    mp_sub(det, prod, det, nwords);
    mp_mod_2exp(det, exp, nwords);
    mp_copy((*inv)[0][0], (*x)[1][1], nwords);
    mp_copy((*inv)[1][1], (*x)[0][0], nwords);
    mp_copy((*inv)[1][0], (*x)[1][0], nwords);
    mp_copy((*inv)[0][1], (*x)[0][1], nwords);
    mp_neg((*inv)[0][1], nwords);
    mp_neg((*inv)[1][0], nwords);
    return (mp_is_zero(det, nwords));
}

void
scalar_mtx_2x2_mul_mod_2exp(scalar_mtx_2x2_t *prod,
                            const scalar_mtx_2x2_t *a,
                            const scalar_mtx_2x2_t *b,
                            int exp,
                            int nwords)
{
    scalar_mtx_2x2_t work;
    scalar_t tmp, temp;
    memset(tmp, 0, nwords * sizeof(digit_t));
    memset(temp, 0, nwords * sizeof(digit_t));
    memset(work, 0, 4 * nwords * sizeof(digit_t));
    mp_mul(tmp, (*a)[0][0], (*b)[0][0], nwords);
    mp_mul(temp, (*a)[0][1], (*b)[1][0], nwords);
    mp_add(work[0][0], temp, tmp, nwords);
    mp_mul(tmp, (*a)[1][0], (*b)[0][1], nwords);
    mp_mul(temp, (*a)[1][1], (*b)[1][1], nwords);
    mp_add(work[1][1], temp, tmp, nwords);

    mp_mul(tmp, (*a)[1][0], (*b)[0][0], nwords);
    mp_mul(temp, (*a)[1][1], (*b)[1][0], nwords);
    mp_add(work[1][0], temp, tmp, nwords);
    mp_mul(tmp, (*a)[0][0], (*b)[0][1], nwords);
    mp_mul(temp, (*a)[0][1], (*b)[1][1], nwords);
    mp_add(work[0][1], temp, tmp, nwords);

    mp_copy((*prod)[0][0], work[0][0], nwords);
    mp_copy((*prod)[0][1], work[0][1], nwords);
    mp_copy((*prod)[1][0], work[1][0], nwords);
    mp_copy((*prod)[1][1], work[1][1], nwords);

    mp_mod_2exp((*prod)[0][0], exp, nwords);
    mp_mod_2exp((*prod)[0][1], exp, nwords);
    mp_mod_2exp((*prod)[1][0], exp, nwords);
    mp_mod_2exp((*prod)[1][1], exp, nwords);
}

int
protocols_verify(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l)
{
    int verify = 1;
    ec_basis_t B_pk_can, B_aux, B_com, B_com_pushed;
    ec_curve_t E_aux, E_pk, E_com;
    ec_isom_t isom_com;
    scalar_mtx_2x2_t M_chall, inv, prod;
    scalar_t test;
    ec_curve_init(&E_com);
    ec_curve_init(&E_aux);
    memset(&M_chall, 0, 4 * NWORDS_ORDER * sizeof(digit_t));
    memset(&inv, 0, 4 * NWORDS_ORDER * sizeof(digit_t));
    memset(&prod, 0, 4 * NWORDS_ORDER * sizeof(digit_t));
    memset(&test, 0, NWORDS_ORDER * sizeof(digit_t));
    ec_curve_init(&E_pk);
    copy_point(&E_pk.A24, &pk->curve.A24);
    fp2_copy(&E_pk.A, &pk->curve.A);
    fp2_copy(&E_pk.C, &pk->curve.C);
    E_pk.is_A24_computed_and_normalized = pk->curve.is_A24_computed_and_normalized;

    // 0. get the canonical basis from pk
    if (!ec_curve_to_basis_2f_from_hint(&B_pk_can, &E_pk, TORSION_EVEN_POWER, pk->hint_pk)) {
        debug_print("Bad pk hint");
        return 0;
    }
    // assert(ec_is_equal(&E_pk.A24, &pk->curve.A24));

    // 2. get the higher dimensional isogeny from the signature

    // 2.a get the canonical basis for the aux curve
    ec_curve_init_from_A(&E_aux, &sig->E_aux_A);
    if (!ec_curve_to_basis_2f_from_hint(&B_aux, &E_aux, TORSION_EVEN_POWER, sig->hint_aux)) {
        debug_print("Bad aux hint");
        return 0;
    }

    // Multiply B_pk_can and B_aux by cofactor
    ec_dbl_iter_basis(&B_pk_can, TORSION_EVEN_POWER - sig->resp_length - HD_extra_torsion, &B_pk_can, &E_pk);
    ec_dbl_iter_basis(&B_aux, TORSION_EVEN_POWER - sig->resp_length - HD_extra_torsion, &B_aux, &E_aux);

    // Adjust B_aux by the transformation matrix
    if (!matrix_scalar_application_even_basis(
            &B_aux, &E_aux, &sig->mat_B_aux_can_to_B_aux, sig->resp_length + HD_extra_torsion)) {
        debug_print("Cannot apply matrix to aux");
        return 0;
    }

    // Prepare basis to be pushed through isogeny
    ec_dbl_iter_basis(&B_com_pushed, sig->resp_length + HD_extra_torsion - TORSION_challenge_torsion, &B_pk_can, &E_pk);

    // 2.b get the isogeny and the codomain curve
    if (!compute_commitment_curve_verify(
            &E_com, &B_com_pushed, &B_com_pushed, &B_pk_can, &B_aux, &E_pk, &E_aux, sig->resp_length)) {
        debug_print("Cannot compute HD isogeny");
        return 0;
    }

    // standardize commitment
    ec_curve_standard(&E_com, &isom_com, &E_com);
    ec_iso_eval(&B_com_pushed.P, &isom_com);
    ec_iso_eval(&B_com_pushed.Q, &isom_com);
    ec_iso_eval(&B_com_pushed.PmQ, &isom_com);

    // Get canonical basis on E_com
    if (!ec_curve_to_basis_2f_from_hint(&B_com, &E_com, TORSION_EVEN_POWER, sig->hint_com)) {
        debug_print("Bad com hint");
        return 0;
    }
    ec_dbl_iter_basis(&B_com, TORSION_EVEN_POWER - TORSION_challenge_torsion, &B_com, &E_com);

    // 5. hash to get the challenge M_chall
    hash_to_challenge(&M_chall, pk, m, l, &E_com, TORSION_challenge_torsion);

    // 6. Apply M_chall to B_com
    if (!matrix_scalar_application_even_basis(&B_com, &E_com, &M_chall, TORSION_challenge_torsion)) {
        debug_print("Cannot apply challenge matrix");
        return 0;
    }

    // Multiply B_com by signature's scalar
    ec_mul(&B_com.P, sig->scalar, TORSION_challenge_torsion, &B_com.P, &E_com);
    ec_mul(&B_com.Q, sig->scalar, TORSION_challenge_torsion, &B_com.Q, &E_com);
    ec_mul(&B_com.PmQ, sig->scalar, TORSION_challenge_torsion, &B_com.PmQ, &E_com);

    // Compare B_com_pushed and B_com: they must be equal
    if (!ec_is_equal(&B_com.P, &B_com_pushed.P)
        || !ec_is_equal(&B_com.Q, &B_com_pushed.Q)
        || !ec_is_equal(&B_com.PmQ, &B_com_pushed.PmQ)) {
        debug_print("Response does not match challenge");
        return 0;
    }

    return (verify);
}
