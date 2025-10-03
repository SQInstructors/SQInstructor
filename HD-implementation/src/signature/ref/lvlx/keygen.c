#include <signature.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>
#include <intbig.h>
#include <id2iso.h>
#include <torsion_constants.h>
#include <endomorphism_action.h>
#include <encoded_sizes.h>

void
secret_key_init(secret_key_t *sk)
{
    secret_commitment_init(sk);
}

void
secret_key_finalize(secret_key_t *sk)
{
    secret_commitment_finalize(sk);
}

void
secret_commitment_init(secret_commitment_t *sk)
{
    quat_left_ideal_init(&(sk->secret_ideal));
    ibz_mat_2x2_init(&(sk->secret_action));
    ec_curve_init(&sk->curve);
}

void
secret_commitment_finalize(secret_commitment_t *sk)
{
    quat_left_ideal_finalize(&(sk->secret_ideal));
    ibz_mat_2x2_finalize(&(sk->secret_action));
}

// assumes torsion to be a power of 2
// computes M such that aM = b
void
change_of_basis_from_lattices(quat_lattice_t *res, const quat_lattice_t *a, const quat_lattice_t *b)
{
    quat_lattice_t change_of_basis;
    quat_lattice_init(&change_of_basis);
    // get inverse of a's matrix
    ibz_mat_4x4_inv_with_det_as_denom(&(change_of_basis.basis), &(change_of_basis.denom), &(a->basis));
    // adjust for denominator of order
    ibz_mat_4x4_scalar_mul(&(change_of_basis.basis), (&(a->denom)), &(change_of_basis.basis));
    // multiply by basis of b
    ibz_mat_4x4_mul(&(change_of_basis.basis), &(change_of_basis.basis), &(b->basis));
    // Adjust for denominator
    ibz_mul(&(change_of_basis.denom), &(change_of_basis.denom), &(b->denom));
    // Reduce denominator of change of basis matrix
    quat_lattice_reduce_denom(res, &change_of_basis);

    quat_lattice_finalize(&change_of_basis);

}

// res = basis*change mod torsion
void
apply_change_of_basis_mod_torsion(ibz_mat_4x4_t *res,
                                  const quat_lattice_t *change_of_basis,
                                  const ibz_mat_4x4_t *basis,
                                  const ibz_t *torsion)
{
    ibz_t factor;
    ibz_init(&factor);
    ibz_mat_4x4_mul(res, basis, &(change_of_basis->basis));
    ibz_invmod(&factor, &(change_of_basis->denom), torsion);
    ibz_mat_4x4_scalar_mul(res, &factor, res);
    ibz_mat_4x4_mod(res, res, torsion);
    ibz_finalize(&factor);
}

uint16_t
common_part_of_gen(secret_common_t *sk)
{
    uint16_t found = 0;
    ec_basis_t B_0_two;

    // iterating until a solution has been found
    while (!found) {

        found = quat_sampling_random_ideal_O0_given_norm(
            &sk->secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);

        // replacing the secret ideal by a shorter equivalent one for efficiency
        found = found && quat_lideal_prime_norm_reduced_equivalent(
                             &sk->secret_ideal, &QUATALG_PINFTY, QUAT_primality_num_iter, QUAT_equiv_bound_coeff);

        // ideal to isogeny clapotis

        found = found && dim2id2iso_arbitrary_isogeny_evaluation(&B_0_two, &sk->image_curve, &sk->secret_ideal);
    }

    // Assert the isogeny was found and images have the correct order
    assert(test_basis_order_twof(&B_0_two, &sk->image_curve, TORSION_EVEN_POWER));

    // Normalize the curve model
    ec_curve_standard(&sk->curve, &sk->isom, &sk->image_curve);
    ec_iso_eval(&B_0_two.P, &sk->isom);
    ec_iso_eval(&B_0_two.Q, &sk->isom);
    ec_iso_eval(&B_0_two.PmQ, &sk->isom);
    
    // Compute a deterministic basis with a hint to speed up verification
    uint8_t hint = ec_curve_to_basis_2f_to_hint(&sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER);

    ec_dbl_iter_basis(
        &sk->canonical_basis_tor, TORSION_EVEN_POWER - TORSION_challenge_torsion, &sk->canonical_basis, &sk->curve);

    // Assert the deterministic basis we computed has the correct order
    assert(test_basis_order_twof(&sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER));
    assert(test_basis_order_twof(&sk->canonical_basis_tor, &sk->curve, TORSION_challenge_torsion));

    // Compute the action of the secret isogeny on the basis
    change_of_basis_matrix_tate(
        &sk->secret_action, &B_0_two, &sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER);
    ibz_mat_2x2_transpose(&sk->secret_action);
    
#ifndef NDEBUG
    ec_point_t P;
    digit_t a[NWORDS_ORDER], b[NWORDS_ORDER];
    ibz_to_digit_array(a, &sk->secret_action[0][0]);
    ibz_to_digit_array(b, &sk->secret_action[0][1]);
    ec_biscalar_mul(&P, a, b,
                    TORSION_EVEN_POWER,
		    &sk->canonical_basis,
		    &sk->curve);
    assert(ec_is_equal(&P, &B_0_two.P));
    
    ibz_to_digit_array(a, &sk->secret_action[1][0]);
    ibz_to_digit_array(b, &sk->secret_action[1][1]);
    ec_biscalar_mul(&P, a, b,
                    TORSION_EVEN_POWER,
		    &sk->canonical_basis,
		    &sk->curve);
    assert(ec_is_equal(&P, &B_0_two.Q));
#endif

    return ((found << 8) | ((uint16_t)hint));
}

int
protocols_keygen(public_key_t *pk, secret_key_t *sk)
{
    int found = 0;

    uint16_t out = common_part_of_gen(sk);
    found = (int)(out >> 8);
    pk->hint_pk = (uint8_t)(out);
    copy_curve(&pk->curve, &sk->curve);

    assert(fp2_is_one(&pk->curve.C) == 0xFFFFFFFF);

    return found;
}

int
protocols_commitment_gen(signature_t *sig, secret_commitment_t *sk_cmt)
{
    int found = 0;

    uint16_t out = common_part_of_gen(sk_cmt);
    found = (int)(out >> 8);

    if (found)
      sig->hint_com = out & 0xff;

    return found;
}
