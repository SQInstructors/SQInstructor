#include "intbig.h"
#include <signature.h>
#include <tools.h>
#include <quaternion_data.h>
#include <id2iso.h>
#include <torsion_constants.h>
#include <encoded_sizes.h>
#include <mp.h>
#include <endomorphism_action.h>

void
ibz_mat_2x2_to_digit_array(scalar_mtx_2x2_t *a, const ibz_mat_2x2_t *x)
{

    ibz_to_digit_array((*a)[0][0], &((*x)[0][0]));
    ibz_to_digit_array((*a)[1][0], &((*x)[1][0]));
    ibz_to_digit_array((*a)[0][1], &((*x)[0][1]));
    ibz_to_digit_array((*a)[1][1], &((*x)[1][1]));
}

void
ibz_mat_2x2_from_digit_array(ibz_mat_2x2_t *a, const scalar_mtx_2x2_t *x)
{

    ibz_copy_digit_array(&((*a)[0][0]), (*x)[0][0]);
    ibz_copy_digit_array(&((*a)[1][0]), (*x)[1][0]);
    ibz_copy_digit_array(&((*a)[0][1]), (*x)[0][1]);
    ibz_copy_digit_array(&((*a)[1][1]), (*x)[1][1]);
}

static void
compute_and_set_basis_change_matrix(ibz_mat_2x2_t *change,
                                    const ec_basis_t *B_aux_mapped,
                                    const ec_basis_t *B_aux_can,
                                    ec_curve_t *E_aux,
                                    int f)
{
#ifndef NDEBUG
    {
        // Ensure all points have the desired order
        assert(test_basis_order_twof(B_aux_can, E_aux, TORSION_EVEN_POWER));
        assert(test_basis_order_twof(B_aux_mapped, E_aux, f));
        fp2_t w0;
        weil(&w0, f, &B_aux_mapped->P, &B_aux_mapped->Q, &B_aux_mapped->PmQ, E_aux);
    }
#endif

    // compute the matrix to go from B_aux_mapped to B_aux_can
    change_of_basis_matrix_tate(change, B_aux_mapped, B_aux_can, E_aux, f);
    ibz_mat_2x2_transpose(change);

    // Assert all values in the matrix are of the expected size for packing
    assert(ibz_bitsize(&(*change)[0][0]) <= SQIsign_response_length + HD_extra_torsion);
    assert(ibz_bitsize(&(*change)[0][1]) <= SQIsign_response_length + HD_extra_torsion);
    assert(ibz_bitsize(&(*change)[1][0]) <= SQIsign_response_length + HD_extra_torsion);
    assert(ibz_bitsize(&(*change)[1][1]) <= SQIsign_response_length + HD_extra_torsion);
}

// if result is 0, no solution was found
// NOTE: do not call this function with I_solutions
// equal to another input
int
compute_solutions(quat_left_ideal_t *I_solutions,
                  quat_lattice_t *I_solutions_left_order,
                  quat_lattice_t *lattice_solutions,
                  const ibz_mat_2x2_t *chall,
                  const secret_commitment_t *commitment,
                  const secret_key_t *sk)
{
    // I_solutions = dual(I_sk) * I_com
    {
        quat_left_ideal_t dual_I_sk_ideal;
        quat_left_ideal_init(&dual_I_sk_ideal);
        quat_lideal_conjugate_without_hnf(
            &dual_I_sk_ideal, I_solutions_left_order, &(sk->secret_ideal), &QUATALG_PINFTY);
        quat_lideal_lideal_mul(I_solutions, &dual_I_sk_ideal, &(commitment->secret_ideal), &QUATALG_PINFTY);
        quat_left_ideal_finalize(&dual_I_sk_ideal);
    }

    // Compute the target matrix
    //   M_com · adj(M_chal) · adj(M_sk)
    ibz_mat_2x2_t M_target;
    ibz_t torsion;
    ibz_mat_2x2_init(&M_target);
    ibz_init(&torsion);
    ibz_pow(&torsion, &ibz_const_two, TORSION_challenge_torsion);

    {
        ibz_mat_2x2_t M_chall_adj;
        ibz_mat_2x2_init(&M_chall_adj);
        ibz_mat_2x2_adjugate(&M_chall_adj, chall);
        ibz_mat_2x2_adjugate(&M_target, &sk->secret_action);
        ibz_mat_2x2_mul_mod(&M_target, &M_chall_adj, &M_target, &torsion);
        ibz_mat_2x2_mul_mod(&M_target, &commitment->secret_action, &M_target, &torsion);
        ibz_mat_2x2_finalize(&M_chall_adj);
    }

    // Solve linear system
    quat_alg_elem_t i_zero;
    quat_alg_elem_init(&i_zero);
    {
        // 1. Put M_target in a column vector
        // Attention: precomupted matrices ACTION_X are transposed!
        ibz_vec_4_t target;
        ibz_vec_4_init(&target);
        ibz_copy(&target[0], &M_target[0][0]);
        ibz_copy(&target[1], &M_target[1][0]);
        ibz_copy(&target[2], &M_target[0][1]);
        ibz_copy(&target[3], &M_target[1][1]);
        // 2. Put the action of the basis of B_p,∞ in a matrix
        ibz_mat_4x4_t O_action;
        ibz_mat_4x4_init(&O_action);
        ibz_set(&O_action[0][0], 1);
        ibz_set(&O_action[3][0], 1);
        ibz_copy(&O_action[0][1], &ACTION_I[0][0]);
        ibz_copy(&O_action[1][1], &ACTION_I[0][1]);
        ibz_copy(&O_action[2][1], &ACTION_I[1][0]);
        ibz_copy(&O_action[3][1], &ACTION_I[1][1]);
        ibz_copy(&O_action[0][2], &ACTION_J[0][0]);
        ibz_copy(&O_action[1][2], &ACTION_J[0][1]);
        ibz_copy(&O_action[2][2], &ACTION_J[1][0]);
        ibz_copy(&O_action[3][2], &ACTION_J[1][1]);
        ibz_copy(&O_action[0][3], &ACTION_K[0][0]);
        ibz_copy(&O_action[1][3], &ACTION_K[0][1]);
        ibz_copy(&O_action[2][3], &ACTION_K[1][0]);
        ibz_copy(&O_action[3][3], &ACTION_K[1][1]);
        // 3. Multiply by the lattice basis to get the action of its
        // endomorphisms
        ibz_mat_4x4_mul(&O_action, &O_action, &I_solutions->lattice.basis);
        int divisible = ibz_mat_4x4_scalar_div(&O_action, &I_solutions->lattice.denom, &O_action);
        assert(divisible);
        (void)divisible;
        ibz_mat_4x4_mod(&O_action, &O_action, &torsion);
        // 4. Solve system
        ibz_mat_4x4_solve_system_mod(&target, &O_action, &target, &torsion);
        // 5. Convert solution to standard basis
        quat_lattice_to_standard_basis(&i_zero, &(I_solutions->lattice), &target);

        ibz_mat_4x4_finalize(&O_action);
        ibz_vec_4_finalize(&target);
    }

    // Replace I_solutions with equivalent ideal
    // TODO: kinda broken when i_zero ∈ ℤ
    if (!ibz_is_zero(&i_zero.coord[1]) || !ibz_is_zero(&i_zero.coord[2]) || !ibz_is_zero(&i_zero.coord[3]))
        quat_lideal_equivalent(I_solutions, I_solutions, &i_zero, &QUATALG_PINFTY);

    // Compute the lattice of equivalent solutions:
    //   I ∩ O(E, S) = N(I)Z + Ntor*I
    ibz_mat_4x4_scalar_mul(&(lattice_solutions->basis), &torsion, &(I_solutions->lattice.basis));
    ibz_copy(&(lattice_solutions->denom), &(I_solutions->lattice.denom));
    quat_lattice_hnf(lattice_solutions);
    assert(ibz_is_zero(&lattice_solutions->basis[1][0]));
    assert(ibz_is_zero(&lattice_solutions->basis[2][0]));
    assert(ibz_is_zero(&lattice_solutions->basis[3][0]));
    ibz_mul(&(lattice_solutions->basis[0][0]), &(I_solutions->norm), &(lattice_solutions->denom));
    quat_lattice_reduce_denom(lattice_solutions, lattice_solutions);
    // 11: Λ ← I ∩ O(E, S); seems useless
    // quat_lattice_intersect(lattice_solutions, &(I_solutions->lattice), lattice_solutions);
    ibz_mat_2x2_finalize(&M_target);
    quat_alg_elem_finalize(&i_zero);
    ibz_finalize(&torsion);
    return (1);
}

void
set_radius(ibz_t *radius, const ibz_t *Ntors, const quat_alg_t *alg)
{
    ibz_pow(radius, &ibz_const_two, SQIsign_response_length);
    ibz_sub(radius, radius, &ibz_const_one);
#ifndef NDEBUG
    // Radius >= sqrt(pN³)
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_pow(&tmp, Ntors, 3);
    ibz_mul(&tmp, &tmp, &alg->p);
    ibz_sqrt_floor(&tmp, &tmp);
    assert(ibz_cmp(radius, &tmp) >= 0);
    ibz_finalize(&tmp);
#endif
}

void
hash_to_chall_wrapper(ibz_mat_2x2_t *chall,
                      const public_key_t *pk,
                      const unsigned char *message,
                      int message_length,
                      ec_curve_t *commitment)
{
    scalar_mtx_2x2_t chall_internal;
    hash_to_challenge(&chall_internal, pk, message, message_length, commitment, TORSION_challenge_torsion);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_copy_digit_array(&((*chall)[i][j]), chall_internal[i][j]);
        }
    }
}

int
protocols_sign(signature_t *sig, const public_key_t *pk, secret_key_t *sk, const unsigned char *m, size_t l)
{
    int res = 1;
    secret_commitment_t com;
    ibz_mat_2x2_t chall, change;
    quat_left_ideal_t I_solutions;
    quat_lattice_t I_solutions_left_order;
    quat_left_ideal_t I_resp;
    quat_left_ideal_t I_aux;
    ibz_t radius, norm, tmp, torsion;
    quat_alg_elem_t elem;
    ibz_vec_4_t elem_coords;
    quat_lattice_t lattice_solutions;
    ec_basis_t B_mapped;
    ec_curve_t E_aux;
    ec_curve_init(&E_aux);
    quat_lattice_init(&lattice_solutions);
    quat_alg_elem_init(&elem);
    ibz_vec_4_init(&elem_coords);
    ibz_init(&radius);
    ibz_init(&norm);
    ibz_init(&tmp);
    ibz_init(&torsion);
    quat_left_ideal_init(&I_resp);
    quat_left_ideal_init(&I_solutions);
    quat_left_ideal_init(&I_aux);
    quat_lattice_init(&I_solutions_left_order);
    secret_commitment_init(&com);
    ibz_mat_2x2_init(&chall);
    ibz_mat_2x2_init(&change);
    ibz_pow(&torsion, &ibz_const_two, TORSION_challenge_torsion);

    // Commitment generation
    // Stores the commitment hint in the signature
    res = protocols_commitment_gen(sig, &com);
    if (!res)
        goto fin;

    // Hash to Challenge
    hash_to_chall_wrapper(&chall, pk, m, l, &com.curve);

    // this computes a solution I_solutions, then the lattice_solutions that can be
    // used to sample any solution ideal
    res = res && compute_solutions(&I_solutions, &I_solutions_left_order, &lattice_solutions, &chall, &com, sk);
    if (!res)
        goto fin;

    // Here we perform a Lattice-Ball intersection on Lambda with suitable radius
    // from this we get the response ideal I_resp
    set_radius(&radius, &torsion, &QUATALG_PINFTY);
    // adjust radius by ideal norm
    ibz_mul(&radius, &radius, &I_solutions.norm);
    for (int i = 0;; i++) {
        res = res && quat_lattice_sample_from_ball(&elem, &elem_coords, &lattice_solutions, &QUATALG_PINFTY, &radius);
        if (!res)
            goto fin;
        quat_alg_norm(&norm, &tmp, &elem, &QUATALG_PINFTY);
        assert(ibz_is_one(&tmp));

        if (ibz_is_odd(&norm))
            break;
        if (i >= 256) {
            debug_print("Too many rejections, restarting signing.");
            res = 0;
            goto fin;
        }
    }
    // get I_resp from the solutions element
    quat_lideal_equivalent(&I_resp, &I_solutions, &elem, &QUATALG_PINFTY);

    // Compute the scalar matrix corresponding to the action of I_resp
    ibz_invmod(&norm, &com.secret_ideal.norm, &torsion);
    ibz_mat_2x2_det(&tmp, &com.secret_action);
    ibz_mod(&tmp, &tmp, &torsion);
    ibz_mul(&tmp, &tmp, &norm);
    ibz_mod(&tmp, &tmp, &torsion);
    ibz_mul(&tmp, &tmp, &elem_coords[0]);
    ibz_mod(&tmp, &tmp, &torsion);
    ibz_to_digit_array(sig->scalar, &tmp);
    
    // 8. Id2iso on Iresp
    //  create an O0-ideal I_aux of norm &norm = 2**a - q
    //  with q = norm(I_resp) and a = resp_len = ceil(log_2(q))
    //  ISSUE this should be uint16_t
    sig->resp_length = (uint16_t)ibz_bitsize(&(I_resp.norm));
    assert(sig->resp_length <= SQIsign_response_length);
    ibz_pow(&norm, &ibz_const_two, sig->resp_length);
    ibz_sub(&norm, &norm, &(I_resp.norm));
    assert(ibz_cmp(&norm, &ibz_const_zero) > 0);
    assert(ibz_cmp(&norm, &QUATALG_PINFTY.p) < 0);
    res = res && quat_sampling_random_ideal_O0_given_norm(
                     &(I_aux), &norm, 0, &QUAT_represent_integer_params, &QUAT_prime_cofactor);
    if (!res)
        goto fin;

    // save norm for later
    ibz_copy(&norm, &I_resp.norm);

    // Push forward I_aux through I_sk*I_resp, then compute I_sk*I_resp*I_aux
    quat_lideal_lideal_mul(&I_resp, &sk->secret_ideal, &I_resp, &QUATALG_PINFTY);
    quat_lattice_intersect(&I_aux.lattice, &I_aux.lattice, &I_resp.lattice);
    quat_lideal_inverse_lattice_without_hnf(&I_solutions.lattice, &(I_resp), &QUATALG_PINFTY);
    quat_lattice_mul(&I_aux.lattice, &I_solutions.lattice, &I_aux.lattice, &QUATALG_PINFTY);
    quat_lattice_mul(&I_resp.lattice, &I_resp.lattice, &I_aux.lattice, &QUATALG_PINFTY);
    ibz_mul(&I_resp.norm, &I_aux.norm, &I_resp.norm);

    // Translate final I_resp to isogeny
    res = res && dim2id2iso_arbitrary_isogeny_evaluation(&B_mapped, &E_aux, &I_resp);
    if (!res)
        goto fin;

    // normalize  Montgomery A-coefficient and store
    ec_normalize_curve_and_A24(&E_aux);
    fp2_copy(&(sig->E_aux_A), &(E_aux.A));

    // multiply B_mapped by cofactor 2**(e - sig->resp_length - HD_extra_torsion)
    ec_dbl_iter_basis(&B_mapped, TORSION_EVEN_POWER - sig->resp_length - HD_extra_torsion, &B_mapped, &E_aux);

    // get the hint and torsion basis for E_aux
    ec_basis_t B_aux_can;
    sig->hint_aux = ec_curve_to_basis_2f_to_hint(&B_aux_can, &E_aux, TORSION_EVEN_POWER);

    // compute change of matrix to get from B_aux_can the image B_mapped
    compute_and_set_basis_change_matrix(&change, &B_mapped, &B_aux_can, &E_aux, sig->resp_length + HD_extra_torsion);

    // Adjust the matrix to get the image of B_pk from B_aux_can
    ibz_pow(&torsion, &ibz_const_two, sig->resp_length + HD_extra_torsion);
    ibz_mul(&chall[0][0], &norm, &sk->secret_action[0][0]);
    ibz_mul(&chall[0][1], &norm, &sk->secret_action[0][1]);
    ibz_mul(&chall[1][0], &norm, &sk->secret_action[1][0]);
    ibz_mul(&chall[1][1], &norm, &sk->secret_action[1][1]);
    ibz_mat_2x2_inv_mod(&chall, &chall, &torsion);
    ibz_mat_2x2_mul_mod(&change, &chall, &change, &torsion);
    ibz_mat_2x2_to_digit_array(&(sig->mat_B_aux_can_to_B_aux), &change);

// Finalizations
fin:;
    quat_lattice_finalize(&lattice_solutions);
    ibz_finalize(&radius);
    ibz_finalize(&norm);
    ibz_finalize(&tmp);
    ibz_finalize(&torsion);
    ibz_vec_4_finalize(&elem_coords);
    quat_alg_elem_finalize(&elem);
    ibz_mat_2x2_finalize(&chall);
    ibz_mat_2x2_finalize(&change);
    quat_left_ideal_finalize(&I_solutions);
    quat_lattice_finalize(&I_solutions_left_order);
    quat_left_ideal_finalize(&I_resp);
    quat_left_ideal_finalize(&I_aux);
    secret_commitment_finalize(&com);
    return (res);
}
