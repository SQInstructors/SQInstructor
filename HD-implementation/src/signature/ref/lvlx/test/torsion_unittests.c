#include <signature.h>
#include <id2iso.h>
#include <encoded_sizes.h>

int
signing_torsion_unittest_change_of_basis_from_lattices()
{
    int res = 0;

    quat_lattice_t a, b, change, prod;
    quat_lattice_init(&a);
    quat_lattice_init(&b);
    quat_lattice_init(&change);
    quat_lattice_init(&prod);

    ibz_set(&(a.basis[0][0]), 1);
    ibz_set(&(a.basis[1][1]), 1);
    ibz_set(&(a.basis[2][2]), 1);
    ibz_set(&(a.basis[3][3]), 1);
    ibz_set(&(a.basis[0][2]), 17);
    ibz_set(&(a.denom), 5);

    ibz_mat_4x4_copy(&b.basis, &a.basis);
    ibz_mat_4x4_scalar_mul(&b.basis, &(a.denom), &b.basis);
    ibz_set(&(b.basis[0][2]), 0);
    ibz_set(&(b.basis[0][1]), 11);
    ibz_set(&(b.basis[1][1]), -7);
    ibz_set(&(b.basis[1][2]), 19);
    ibz_set(&(b.denom), 3);

    change_of_basis_from_lattices(&change, &a, &b);
    ibz_mat_4x4_mul(&prod.basis, &a.basis, &change.basis);
    ibz_mul(&prod.denom, &a.denom, &change.denom);
    res = res || !quat_lattice_equal(&prod, &b);

    if (res) {
        printf("Torsion signing unit test change_of_basis_from_lattices failed\n");
    }
    quat_lattice_finalize(&a);
    quat_lattice_finalize(&b);
    quat_lattice_finalize(&change);
    quat_lattice_finalize(&prod);
    return (res);
}

int
signing_torsion_unittest_apply_change_of_basis_mod_torsion()
{
    int res = 0;

    quat_lattice_t a, b, change;
    ibz_t m, f;
    ibz_mat_4x4_t cmp, prod, inpt;
    quat_lattice_init(&a);
    quat_lattice_init(&b);
    quat_lattice_init(&change);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&inpt);
    ibz_mat_4x4_init(&cmp);
    ibz_init(&m);
    ibz_init(&f);

    ibz_set(&(a.basis[0][0]), 1);
    ibz_set(&(a.basis[1][1]), 1);
    ibz_set(&(a.basis[2][2]), 1);
    ibz_set(&(a.basis[3][3]), 1);
    ibz_set(&(a.basis[0][2]), 17);
    ibz_set(&(a.denom), 5);

    ibz_mat_4x4_copy(&b.basis, &a.basis);
    ibz_mat_4x4_scalar_mul(&b.basis, &(a.denom), &b.basis);
    ibz_set(&(b.basis[0][2]), 0);
    ibz_set(&(b.basis[0][1]), 11);
    ibz_set(&(b.basis[1][1]), -7);
    ibz_set(&(b.basis[1][2]), 19);
    ibz_set(&(b.denom), 3);

    ibz_set(&m, 16);
    change_of_basis_from_lattices(&change, &a, &b);
    ibz_invmod(&f, &(a.denom), &m);
    ibz_mat_4x4_scalar_mul(&inpt, &f, &(a.basis));
    ibz_invmod(&f, &(b.denom), &m);
    ibz_mat_4x4_scalar_mul(&cmp, &f, &(b.basis));
    apply_change_of_basis_mod_torsion(&prod, &change, &inpt, &m);
    ibz_mat_4x4_mod(&prod, &prod, &m);
    ibz_mat_4x4_mod(&cmp, &cmp, &m);
    res = res || !ibz_mat_4x4_equal(&cmp, &prod);

    if (res) {
        printf("Torsion signing unit test apply_change_of_basis_mod_torsion failed\n");
    }
    ibz_finalize(&m);
    ibz_finalize(&f);
    quat_lattice_finalize(&a);
    quat_lattice_finalize(&b);
    quat_lattice_finalize(&change);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&inpt);
    ibz_mat_4x4_finalize(&cmp);
    return (res);
}

int
signing_torsion_unittest_common_gen()
{
    int res = 0;
    secret_commitment_t skc;
    secret_commitment_init(&skc);
    secret_key_t skk;
    secret_key_init(&skk);
    ec_curve_t curve;
    ec_basis_t basis;
    fp2_t j, j2;
    ec_curve_init(&curve);

    res = res || !common_part_of_gen(&skc);
    res = res || !test_basis_order_twof(&(skc.canonical_basis_tor), &(skc.curve), TORSION_challenge_torsion);
    res = res || !test_basis_order_twof(&(skc.canonical_basis), &(skc.curve), TORSION_EVEN_POWER);
    res = res || !dim2id2iso_arbitrary_isogeny_evaluation(&basis, &curve, &skc.secret_ideal);
    ec_j_inv(&j, &skc.curve);
    ec_j_inv(&j2, &curve);
    res = res || !fp2_is_equal(&j, &j2);
    matrix_application_even_basis(&basis, &curve, &skc.secret_action, TORSION_EVEN_POWER);
    res = res || !test_basis_order_twof(&basis, &curve, TORSION_EVEN_POWER);

    res = res || !common_part_of_gen(&skk);
    res = res || !test_basis_order_twof(&(skk.canonical_basis_tor), &(skk.curve), TORSION_challenge_torsion);
    res = res || !test_basis_order_twof(&(skk.canonical_basis), &(skk.curve), TORSION_EVEN_POWER);
    res = res || !dim2id2iso_arbitrary_isogeny_evaluation(&basis, &curve, &skk.secret_ideal);
    ec_j_inv(&j, &skk.curve);
    ec_j_inv(&j2, &curve);
    res = res || !fp2_is_equal(&j, &j2);
    matrix_application_even_basis(&basis, &curve, &skk.secret_action, TORSION_EVEN_POWER);
    res = res || !test_basis_order_twof(&basis, &curve, TORSION_EVEN_POWER);

    if (res) {
        printf("Torsion signing unit test common_gen failed\n");
    }
    secret_commitment_finalize(&skc);
    secret_key_finalize(&skk);
    return (res);
}

int
signing_torsion_unittest_commitment_gen()
{
    int res = 0;
    secret_commitment_t sk;
    secret_commitment_init(&sk);
    ec_curve_t curve;
    ec_basis_t basis;
    fp2_t j, js, jp;
    ec_curve_init(&curve);

    res = res || !common_part_of_gen(&sk);
    res = res || !test_basis_order_twof(&(sk.canonical_basis_tor), &(sk.curve), TORSION_challenge_torsion);
    res = res || !test_basis_order_twof(&(sk.canonical_basis), &(sk.curve), TORSION_EVEN_POWER);

    res = res || !dim2id2iso_arbitrary_isogeny_evaluation(&basis, &curve, &sk.secret_ideal);
    ec_j_inv(&js, &sk.curve);
    ec_j_inv(&j, &curve);
    ec_j_inv(&jp, &curve);
    res = res || !fp2_is_equal(&j, &jp);
    res = res || !fp2_is_equal(&j, &js);
    matrix_application_even_basis(&basis, &curve, &sk.secret_action, TORSION_EVEN_POWER);
    res = res || !test_basis_order_twof(&basis, &curve, TORSION_EVEN_POWER);

    if (res) {
        printf("Torsion signing unit test commitment_gen failed\n");
    }
    secret_commitment_finalize(&sk);
    return (res);
}

int
signing_torsion_unittest_action_of_basis()
{
    // here we want to test the function compute_order_action_on_torsion_basis
    // the idea of the function is to test it on orders we know the action of
    // the endomorphisms on the basis
    int res = 1; // TODO: put 0 when there is a test

    // first test with the identity matrix
    ibz_mat_4x4_t order_action_on_torsion_basis;
    ibz_mat_4x4_init(&order_action_on_torsion_basis);
    ibz_mat_2x2_t mat_identity;
    ibz_mat_2x2_init(&mat_identity);
    ibz_mat_2x2_set(&mat_identity, 1, 0, 0, 1);
    /* compute_order_action_on_torsion_basis(&order_action_on_torsion_basis, &mat_identity); */
    /* // now we need to check that we are getting the standard action */
    /* // that is given by  */
    /* // #define ACTION_GEN2 (CURVES_WITH_ENDOMORPHISMS->action_gen2) */
    /* // #define ACTION_GEN3 (CURVES_WITH_ENDOMORPHISMS->action_gen3) */
    /* // #define ACTION_GEN4 (CURVES_WITH_ENDOMORPHISMS->action_gen4) */

    /* // TODO: do comparison modulo the torsion */
    /* // problem : I need to reduce and do it in a new constant */
    /* ibz_t res_mod; */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[0][1], &ACTION_GEN2[0][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[1][1], &ACTION_GEN2[0][1])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[2][1], &ACTION_GEN2[1][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[3][1], &ACTION_GEN2[1][1])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[0][2], &ACTION_GEN3[0][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[1][2], &ACTION_GEN3[0][1])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[2][2], &ACTION_GEN3[1][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[3][2], &ACTION_GEN3[1][1])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[0][3], &ACTION_GEN4[0][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[1][3], &ACTION_GEN4[0][1])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[2][3], &ACTION_GEN4[1][0])); */
    /* res = res || (0 != ibz_cmp(&order_action_on_torsion_basis[3][3], &ACTION_GEN4[1][1])); */

    if (res) {
        printf("Torsion signing unit test action_of_basis failed\n");
    }
    return res;
}

int
signing_torsion_unittest_keygen()
{
    int res = 0;
    secret_key_t sk;
    public_key_t pk;
    secret_key_init(&sk);
    public_key_init(&pk);
    ec_curve_t curve;
    ec_basis_t basis;
    ec_curve_init(&curve);
    fp2_t j, js, jp;

    res = res || !protocols_keygen(&pk, &sk);
    res = res || !test_basis_order_twof(&(sk.canonical_basis_tor), &(sk.curve), TORSION_challenge_torsion);
    res = res || !test_basis_order_twof(&(sk.canonical_basis), &(sk.curve), TORSION_EVEN_POWER);
    res = res || !test_basis_order_twof(&(sk.canonical_basis_tor), &(pk.curve), TORSION_challenge_torsion);
    res = res || !test_basis_order_twof(&(sk.canonical_basis), &(pk.curve), TORSION_EVEN_POWER);

    res = res || !dim2id2iso_arbitrary_isogeny_evaluation(&basis, &curve, &sk.secret_ideal);
    ec_j_inv(&js, &sk.curve);
    ec_j_inv(&j, &curve);
    ec_j_inv(&jp, &curve);
    res = res || !fp2_is_equal(&j, &jp);
    res = res || !fp2_is_equal(&j, &js);
    matrix_application_even_basis(&basis, &curve, &sk.secret_action, TORSION_EVEN_POWER);
    res = res || !test_basis_order_twof(&basis, &curve, TORSION_EVEN_POWER);

    if (res) {
        printf("Torsion signing unit test keygen failed\n");
    }
    secret_key_finalize(&sk);
    return (res);
}

int
signing_torsion_unittest_set_radius()
{
    int res = 0;
    quat_left_ideal_t lideal;
    quat_alg_elem_t elem;
    quat_lattice_t o;
    ibz_t n, r;
    quat_lattice_init(&o);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&elem);
    ibz_init(&n);
    ibz_init(&r);

    quat_alg_elem_set(&elem, 1, 0, 3, 2, -7);
    ibz_set(&n, 4 * 17);
    quat_lideal_create(&lideal, &elem, &n, &MAXORD_O0, &QUATALG_PINFTY);
    quat_lideal_right_order(&o, &lideal, &QUATALG_PINFTY);
    quat_alg_elem_set(&elem, 1, 9, 3, 1, -7);
    quat_lattice_to_standard_basis(&elem, &o, &(elem.coord));
    ibz_set(&n, 9 * 5);
    quat_lideal_create(&lideal, &elem, &n, &o, &QUATALG_PINFTY);

    assert(0 == ibz_cmp(&(lideal.norm), &n));
    set_radius(&n, &n, &QUATALG_PINFTY);
    ibz_div(&n, &r, &n, &QUATALG_PINFTY.p);

    res = res || (!ibz_is_zero(&r));
    res = res || (0 != ibz_cmp(&(lideal.norm), &n));

    if (res) {
        printf("Torsion signing unit test set_radius failed\n");
    }
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&o);
    quat_alg_elem_finalize(&elem);
    ibz_finalize(&n);
    ibz_finalize(&r);
    return (res);
}

static int
equal_keypair(const secret_key_t *sk1, const public_key_t *pk1, const secret_key_t *sk2, const public_key_t *pk2)
{
    int res = 1;
    ec_curve_t c1, c2;
    
    ec_curve_init(&c1);
    ec_curve_init(&c2);
    copy_point(&c1.A24, &pk1->curve.A24);
    fp2_copy(&c1.A, &pk1->curve.A);
    fp2_copy(&c1.C, &pk1->curve.C);
    c1.is_A24_computed_and_normalized = pk1->curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c1);
    copy_point(&c2.A24, &pk2->curve.A24);
    fp2_copy(&c2.A, &pk2->curve.A);
    fp2_copy(&c2.C, &pk2->curve.C);
    c2.is_A24_computed_and_normalized = pk2->curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c2);
    res = res && ec_is_equal(&c1.A24, &c2.A24);

    copy_point(&c1.A24, &sk1->curve.A24);
    fp2_copy(&c1.A, &sk1->curve.A);
    fp2_copy(&c1.C, &sk1->curve.C);
    c1.is_A24_computed_and_normalized = sk1->curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c1);
    copy_point(&c2.A24, &sk2->curve.A24);
    fp2_copy(&c2.A, &sk2->curve.A);
    fp2_copy(&c2.C, &sk2->curve.C);
    c2.is_A24_computed_and_normalized = sk2->curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c2);
    res = res && ec_is_equal(&c1.A24, &c2.A24);

    copy_point(&c1.A24, &sk1->image_curve.A24);
    fp2_copy(&c1.A, &sk1->image_curve.A);
    fp2_copy(&c1.C, &sk1->image_curve.C);
    c1.is_A24_computed_and_normalized = sk1->image_curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c1);
    copy_point(&c2.A24, &sk2->image_curve.A24);
    fp2_copy(&c2.A, &sk2->image_curve.A);
    fp2_copy(&c2.C, &sk2->image_curve.C);
    c2.is_A24_computed_and_normalized = sk2->image_curve.is_A24_computed_and_normalized;
    ec_normalize_curve_and_A24(&c2);
    res = res && ec_is_equal(&c1.A24, &c2.A24);

    res = res && ibz_mat_4x4_equal(&sk1->secret_ideal.lattice.basis, &sk2->secret_ideal.lattice.basis);
    res = res && (0 == ibz_cmp(&sk1->secret_ideal.lattice.denom, &sk2->secret_ideal.lattice.denom));
    res = res && (0 == ibz_cmp(&sk1->secret_ideal.norm, &sk2->secret_ideal.norm));
    res = res && (sk1->secret_ideal.parent_order == sk2->secret_ideal.parent_order);
    res = res && (0 == ibz_cmp(&sk1->secret_action[0][0], &sk2->secret_action[0][0]));
    res = res && (0 == ibz_cmp(&sk1->secret_action[0][1], &sk2->secret_action[0][1]));
    res = res && (0 == ibz_cmp(&sk1->secret_action[1][0], &sk2->secret_action[1][0]));
    res = res && (0 == ibz_cmp(&sk1->secret_action[1][1], &sk2->secret_action[1][1]));
    res = res && fp2_is_equal(&sk1->isom.D, &sk2->isom.D);
    res = res && fp2_is_equal(&sk1->isom.Nx, &sk2->isom.Nx);
    res = res && fp2_is_equal(&sk1->isom.Nz, &sk2->isom.Nz);
    res = res && ec_is_equal(&sk1->canonical_basis.P, &sk2->canonical_basis.P);
    res = res && ec_is_equal(&sk1->canonical_basis.Q, &sk2->canonical_basis.Q);
    res = res && ec_is_equal(&sk1->canonical_basis.PmQ, &sk2->canonical_basis.PmQ);
    res = res && ec_is_equal(&sk1->canonical_basis_tor.P, &sk2->canonical_basis_tor.P);
    res = res && ec_is_equal(&sk1->canonical_basis_tor.Q, &sk2->canonical_basis_tor.Q);
    res = res && ec_is_equal(&sk1->canonical_basis_tor.PmQ, &sk2->canonical_basis_tor.PmQ);
    return (res);
}

int
signing_torsion_unittest_sk_encoding()
{
    int res = 0;
    public_key_t pk, pkt, pkd;
    secret_key_t sk, skt, skd;
    secret_key_init(&sk);
    public_key_init(&pk);
    secret_key_init(&skt);
    public_key_init(&pkt);
    secret_key_init(&skd);
    public_key_init(&pkd);
    protocols_keygen(&pk, &sk);
    copy_basis(&skt.canonical_basis, &sk.canonical_basis);
    copy_basis(&skt.canonical_basis_tor, &sk.canonical_basis_tor);
    copy_curve(&skt.curve, &sk.curve);
    copy_curve(&pkt.curve, &pk.curve);
    pkt.hint_pk = pk.hint_pk;
    copy_curve(&skt.image_curve, &sk.image_curve);
    fp2_copy(&skt.isom.D, &sk.isom.D);
    fp2_copy(&skt.isom.Nx, &sk.isom.Nx);
    fp2_copy(&skt.isom.Nz, &sk.isom.Nz);
    ibz_mat_4x4_copy(&skt.secret_ideal.lattice.basis, &sk.secret_ideal.lattice.basis);
    ibz_copy(&skt.secret_ideal.lattice.denom, &sk.secret_ideal.lattice.denom);
    ibz_copy(&skt.secret_ideal.norm, &sk.secret_ideal.norm);
    skt.secret_ideal.parent_order = sk.secret_ideal.parent_order;
    ibz_mat_2x2_copy(&skt.secret_action, &sk.secret_action);

    byte_t enc[SECRETKEY_BYTES];
    secret_key_to_bytes(enc, &sk, &pk);
    secret_key_from_bytes(&skd, &pkd, enc);
    res = res | !equal_keypair(&sk, &pk, &skd, &pkd);
    res = res | !equal_keypair(&sk, &pk, &skt, &pkt);

    if (res) {
        printf("Torsion signing unit test sk_encoding failed\n");
    }

    secret_key_finalize(&sk);
    // public_key_finalize(&pk);
    secret_key_finalize(&skt);
    // public_key_finalize(&pkt);
    secret_key_finalize(&skd);
    // public_key_finalize(&pkd);
    return (res);
}

int
main()
{
    printf("Running torsion signing unit tests:\n\n");
    int res = 0;
    // keygen
    res = res | signing_torsion_unittest_change_of_basis_from_lattices();
    res = res | signing_torsion_unittest_apply_change_of_basis_mod_torsion();
    res = res | signing_torsion_unittest_common_gen();
    res = res | signing_torsion_unittest_commitment_gen();
    res = res | signing_torsion_unittest_keygen();
    // res = res | signing_torsion_unittest_action_of_basis();
    //  signing
    // res = res | signing_torsion_unittest_set_radius();
    res = res | signing_torsion_unittest_sk_encoding();
    if (res) {
        printf("\nSome tests failed!\n");
    }
    return (res);
}
