/** @file
 *
 * @brief The key generation and signature protocols
 */

#ifndef SIGNATURE_H
#define SIGNATURE_H

#include <sqisign_namespace.h>
#include <ec.h>
#include <quaternion.h>
#include <verification.h>

/** @defgroup signature SQIsignHD key generation and signature protocols
 * @{
 */
/** @defgroup signature_t Types for SQIsignHD key generation and signature protocols
 * @{
 */

/** @brief Secret data associated to a public curve (commitment or public key)
 *
 * @typedef secret_common_t
 *
 * @struct secret_common
 *
 */
typedef struct secret_common
{
    ec_curve_t image_curve;             // the image curve of φ
    ec_isom_t isom;                     // the isomorphism to the standardized curve
    ec_curve_t curve;                   // the standardized public curve
    quat_left_ideal_t secret_ideal;     // the ideal I_φ, corresponding to an isogeny φ
    ec_basis_t canonical_basis;         // the canonical basis of the public curve (order 2^e)
    ec_basis_t canonical_basis_tor;     // the canonical basis of the torsion of the public curve
                                        // (order 2^etor),should always be just 2**(e - e_tor) * canonical_basis
    ibz_mat_2x2_t secret_action;        // the action of the secret ideal/isogeny on the
                                        // canonical bases, defined as
                                        // secret_action · canonical_basis = φ(B₀)
                                        // where B₀ is the canonical basis of E₀
} secret_common_t;

typedef secret_common_t secret_commitment_t;
typedef secret_common_t secret_key_t;
/** @}
 */

/*************************** Functions *****************************/

void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);

void secret_commitment_init(secret_commitment_t *sk);
void secret_commitment_finalize(secret_commitment_t *sk);


/** @brief Convert to scalar_mtx_2x2_t type
 *
 * Only works for small enough and positie coefficients
 *
 * @param a Output: scalar_mtx_2x2_t type matrix containing the same coefficients as x
 * @param x Must have small positie coefficients
 */
void ibz_mat_2x2_to_digit_array(scalar_mtx_2x2_t *a, const ibz_mat_2x2_t *x);

/** @brief Convert from scalar_mtx_2x2_t type
 *
 * Only works for small enough and positive coefficients
 *
 * @param a Output: matrix containing the same coefficients as x
 * @param x
 */
void ibz_mat_2x2_from_digit_array(ibz_mat_2x2_t *a, const scalar_mtx_2x2_t *x);

/**
 * @brief Key generation
 *
 * @param pk Output: will contain the public key
 * @param sk Output: will contain the secret key
 * @returns 1 if success, 0 otherwise
 */
int protocols_keygen(public_key_t *pk, secret_key_t *sk);

/**
 * @brief Key generation for the commitment
 *
 * @param cmt Output: will contain the public key
 * @param sk_cmt Output: will contain the secret key
 * @returns 1 if success, 0 otherwise
 */
int protocols_commitment_gen(signature_t *sig, secret_commitment_t *sk_cmt);

/**
 * @brief Signature computation
 *
 * @param sig Output: will contain the signature
 * @param sk secret key
 * @param pk public key
 * @param m message
 * @param l size
 * @returns 1 if success, 0 otherwise
 */
int protocols_sign(signature_t *sig, const public_key_t *pk, secret_key_t *sk, const unsigned char *m, size_t l);

/*************************** Encoding *****************************/

/** @defgroup encoding Encoding and decoding functions
 * @{
 */

/**
 * @brief Encodes a secret key as a byte array
 *
 * @param enc : Byte array to encode the secret key (including public key) in
 * @param sk : Secret key to encode
 * @param pk : Public key to encode
 */
void secret_key_to_bytes(unsigned char *enc, const secret_key_t *sk, const public_key_t *pk);

/**
 * @brief Decodes a secret key (and public key) from a byte array
 *
 * @param sk : Structure to decode the secret key in
 * @param pk : Structure to decode the public key in
 * @param enc : Byte array to decode
 */
void secret_key_from_bytes(secret_key_t *sk, public_key_t *pk, const unsigned char *enc);

/** @}
 */

/** @defgroup signature_torsion Torsion SQIsign headers and constants
 * @{
 */


// determine the actual radius of the ball to sample from
// for scalar it is sqrt(pN³)
void set_radius(ibz_t *radius, const ibz_t *Ntors, const quat_alg_t *alg);

/** @defgroup compute_solutions_subfunctions Subfunctions of compute_solutions, to be moved to a bettr place
 * @{
 */

void lattice_solutions_for_scalar(quat_lattice_t *lattice_solutions,
                                  const quat_left_ideal_t *I_solutions,
                                  const ibz_t *torsion);

/** @}
 */

/** @defgroup keygen_subfunctions Subfunctions of keygen (and sometime signature)
 * @{
 */
// assumes torsion to be a power of 2
void change_of_basis_from_lattices(quat_lattice_t *res,
                                   const quat_lattice_t *I_solutions_lattice,
                                   const quat_lattice_t *sk_order);

void apply_change_of_basis_mod_torsion(ibz_mat_4x4_t *res,
                                       const quat_lattice_t *change_of_basis,
                                       const ibz_mat_4x4_t *basis,
                                       const ibz_t *torsion);

uint16_t common_part_of_gen(secret_common_t *sk);

/** @}
 */

/** @}
 */
/** @}
 */

#endif
