/** @file
 *
 * @brief The verification protocol
 */

#ifndef VERIFICATION_H
#define VERIFICATION_H

#include <sqisign_namespace.h>
#include <ec.h>

/** @defgroup verification SQIsignHD verification protocol
 * @{
 */

/** @defgroup verification_t Types for SQIsignHD verification protocol
 * @{
 */
typedef unsigned char byte_t;
typedef digit_t scalar_t[NWORDS_ORDER];
typedef scalar_t scalar_mtx_2x2_t[2][2];

/** @brief Type for the signature
 *
 * @typedef signature_t
 *
 * @struct signature
 *
 */
typedef struct signature
{
    fp2_t E_aux_A; // the Montgomery A-coefficient for the auxiliary curve
    uint8_t hint_aux;
    scalar_mtx_2x2_t mat_B_aux_can_to_B_aux; /// the matrix of the desired basis on the aux curve
    uint8_t hint_com;                        // hint for E_com
    scalar_t scalar;                         // correction scalar response isogeny
    uint16_t resp_length;                    /// length of the response as 22 isogenies
} signature_t;

/** @brief Type for the public keys
 *
 * @typedef public_key_t
 *
 * @struct public_key
 *
 */
typedef struct public_key
{
    ec_curve_t curve; /// the normalized A-coefficient of the Montgomery curve
    uint8_t hint_pk;
    // memo: the torsion basis will be given by multiplying the canonical basis by 2**(e - e_tor)
} public_key_t;

/** @}
 */

/*************************** Functions *****************************/

void public_key_init(public_key_t *pk);

void hash_to_challenge(scalar_mtx_2x2_t *chall,
                       const public_key_t *pk,
                       const unsigned char *message,
                       int message_length,
                       const ec_curve_t *commitment,
                       int torsion_exponent);

byte_t *ec_curve_to_bytes(byte_t *enc, const ec_curve_t *curve);

const byte_t *ec_curve_from_bytes(ec_curve_t *curve, const byte_t *enc);

byte_t *ec_basis_to_bytes(byte_t *enc, const ec_basis_t *basis);

const byte_t *ec_basis_from_bytes(ec_basis_t *basis, const byte_t *enc);

void scalar_mtx_2x2_mod_2exp(scalar_mtx_2x2_t *x, int exp, int nwords);

// return 1 if invertible, 0 else
// cannot handle negative numbers in inputs
int scalar_mtx_2x2_inv_mod_2exp_up_to_scalar(scalar_mtx_2x2_t *inv, const scalar_mtx_2x2_t *x, int exp, int nwords);

// cannot handle negative numbers in inputs
void scalar_mtx_2x2_mul_mod_2exp(scalar_mtx_2x2_t *prod,
                                 const scalar_mtx_2x2_t *a,
                                 const scalar_mtx_2x2_t *b,
                                 int exp,
                                 int nwords);

/**
 * @brief Verification
 *
 * @param sig signature
 * @param pk public key
 * @param m message
 * @param l size
 * @returns 1 if the signature verifies, 0 otherwise
 */
int protocols_verify(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l);

/*************************** Encoding *****************************/

/** @defgroup encoding Encoding and decoding functions
 * @{
 */

/**
 * @brief Encodes a signature as a byte array
 *
 * @param enc : Byte array to encode the signature in
 * @param sig : Signature to encode
 */
void signature_to_bytes(unsigned char *enc, const signature_t *sig);

/**
 * @brief Decodes a signature from a byte array
 *
 * @param sig : Structure to decode the signature in
 * @param enc : Byte array to decode
 */
void signature_from_bytes(signature_t *sig, const unsigned char *enc);

/**
 * @brief Encodes a public key as a byte array
 *
 * @param enc : Byte array to encode the public key in
 * @param pk : Public key to encode
 */
unsigned char *public_key_to_bytes(unsigned char *enc, const public_key_t *pk);

/**
 * @brief Decodes a public key from a byte array
 *
 * @param pk : Structure to decode the public key in
 * @param enc : Byte array to decode
 */
const unsigned char *public_key_from_bytes(public_key_t *pk, const unsigned char *enc);

/** @}
 */

/** @}
 */

#endif
