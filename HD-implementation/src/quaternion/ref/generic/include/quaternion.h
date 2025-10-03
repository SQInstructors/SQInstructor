/** @file
 *
 * @authors Luca De Feo, Sina Schaeffler
 *
 * @brief Declarations for quaternion algebra operations
 */

#ifndef QUATERNION_H
#define QUATERNION_H

// #include <rng.h>
#include <sqisign_namespace.h>
#include "intbig.h"
#include <assert.h>

/** @defgroup quat_quat Quaternion algebra
 * @{
 */

/** @defgroup quat_vec_t Types for integer vectors and matrices
 * @{
 */

/** @brief Type for vector of 2 big integers
 *
 * @typedef ibz_vec_2_t
 */
typedef ibz_t ibz_vec_2_t[2];

/** @brief Type for vectors of 4 integers
 *
 * @typedef ibz_vec_4_t
 *
 * Represented as a vector of 4 ibz_t (big integer) elements
 */
typedef ibz_t ibz_vec_4_t[4];

/** @brief Type for 2 by 2 matrices of integers
 *
 * @typedef ibz_mat_2x2_t
 *
 * Represented as a matrix of 2 vectors of 2 ibz_t (big integer) elements
 */
typedef ibz_t ibz_mat_2x2_t[2][2];

/** @brief Type for 4 by 4 matrices of integers
 *
 * @typedef ibz_mat_4x4_t
 *
 * Represented as a matrix of 4 vectors of 4 ibz_t (big integer) elements
 */
typedef ibz_t ibz_mat_4x4_t[4][4];
/**
 * @}
 */

/** @defgroup quat_quat_t Types for quaternion algebras
 * @{
 */

/** @brief Type for quaternion algebras
 *
 * @typedef quat_alg_t
 *
 * @struct quat_alg
 *
 * The quaternion algebra ramified at p = 3 mod 4 and ∞.
 */
typedef struct quat_alg
{
    ibz_t p; ///< Prime number, must be = 3 mod 4.
} quat_alg_t;

/** @brief Type for quaternion algebra elements
 *
 * @typedef quat_alg_elem_t
 *
 * @struct quat_alg_elem
 *
 * Represented as a array *coord* of 4 ibz_t integers and a common ibz_t denominator *denom*.
 *
 * The representation is not necessarily normalized, that is, gcd(denom, content(coord)) might not
 * be 1. For getting a normalized representation, use the quat_alg_normalize function
 *
 * The elements are always represented in basis (1,i,j,ij) of the quaternion algebra, with i^2=-1
 * and j^2 = -p
 */
typedef struct quat_alg_elem
{
    ibz_t denom;       ///< Denominator by which all coordinates are divided (big integer, must not be 0)
    ibz_vec_4_t coord; ///< Numerators of the 4 coordinates of the quaternion algebra element in basis (1,i,j,ij)
} quat_alg_elem_t;

/** @brief Type for lattices in dimension 4
 *
 * @typedef quat_lattice_t
 *
 * @struct quat_lattice
 *
 * Represented as a rational (`frac`) times an integreal lattice (`basis`)
 *
 * The basis is such that its columns divided by its denominator are elements of
 * the quaternion algebra, represented in basis (1,i,j,ij) where i^2 = -1, j^2 = -p.
 *
 * All lattices must have full rank (4)
 */
typedef struct quat_lattice
{
    ibz_t denom;         ///< Denominator by which the basis is divided (big integer, must not be 0)
    ibz_mat_4x4_t basis; ///< Integer basis of the lattice  (its columns divided by denom are
                         ///< algebra elements in the usual basis)
} quat_lattice_t;

/** @brief Type for left ideals of maximal orders in quaternion algebras
 *
 * @typedef quat_left_ideal_t
 *
 * @struct quat_left_ideal
 *
 * The basis of the lattice representing it is such that its columns divided by its denominator are
 * elements of the quaternion algebra, represented in basis (1,i,j,ij) where i^2 = -1, j^2 = -p.
 */
typedef struct quat_left_ideal
{
    quat_lattice_t lattice;             ///< lattice representing the ideal
    ibz_t norm;                         ///< norm of the lattice
    const quat_lattice_t *parent_order; ///< should be a maximal order
} quat_left_ideal_t;
/** @}
 */

/** @brief Type for extremal maximal orders
 *
 * @typedef quat_p_extremal_maximal_order_t
 *
 * @struct quat_p_extremal_maximal_order
 *
 * The basis of the order representing it is in hermite normal form, and its columns divid
ed by its denominator are elements of the quaternion algebra, represented in basis (1,z,t,
tz) where z^2 = -q, t^2 = -p.
*/
typedef struct quat_p_extremal_maximal_order
{
    quat_lattice_t order; ///< the order represented as a lattice
    quat_alg_elem_t z;    ///< the element of small discriminant
    quat_alg_elem_t t;    ///< the element of norm p orthogonal to z
    uint32_t q;           ///< the absolute value of the square of z
} quat_p_extremal_maximal_order_t;

/** @brief Type for represent integer parameters
 *
 * @typedef quat_p_extremal_maximal_order_t
 *
 * @struct quat_p_extremal_maximal_order
 *
 */
typedef struct quat_represent_integer_params
{
    int primality_test_iterations;                ///< Primality test iterations
    const quat_p_extremal_maximal_order_t *order; ///< The standard extremal maximal order
    const quat_alg_t *algebra;                    ///< The quaternion algebra
} quat_represent_integer_params_t;

/*************************** Functions *****************************/

/** @defgroup quat_c Constructors and Destructors
 * @{
 */
void quat_alg_init_set(quat_alg_t *alg, const ibz_t *p);
void quat_alg_finalize(quat_alg_t *alg);

void quat_alg_elem_init(quat_alg_elem_t *elem);
void quat_alg_elem_finalize(quat_alg_elem_t *elem);

void ibz_vec_2_init(ibz_vec_2_t *vec);
void ibz_vec_2_finalize(ibz_vec_2_t *vec);

void ibz_vec_4_init(ibz_vec_4_t *vec);
void ibz_vec_4_finalize(ibz_vec_4_t *vec);

void ibz_mat_2x2_init(ibz_mat_2x2_t *mat);
void ibz_mat_2x2_finalize(ibz_mat_2x2_t *mat);

void ibz_mat_4x4_init(ibz_mat_4x4_t *mat);
void ibz_mat_4x4_finalize(ibz_mat_4x4_t *mat);

void quat_lattice_init(quat_lattice_t *lat);
void quat_lattice_finalize(quat_lattice_t *lat);

void quat_left_ideal_init(quat_left_ideal_t *lideal);
void quat_left_ideal_finalize(quat_left_ideal_t *lideal);
/** @}
 */

/** @defgroup quat_printers Print functions for types from the quaternion module
 * @{
 */
void ibz_mat_2x2_print(const ibz_mat_2x2_t *mat);
void ibz_mat_4x4_print(const ibz_mat_4x4_t *mat);
void ibz_vec_2_print(const ibz_vec_2_t *vec);
void ibz_vec_4_print(const ibz_vec_4_t *vec);

void quat_lattice_print(const quat_lattice_t *lat);
void quat_alg_print(const quat_alg_t *alg);
void quat_alg_elem_print(const quat_alg_elem_t *elem);
void quat_left_ideal_print(const quat_left_ideal_t *lideal);

/** @}
 */

/** @defgroup quat_int Integer functions for quaternion algebra
 * @{
 */

/** @defgroup quat_int_mat Integer matrix and vector functions
 * @{
 */

/** @brief Set matrix coefficients to the given integers
 *
 * @param mat Output: Matrix
 * @param a00
 * @param a01
 * @param a10
 * @param a11
 */
void ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a00, int a01, int a10, int a11); // test/dim2

/** @brief Copy matrix
 *
 * @param copy Output: Matrix into which copied will be copied
 * @param copied
 */
void ibz_mat_2x2_copy(ibz_mat_2x2_t *copy, const ibz_mat_2x2_t *copied);

/** @brief in-place 2x2 matrix transposition
 *
 * @param mat Matrix, modified in place
 */
inline void
ibz_mat_2x2_transpose(ibz_mat_2x2_t *mat) {
  ibz_swap(&(*mat)[1][0], &(*mat)[0][1]);
}

/** @brief 2x2 matrix adjugation
 *
 * @param adj adjugate matrix
 * @param mat matrix
 */
void
ibz_mat_2x2_adjugate(ibz_mat_2x2_t *adj, const ibz_mat_2x2_t *mat);

/**
 * @brief determinant of a 2x2 matrix
 *
 * @param det Output integer
 * @param mat Input matrix
 */
void ibz_mat_2x2_det(ibz_t *det, const ibz_mat_2x2_t *mat);

/**
 * @brief a*b for 2x2 integer matrices modulo m
 *
 * @param prod Output matrix
 * @param mat_a Input matrix
 * @param mat_b Input matrix
 * @param m Integer modulo
 */
void ibz_mat_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b, const ibz_t *m);

/**
 * @brief Inverse of 2x2 integer matrices modulo m
 *
 * @param inv Output matrix
 * @param mat Input matrix
 * @param m Integer modulo
 * @return 1 if inverse exists 0 otherwise
 */
int ibz_mat_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);

/** @brief mat*vec in dimension 2 for integers
 *
 * @param res Output vector
 * @param mat Input matrix
 * @param vec Input vector
 */
void ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);

/** @brief Set a vector of 4 integers to given values
 *
 * @param vec Output: is set to given coordinates
 * @param coord0
 * @param coord1
 * @param coord2
 * @param coord3
 */
void ibz_vec_4_set(ibz_vec_4_t *vec, int32_t coord0, int32_t coord1, int32_t coord2, int32_t coord3);

/** @brief set res to values coord0,coord1,coord2,coord3
 *
 * @param res Output: Will contain vector (coord0,coord1,coord2,coord3)
 * @param coord0
 * @param coord1
 * @param coord2
 * @param coord3
 */
void ibz_vec_4_copy_ibz(ibz_vec_4_t *res,
                        const ibz_t *coord0,
                        const ibz_t *coord1,
                        const ibz_t *coord2,
                        const ibz_t *coord3);

/** @brief reduce all coefficients of vec modulo m
 *
 * @param res Output: vec mod m
 * @param vec Input vector
 * @param m Positive integer modulus (m>0)
 */
void ibz_vec_4_mod(ibz_vec_4_t *res, const ibz_vec_4_t *vec, const ibz_t *m);

/** @brief Equality test for 4x4 integer matrices
 *
 * @returns 1 if equal, 0 otherwise
 * @param mat1
 * @param mat2
 */
int ibz_mat_4x4_equal(const ibz_mat_4x4_t *mat1, const ibz_mat_4x4_t *mat2);

/** @brief multiplies all values in vector by same scalar
 *
 * @param prod Output
 * @param scalar
 * @param vec
 */
void ibz_vec_4_scalar_mul(ibz_vec_4_t *prod, const ibz_t *scalar, const ibz_vec_4_t *vec);

/** @brief Copies all values from a 4x4 integer matrix to another one
 *
 * @param new Output: matrix which will have its entries set to mat's entries
 * @param mat Input matrix
 */
void ibz_mat_4x4_copy(ibz_mat_4x4_t *new, const ibz_mat_4x4_t *mat);

/** @brief Reduce all entries of mat modulo m
 *
 * @param new Output: matrix which will have its entries set to mat's entries modulo m
 * @param mat Input matrix
 * @param m Positive integer modulus (m>0)
 */
void ibz_mat_4x4_mod(ibz_mat_4x4_t *res, const ibz_mat_4x4_t *mat, const ibz_t *m);

/** @brief transpose a 4x4 integer matrix
 *
 * @param transposed Output: is set to the transposition of mat
 * @param mat Input matrix
 */
void ibz_mat_4x4_transpose(ibz_mat_4x4_t *transposed, const ibz_mat_4x4_t *mat);

/** @brief a*b for a,b integer 4x4 matrices
 *
 * Naive implementation
 *
 * @param res Output: A 4x4 integer matrix
 * @param a
 * @param b
 */
void ibz_mat_4x4_mul(ibz_mat_4x4_t *res, const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b);

/** @brief Matrix by integer multiplication
 *
 * @param prod Output
 * @param scalar
 * @param mat
 */
void ibz_mat_4x4_scalar_mul(ibz_mat_4x4_t *prod, const ibz_t *scalar, const ibz_mat_4x4_t *mat);

/** @brief divides all values in matrix by same scalar
 *
 * @returns 1 if scalar divided all values in mat, 0 otherwise (division is performed in both cases)
 * @param quot Output
 * @param scalar
 * @param mat
 */
int ibz_mat_4x4_scalar_div(ibz_mat_4x4_t *quot, const ibz_t *scalar, const ibz_mat_4x4_t *mat);

/**
 * @brief mat*vec
 *
 *
 * @param res Output: coordinate vector
 * @param mat Integer 4x4 matrix
 * @param vec Integer vector (coordinate vector)
 *
 * Multiplies 4x4 integer matrix mat by a 4-integers column vector vec
 */
void ibz_mat_4x4_eval(ibz_vec_4_t *res, const ibz_mat_4x4_t *mat, const ibz_vec_4_t *vec);

/**
 * @brief vec*mat
 *
 *
 * @param res Output: coordinate vector.
 * @param vec Integer vector (coordinate vector)
 * @param mat Integer 4x4 matrix
 *
 * Multiplies 4x4 integer matrix mat by a 4-integers row vector vec (on the left)
 */
void ibz_mat_4x4_eval_t(ibz_vec_4_t *res, const ibz_vec_4_t *vec, const ibz_mat_4x4_t *mat);

/** @brief Matrix determinant and a matrix inv such that inv/det is the inverse matrix of the input
 *
 * Implemented following the methof of 2x2 minors explained at Method from
 * https://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf (visited on 3rd of May
 * 2023, 16h15 CEST)
 *
 * @returns 1 if the determinant of mat is not 0 and an inverse was computed, 0 otherwise
 * @param inv Output: Will contain an integer matrix which, dividet by det, will yield the rational
 * inverse of the matrix if it exists, can be NULL
 * @param det Output: Will contain the determinant of the input matrix, can be NULL
 * @param mat Matrix of which the inverse will be computed
 */
int ibz_mat_4x4_inv_with_det_as_denom(ibz_mat_4x4_t *inv, ibz_t *det, const ibz_mat_4x4_t *mat);

/** @brief Solves the linear system given by the matrix system anf the vector target, modulo mod
 *
 *
 * @param res Output: The solution of the linear system
 * @param system Matrix giving the linear system. Must be invertible modulo mod
 * @param target The vector such that system*res = target (modulo mod)
 * @param mod The modulus
 */
void ibz_mat_4x4_solve_system_mod(ibz_vec_4_t *res,
                                  const ibz_mat_4x4_t *system,
                                  const ibz_vec_4_t *target,
                                  const ibz_t *mod);

/** @}
 */

/** @defgroup quat_integer Higher-level integer functions for quaternion algebra
 * @{
 */

/**
 * @brief Generates a random prime
 *
 * A number is accepted as prime if it passes a 30-round Miller-Rabin test.
 * This function is fairly inefficient and mostly meant for tests.
 *
 * @returns 1 if a prime is found, 0 otherwise
 * @param p Output: The prime (if found)
 * @param is3mod4 If 1, the prime is required to be 3 mod 4, if 0 no congruence condition is imposed
 * @param bitsize Maximal size of output prime
 * @param probability_test_iterations Miller-Rabin iteartions for probabilistic primality testing in
 * rejection sampling
 */
int ibz_generate_random_prime(ibz_t *p, int is3mod4, int bitsize, int probability_test_iterations);

/**
 * @brief Find integers x and y such that x^2 + n*y^2 = p
 *
 * Uses Cornacchia's algorithm, should be used  only for prime p
 *
 * @param x Output
 * @param y Output
 * @param n first parameter defining the equation
 * @param p seond parameter defining the equation, must be prime
 * @return 1 if success, 0 otherwise
 */
int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);

/** @}
 */

/** @defgroup quat_qf Quadratic form functions
 * @{
 */

/**
 * @brief Quadratic form evaluation
 *
 * qf and coord must be represented in the same basis.
 *
 * @param res Output: coordinate vector
 * @param qf Quadratic form (4x4 integer matrix)
 * @param coord Integer vector (coordinate vector)
 */
void quat_qf_eval(ibz_t *res, const ibz_mat_4x4_t *qf, const ibz_vec_4_t *coord);
/** @}
 */

/** @}
 */

/** @defgroup quat_quat_f Quaternion algebra functions
 * @{
 */

/** @brief Sets an algebra element to the given integer values, without normalizing it
 *
 * @param elem Output: algebra element of coordinates [coord0,coord1,coord2,coord3] and denominator
 * denom
 * @param denom Denominator, must be non zero
 * @param coord0 Coordinate on 1 (0th vector of standard algebra basis)
 * @param coord1 Coordinate on i (1st vector of standard algebra basis)
 * @param coord2 Coordinate on j (2nd vector of standard algebra basis)
 * @param coord3 Coordinate on ij (3rd vector of standard algebra basis)
 */
void quat_alg_elem_set(quat_alg_elem_t *elem,
                       int32_t denom,
                       int32_t coord0,
                       int32_t coord1,
                       int32_t coord2,
                       int32_t coord3);

/**
 * @brief Copies an algebra element
 *
 * @param copy Output: The element into which another one is copied
 * @param copied Source element copied into copy
 */
void quat_alg_elem_copy(quat_alg_elem_t *copy, const quat_alg_elem_t *copied);

/** @brief Copies the given values into an algebra element, without normalizing it
 *
 * @param elem Output: algebra element of coordinates [coord0,coord1,coord2,coord3] and denominator
 * denom
 * @param denom Denominator, must be non zero
 * @param coord0 Coordinate on 1 (0th vector of standard algebra basis)
 * @param coord1 Coordinate on i (1st vector of standard algebra basis)
 * @param coord2 Coordinate on j (2nd vector of standard algebra basis)
 * @param coord3 Coordinate on ij (3rd vector of standard algebra basis)
 */
void quat_alg_elem_copy_ibz(quat_alg_elem_t *elem,
                            const ibz_t *denom,
                            const ibz_t *coord0,
                            const ibz_t *coord1,
                            const ibz_t *coord2,
                            const ibz_t *coord3);

/** @brief a+b for algebra elements
 *
 * @param res Output
 * @param a Algebra element
 * @param b Algebra element
 */
void quat_alg_add(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);

void quat_alg_mul(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const quat_alg_t *alg);

/** @brief Multiplies algebra element by integer scalar, without normalizing it
 *
 * @param res Output
 * @param scalar Integer
 * @param elem Algebra element
 */
void quat_alg_elem_mul_by_scalar(quat_alg_elem_t *res, const ibz_t *scalar, const quat_alg_elem_t *elem);

/** @brief reduced norm of alg_elem x
 *
 * @param res_num Output: rational which will contain the numerator of the reduced norm of a
 * @param res_denom Output: rational which will contain the denominator of the reduced norm of a (it
 * is 1 if the norm is integer)
 * @param x Algebra element whose norm is computed
 * @param alg The quaternion algebra
 */
void quat_alg_norm(ibz_t *res_num, ibz_t *res_denom, const quat_alg_elem_t *x, const quat_alg_t *alg);

/** @brief Normalize representation of alg_elem x
 *
 * @param x Algebra element whose representation will be normalized
 *
 * Modification of x.
 * Sets coord and denom of x so that gcd(denom, content(coord))=1
 * without changing the value of x = (coord0/denom, coord1/denom, coord2/denom, coord3/denom).
 */
void quat_alg_normalize(quat_alg_elem_t *x);

/**
 * @brief Standard involution in a quaternion algebra
 *
 * @param conj Output: image of x by standard involution of the quaternion algebra alg
 * @param x element of alg whose image is searched
 */
void quat_alg_conj(quat_alg_elem_t *conj, const quat_alg_elem_t *x);

/**
 * @brief Given `x` ∈ `order`, factor it into its primitive and impritive parts
 *
 * Given `x` ∈ `order`, return a coordinate vector `primitive_x` and an integer `content`
 * such that `x` = `content` · Λ `primitive_x`, where Λ is the basis of `order`
 * and `x` / `content` is primitive in `order`.
 *
 * @param primitive_x Output: coordinates of a primitive element of `order` (in `order`'s basis)
 * @param content Output: content of `x`'s coordinate vector in order's basis
 * @param order order of `alg`
 * @param x element of order, must be in `order`
 */
void quat_alg_make_primitive(ibz_vec_4_t *primitive_x,
                             ibz_t *content,
                             const quat_alg_elem_t *x,
                             const quat_lattice_t *order);

// end quat_quat_f
/** @}
 */

/** @defgroup quat_lat_f Lattice functions
 * @{
 */

/**
 * @brief Lattice equality
 *
 * Lattice bases are assumed to be under HNF, but denominators are free.
 *
 * @returns 1 if both lattices are equal, 0 otherwise
 * @param lat1
 * @param lat2
 */
int quat_lattice_equal(const quat_lattice_t *lat1,
                       const quat_lattice_t *lat2); // ideal, lattice, test/lattice, test/ideal

/** @brief Write an algebra element given by its coordinates in a lattice basis in the tandard basis (1,i,j,ij)
 *
 * @param x Output: algebra element given by coord in basis lattice.basis
 * @param coord coordinates of the element in lattice.basis
 * @param lattice lattice whose basis must be the one in which the element is given by coord
 */
void quat_lattice_to_standard_basis(quat_alg_elem_t *x, const quat_lattice_t *lattice, const ibz_vec_4_t *coord);

/**
 * @brief Modifies a lattice to put it in hermite normal form
 *
 * In-place modification of the lattice.
 *
 * @param lat input lattice
 *
 * On a correct lattice this function changes nothing (since it is already in HNF), but it can be
 * used to put a handmade one in correct form in order to use the other lattice functions.
 */
void quat_lattice_hnf(quat_lattice_t *lat);

/** @brief Divides basis and denominator of a lattice by their gcd
 *
 * @param reduced Output
 * @param lat Lattice
 */
void quat_lattice_reduce_denom(quat_lattice_t *reduced, const quat_lattice_t *lat);

void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);

/**
 * @brief Test whether x ∈ lat. If so, compute its coordinates in lat's basis.
 *
 * @param coord Output: Set to the coordinates of x in lat. May be NULL.
 * @param lat The lattice, not necessarily in HNF but full rank
 * @param x An element of the quaternion algebra
 * @return true if x ∈ lat
 */
int quat_lattice_contains(ibz_vec_4_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x);

/**
 * @brief Conjugate of a lattice with basis not in HNF
 *
 * @param conj Output: The lattice conjugate to lat.  ATTENTION: is not under HNF
 * @param lat Input lattice
 */
void quat_lattice_conjugate_without_hnf(quat_lattice_t *conj, const quat_lattice_t *lat);

/**
 * @brief Multiply a lattice and an algebra element
 *
 * The element is multiplied to the right of the lattice
 *
 * @param prod Output: Lattice lat*elem
 * @param lat Input lattice
 * @param elem Algebra element
 * @param alg The quaternion algebra
 */
void quat_lattice_alg_elem_mul(quat_lattice_t *prod,
                               const quat_lattice_t *lat,
                               const quat_alg_elem_t *elem,
                               const quat_alg_t *alg);

/** @brief a*b for lattices
 *
 * @param res Output
 * @param lat1 Lattice
 * @param lat2 Lattice
 * @param alg The quaternion algebra
 */
void quat_lattice_mul(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg);

/** @brief The index of sublat into overlat
 *
 * Assumes inputs are in HNF.
 *
 * @param index Output
 * @param sublat A lattice in HNF, must be sublattice of overlat
 * @param overlat A lattice in HNF, must be overlattice of sublat
 */
void quat_lattice_index(ibz_t *index, const quat_lattice_t *sublat, const quat_lattice_t *overlat);

/**
 * @brief Sample from the intersection of a lattice with a ball
 *
 * Sample a uniform non-zero vector of norm ≤ `radius` from the lattice.
 *
 * @param res Output: sampled quaternion from the lattice
 * @param lattice Input lattice
 * @param alg The quaternion algebra
 * @param radius The ball radius (quaternion norm)
 * @return 0 if an error occurred (ball too small or RNG error), 1 otherwise
 */
int quat_lattice_sample_from_ball(quat_alg_elem_t *res,
                                  const quat_lattice_t *lattice,
                                  const quat_alg_t *alg,
                                  const ibz_t *radius);

// end quat_lat_f
/** @}
 */

/** @defgroup quat_lideal_f Functions for left ideals
 * @{
 */

/** @defgroup quat_lideal_c Creating left ideals
 * @{
 */

/**
 * @brief Left ideal of order, generated by x and N as order*x+order*N
 *
 * @param lideal Output: left ideal
 * @param alg quaternion algebra
 * @param order maximal order of alg whose left ideal is searched
 * @param x generating element. Must be non-zero
 * @param N generating integer
 *
 * Creates the left ideal in order generated by the element x and the integer N.
 * If x is not divisible (inside the order) by any integer divisor n>1 of N,
 * then the norm of the output ideal is N.
 *
 */
void quat_lideal_create(quat_left_ideal_t *lideal,
                        const quat_alg_elem_t *x,
                        const ibz_t *N,
                        const quat_lattice_t *order,
                        const quat_alg_t *alg);

/** @}
 */

/** @defgroup quat_lideal_gen Generators of left ideals
 * @{
 */

/**
 * @brief Generator of 'lideal'
 *
 * @returns 1 if such a generator was found, 0 otherwise
 * @param gen Output: non scalar generator of lideal
 * @param lideal left ideal
 * @param alg the quaternion algebra
 *
 * Ideal is generated by gen and the ideal's norm
 *
 * Bound has as default value QUATERNION_lideal_generator_search_bound
 */
int quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *alg);
/** @}
 */

/** @defgroup quat_lideal_op Operations on left ideals
 * @{
 */

/**
 * @brief Copies an ideal
 *
 * @param copy Output: The ideal into which another one is copied
 * @param copied Source ideal copied into copy. The parent order is not copied (only the pointer).
 */
void quat_lideal_copy(quat_left_ideal_t *copy, const quat_left_ideal_t *copied);

/**
 * @brief  Conjugate of a left ideal (not in HNF)
 *
 * @param conj Output: Ideal conjugate to lideal, with norm and parent order correctly set, but its
 * lattice not in HNF
 * @param new_parent_order Output: Will be set to the right order of lideal, and serve as parent
 * order for conj (so must have at least the lifetime of conj)
 * @param lideal input left ideal (of which conj will be the conjugate)
 * @param alg the quaternion algebra
 */
void quat_lideal_conjugate_without_hnf(quat_left_ideal_t *conj,
                                       quat_lattice_t *new_parent_order,
                                       const quat_left_ideal_t *lideal,
                                       const quat_alg_t *alg);

/**
 * @brief  Intersection of two left ideals
 *
 * @param intersection Output: Left ideal which is the intersection of the 2 inputs
 * @param lideal1 left ideal
 * @param lideal2 left ideal
 * @param alg the quaternion algebra
 */
void quat_lideal_inter(quat_left_ideal_t *intersection,
                       const quat_left_ideal_t *lideal1,
                       const quat_left_ideal_t *lideal2,
                       const quat_alg_t *alg);

/**
 * @brief  Sum of two left ideals
 *
 * @param sum Output: Left ideal which is the sum of the 2 inputs
 * @param lideal1 left ideal
 * @param lideal2 left ideal
 * @param alg the quaternion algebra
 */
void quat_lideal_add(quat_left_ideal_t *sum,
                     const quat_left_ideal_t *lideal1,
                     const quat_left_ideal_t *lideal2,
                     const quat_alg_t *alg);

/**
 * @brief  Left ideal product of left ideal I and element alpha
 *
 * @param product Output: lideal I*alpha, must have integer norm
 * @param lideal left ideal
 * @param alpha element multiplied to lideal to get the product ideal
 * @param alg the quaternion algebra
 *
 * I*alpha where I is a left-ideal and alpha an element of the algebra
 *
 * The resulting ideal must have an integer norm
 *
 */
void quat_lideal_mul(quat_left_ideal_t *product,
                     const quat_left_ideal_t *lideal,
                     const quat_alg_elem_t *alpha,
                     const quat_alg_t *alg);

/** @brief Computes the inverse ideal (for a left ideal of a maximal order) without putting it under
 * HNF
 *
 * This function returns a lattice not under HNF. For careful internal use only
 *
 * Computes the inverse ideal for lideal as conjugate(lideal)/norm(lideal)
 *
 * @param inv Output: lattice which is lattice representation of the inverse ideal of lideal
 * ATTENTION: is not under HNF. hnf computation must be applied before using lattice functions on it
 * @param lideal Left ideal of a maximal order in alg
 * @param alg The quaternion algebra
 */
void quat_lideal_inverse_lattice_without_hnf(quat_lattice_t *inv,
                                             const quat_left_ideal_t *lideal,
                                             const quat_alg_t *alg);

/**
 * @brief  Left ideal product of left ideal I1 and left ideal I2
 *
 * Assumes the inputs are compatible: the right order of I1 matches the left order of I2
 *
 * @param product Output: lideal I1*I2
 * @param lideal1 left ideal
 * @param lideal2 left ideal
 * @param alg the quaternion algebra
 */
void quat_lideal_lideal_mul(quat_left_ideal_t *product,
                            const quat_left_ideal_t *lideal1,
                            const quat_left_ideal_t *lideal2,
                            const quat_alg_t *alg);

/**
 * @brief L2-reduce the basis of the left ideal, without considering its denominator
 *
 * This function reduce the basis of the lattice of the ideal, but it does completely ignore its
 * denominator. So the outputs of this function must still e divided by the appropriate power of
 * lideal.lattice.denom.
 *
 * Implements the L2 Algorithm of Nguyen-Stehlé, also known as fplll:
 * https://iacr.org/archive/eurocrypt2005/34940217/34940217.pdf
 *
 * Parameters are in lll/lll_internals.h
 *
 * @param reduced Output: Lattice defining the ideal, which has its basis in a lll-reduced form.
 * Must be divided by lideal.lattice.denom before usage
 * @param gram Output: Matrix of the quadratic form given by the norm on the basis of the reduced
 * ideal, divided by the norm of the ideal
 * @param lideal ideal whose basis will be reduced
 * @param alg the quaternion algebra
 */
void quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced,
                              ibz_mat_4x4_t *gram,
                              const quat_left_ideal_t *lideal,
                              const quat_alg_t *alg); // replaces lideal_lll

/**
 * @brief Multplies two ideals and L2-reduces the lattice of the result
 *
 * Implements the L2 Algorithm of Nguyen-Stehlé, also known as fplll:
 * https://iacr.org/archive/eurocrypt2005/34940217/34940217.pdf
 *
 * Parameters are in lll/lll_internals.h
 *
 * @param prod Output: The product ideal with its lattice basis being L2-reduced
 * @param gram Output: Gram matrix of the reduced norm (as quadratic but not bilinear form) on the
 * basis of prod, divided by the norm of prod
 * @param lideal1 Ideal at left in the product
 * @param lideal2 Ideal at right in the product
 * @param alg The quaternion algebra
 */
void quat_lideal_lideal_mul_reduced(quat_left_ideal_t *prod,
                                    ibz_mat_4x4_t *gram,
                                    const quat_left_ideal_t *lideal1,
                                    const quat_left_ideal_t *lideal2,
                                    const quat_alg_t *alg);

/**
 * @brief  Right order of a left ideal
 *
 * @param order Output: right order of the given ideal
 * @param lideal left ideal
 * @param alg the quaternion algebra
 */
void quat_lideal_right_order(quat_lattice_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg);

/** @brief Test if order is maximal
 *
 * Checks if the discriminant of the order equals the prime p defining the quaternion algebra.
 *
 * It is not verified whether the order is really an order. The output 1 only means that if it is an
 * order, then it is maximal.
 *
 * @returns 1 if order is maximal (assuming it is an order), 0 otherwise
 * @param order An order of the quaternion algebra (assumes to be an order, this is not tested)
 * @param alg The quaternion algebra
 */
int quat_order_is_maximal(const quat_lattice_t *order,
                          const quat_alg_t *alg); // ideal (only in asserts)

/**
 * @brief res = orig*(conjugate(factor)/norm(orig)) equivalent ideal to orig
 *
 * @returns 1 if the computation succeeded and 0 otherwise
 * @param res Output: An ideal equivalent to orig by factor
 * @param alg The quaternion algebra
 * @param orig Input ideal
 * @param factor Quternion, must be in orig
 */
void quat_lideal_equivalent(quat_left_ideal_t *res,
                            const quat_left_ideal_t *orig,
                            const quat_alg_elem_t *factor,
                            const quat_alg_t *alg);

/**
 * @brief Replaces an ideal by a smaller equivalent one of prime norm
 *
 * @returns 1 if the computation succeeded and 0 otherwise
 * @param lideal In- and Output: Ideal to be replaced
 * @param alg The quaternion algebra
 * @param primality_num_iter number of repetition for primality testing
 * @param equiv_bound_coeff bound on the coefficients for the candidates
 */
int quat_lideal_prime_norm_reduced_equivalent(quat_left_ideal_t *lideal,
                                              const quat_alg_t *alg,
                                              const int primality_num_iter,
                                              const int equiv_bound_coeff);

/** @}
 */

// end quat_lideal_f
/** @}
 */

/** @defgroup quat_normeq Functions specific to special extremal maximal orders
 * @{
 */

/**
 * @brief Representing an integer by the quadratic norm form of a maximal extremal order
 *
 * @returns 1 if the computation succeeded
 * @param gamma Output: a quaternion element
 * @param n_gamma Target norm of gamma. n_gamma must be odd. If n_gamma/(p*params.order->q) <
 * 2^QUAT_repres_bound_input failure is likely
 * @param non_diag If set to 1 (instead of 0) and the order is O0, an additional property is ensured
 * @param params Represent integer parameters specifying the algebra, the special extremal order,
 * the number of trials for finding gamma and the number of iterations of the primality test.
 * Special requirements apply if non-diag is set to 1
 *
 * This algorithm finds a primitive quaternion element gamma of n_gamma inside any maximal extremal
 * order. Failure is possible. Most efficient for the standard order.
 *
 *  If non-diag is set to 1,this algorithm finds a primitive quaternion element gamma with some
 * special properties used in fixed degree isogeny of n_gamma inside any maximal extremal order such
 * that params->order->q=1 mod 4. Failure is possible. Most efficient for the standard order. The
 * most important property is to avoid diagonal isogenies, meaning that the gamma returned by the
 * algorithm must not be contained inside ZZ + 2 O where O is the maximal order params->order When O
 * is the special order O0 corresponding to j=1728, we further need to avoid endomorphisms of E0xE0
 * and there is another requirement
 *
 * If non-diag is set to 1, the number of trials for finding gamma (in params), the number of
 * iterations of the primality test and the value of params->order->q is required to be 1 mod 4
 */
int quat_represent_integer(quat_alg_elem_t *gamma,
                           const ibz_t *n_gamma,
                           int non_diag,
                           const quat_represent_integer_params_t *params);

/** @brief Basis change to (1,i,(i+j)/2,(1+ij)/2) for elements of O0
 *
 * Change the basis in which an element is give from 1,i,j,ij to (1,i,(i+j)/2,(1+ij)/2) the ususal
 * basis of the special maximal order O0 Only for elements of O0
 *
 * @param vec Output: Coordinates of el in basis (1,i,(i+j)/2,(1+ij)/2)
 * @param el Imput: An algebra element in O0
 */
void quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el);

/**
 * @brief Random O0-ideal of given norm
 *
 * Much faster if norm is prime and is_prime is set to 1
 *
 * @param lideal Output: O0-ideal of norm norm
 * @param norm Norm of the ideal to be found
 * @param is_prime Indicates if norm is prime: 1 if it is, 0 otherwise
 * @param params Represent Integer parameters from the level-dependent constants
 * @param prime_cofactor Prime distinct from the prime p defining the algebra but of similar size
 * and coprime to norm. If is_prime is 1, it might be NULL.
 * @returns 1 if success, 0 if no ideal found or randomness failed
 */
int quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,
                                             const ibz_t *norm,
                                             int is_prime,
                                             const quat_represent_integer_params_t *params,
                                             const ibz_t *prime_cofactor);
// end quat_normeq
/** @}
 */
// end quat_quat
/** @}
 */

#endif
