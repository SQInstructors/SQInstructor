# add_executable(sqisign_test_verification_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_verify.c)
# target_link_libraries(sqisign_test_verification_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_MP} ${LIB_VERIFICATION_${SVARIANT_UPPER}} ${LIB_QUATERNION} ${GMP} sqisign_common_test sqisign_common_sys)
# target_include_directories(sqisign_test_verification_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_MP} ${INC_QUATERNION} ${INC_SIGNATURE_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE} ${PROJECT_SOURCE_DIR}/src/mini-gmp)

# add_test(sqisign_test_verification_${SVARIANT_LOWER} sqisign_test_verification_${SVARIANT_LOWER} 3)
