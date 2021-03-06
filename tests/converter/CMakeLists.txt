
# Add tests for conversion of GAMESS to QMCPACK input (via convert4qmc)

FUNCTION(ADD_CONVERTER_TEST test_name)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
    SET(EXE_NAME "${qmcpack_BINARY_DIR}/bin/convert4qmc")
    IF ( USE_MPI )
        SET(EXE_NAME "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${qmcpack_BINARY_DIR}/bin/convert4qmc")
    ENDIF()
    ADD_TEST(NAME converter_${test_name} COMMAND python converter_test.py "${CMAKE_CURRENT_BINARY_DIR}/${test_name}" --exe "${EXE_NAME}")
    SET_TESTS_PROPERTIES(converter_${test_name} PROPERTIES TIMEOUT 120)
    
ENDFUNCTION()

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/converter_test.py" "${CMAKE_CURRENT_BINARY_DIR}")
ADD_CONVERTER_TEST(test_He_sto3g)
ADD_CONVERTER_TEST(test_Be_sto3g)
ADD_CONVERTER_TEST(test_Be_ccd)
ADD_CONVERTER_TEST(test_O_ext)
ADD_CONVERTER_TEST(test_C_sto3g)
ADD_CONVERTER_TEST(test_HCNp)
ADD_CONVERTER_TEST(test_aldet1)
ADD_CONVERTER_TEST(test_aldet5)
