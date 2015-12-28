
set(EXAMPLE_FILES
        README.userweight
        ${TOP_SRC_DIR}/src/exposure_weight_base.h
        fixed_power.cc linear.cc
        root_finding.hh sinc.cc
        tabular.cc templated_power.cc
        variable_power.cc
)

if(ENABLE_OPENCL)
  list(APPEND EXAMPLE_FILES
    generalized_gaussian.cl
    variable_power.cl
  )
ENDIF()

install(FILES ${EXAMPLE_FILES}
    DESTINATION ${DocumentsDirectory}/examples/enfuse
)

# TODO: install Makefile.userweight
