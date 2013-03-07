# Search for OpenCL header and librariy
#
# It will set
#  OPENCL_FOUND       - system has OpenCL
#  OPENCL_INCLUDE_DIR - the OpenCL include directory
#  OPENCL_LIBRARY     - link these to use OpenCL
#

IF(WIN32)
  # search for AMD/nVidia/Intel SDK with OpenCL headers and lib
  SET(CL_INCLUDE_PREFIX "CL")
  FIND_PATH(OPENCL_INCLUDE_DIR
    NAMES
        CL/cl.hpp 
    PATHS
        $ENV{AMDAPPSDKROOT}/include
        $ENV{INTELOCLSDKROOT}/include
        $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/inc
  )

  IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # 32 bit libs
    SET(OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        $ENV{AMDAPPSDKROOT}/lib/x86
        $ENV{INTELOCLSDKROOT}/lib/x86
        $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/lib/Win32
    )
  ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
    # 64 bit libs
    SET(OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        $ENV{AMDAPPSDKROOT}/lib/x86_64
        $ENV{INTELOCLSDKROOT}/lib/x64
        $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/lib/x64
    )
  ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)

  FIND_LIBRARY(
    OPENCL_LIBRARY
    NAMES OpenCL
    PATHS ${OPENCL_LIB_SEARCH_PATH}
  )
ELSE(WIN32)

  IF(APPLE)
    SET(CL_INCLUDE_PREFIX "OpenCL")
  ELSE()
    SET(CL_INCLUDE_PREFIX "CL")
  ENDIF()
  FIND_PATH(OPENCL_INCLUDE_DIR
    NAMES
      ${CL_INCLUDE_PREFIX}/cl.h
    PATHS
      /usr/local/include
      /usr/include
      "/usr/local/cuda/include"
  )
  FIND_LIBRARY(OPENCL_LIBRARY
    NAMES
        OpenCL CL clparser
    PATHS
        /usr/local/lib
        /usr/lib
  )
ENDIF(WIN32)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenCL DEFAULT_MSG OPENCL_LIBRARY OPENCL_INCLUDE_DIR)

IF(OPENCL_FOUND)
  STRING(TOUPPER ${CL_INCLUDE_PREFIX} _INC_PREFIX)
  SET(HAVE_${_INC_PREFIX}_CL_HPP 1)
ENDIF()

MARK_AS_ADVANCED(OPENCL_INCLUDE_DIR OPENCL_LIBRARY)

