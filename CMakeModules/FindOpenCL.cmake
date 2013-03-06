# Search for OpenCL header and librariy
#
# It will set
#  OPENCL_FOUND       - system has OpenCL
#  OPENCL_INCLUDE_DIR - the OpenCL include directory
#  OPENCL_LIBRARY     - link these to use OpenCL
#

IF(WIN32)
# search for AMD/nVidia/Intel SDK with OpenCL headers and lib
  FIND_PATH(OPENCL_INCLUDE_DIR
    NAMES
        CL/cl.h OpenCL/cl.h
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

  foreach(_inc "CL" "OpenCL")
    FIND_PATH(OPENCL_INCLUDE_DIR 
      NAMES
        ${_inc}/cl.h
      PATHS
        /usr/local/include
        /usr/include
        "/usr/local/cuda/include"
    )
    if(OPENCL_INCLUDE_DIR)
      string(TOUPPER ${_inc} _I)
      set(HAVE_${_I}_OPENCL_H 1)
    endif()
  endforeach()
  
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

MARK_AS_ADVANCED(OPENCL_INCLUDE_DIR OPENCL_LIBRARY)

