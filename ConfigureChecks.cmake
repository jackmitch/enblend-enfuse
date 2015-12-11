# This file is part of enblend.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009-2015, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
#

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(CheckSymbolExists)
include(CheckFunctionExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)
include(TestBigEndian)

MACRO(AddHeaderTestToConfig header var)
STRING(CONFIGURE "/* Define if you have the <${header}> header file. */\n#cmakedefine ${var} 1\n" _str)
SET(CMAKE_HEADER_EXISTS "${CMAKE_HEADER_EXISTS}${_str}\n")
UNSET(_str)
ENDMACRO()

MACRO(AddFunctionTestToConfig function var)
STRING(CONFIGURE "/* Define if you have the '${function}' function. */\n#cmakedefine ${var} 1\n" _str)
SET(CMAKE_FUNCTION_EXISTS "${CMAKE_FUNCTION_EXISTS}${_str}\n")
UNSET(_str)
ENDMACRO()


#Check for some include files and set appropriate variables
# e.g. "sys/dir.h" found => "HAVE_SYS_DIR_H" set to 1
foreach(_fl "dirent.h" "fenv.h"
    "inttypes.h" "limits.h"
    "memory.h" "stdint.h" "stdlib.h" "stdbool.h" "strings.h" "string.h"
    "sys/stat.h" "sys/types.h" "unistd.h" "windows.h" "sys/times.h")
  string(REGEX REPLACE "[/\\.]" "_" _var ${_fl})
  string(TOUPPER "${_var}" _FLN)
  check_include_file_cxx("${_fl}" "HAVE_${_FLN}" )
  AddHeaderTestToConfig("${_fl}" "HAVE_${_FLN}")
endforeach()

#Check for functions
set(CMAKE_REQUIRED_LIBRARIES -lm)
foreach(_fc atexit fesetround floor fseeko lrint lrintf
    memset mkstemp pow select sqrt malloc
    strchr strcspn strdup strerror strerror_r strpbrk strrchr strtol strtok_r strerror_r strtoul)
  string(TOUPPER "${_fc}" _FC)
  check_function_exists(${_fc} "HAVE_${_FC}")
  AddFunctionTestToConfig("${_fc}" "HAVE_${_FC}")
endforeach()

if(HAVE_DIRENT_H)
  check_cxx_source_compiles(
    "
    #include <dirent.h>
    DIR *DI;
    int i = closedir(DI);
    int main(){return(0);}
    "
    CLOSEDIR_INT)
endif(HAVE_DIRENT_H)
if (NOT CLOSEDIR_INT)
  set(CLOSEDIR_VOID 1)
endif()

CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
AddFunctionTestToConfig("clock_gettime" "HAVE_CLOCK_GETTIME")

# Check for restrict keyword
# Builds the macro A_C_RESTRICT form automake
set(RESTRICT)
foreach(ac_kw __restrict __restrict__ _Restrict restrict)
  check_cxx_source_compiles(
  "
  int foo (int * ${ac_kw} ip) {
    return ip[0];
  }
  int main(){
    int s[1];
    int * ${ac_kw} t = s;
    t[0] = 0;
    return foo(t);
  }
  "
  RESTRICT_${ac_kw}_FOUND)
  if(RESTRICT_${ac_kw}_FOUND)
    set(RESTRICT ${ac_kw})
    break()
  endif()
endforeach()

if(VIGRA_FOUND AND NOT VIGRA_VERSION_CHECK)
  unset(VIGRA_SETIMAGEINDEX CACHE)
  set(CMAKE_REQUIRED_INCLUDES ${VIGRA_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${VIGRA_LIBRARIES})
  check_cxx_source_compiles(
  "
  #include <vigra/imageinfo.hxx>
  #include <vigra/impexalpha.hxx>

  vigra::ImageImportInfo info(\"image.tif\");
  int main(){
    vigra::BRGBImage image;
    vigra::BImage alpha;
    vigra::impexListFormats();
    info.setImageIndex(99);
    vigra::importImageAlpha(info, destImage(image), destImage(alpha));
    return(0);
  }
  "
  VIGRA_SETIMAGEINDEX)
  if(NOT VIGRA_SETIMAGEINDEX)
    message(FATAL_ERROR "Vigra library too old")
  endif()
endif()

message(STATUS "CMAKE_ERRFILE = ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log")
message(STATUS "CMAKE_LOGFILE = ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log")

check_type_size( ptrdiff_t    HAVE_PTRDIFF_T )
# How to check  if stdbool.h conforms to C99?

check_cxx_source_compiles(
    "
    #include <sys/dir.h>
    DIR *a = 0;
    int main(){return(0);}
    "
    HAVE_SYS_DIR_H)
AddHeaderTestToConfig("sys/dir.h" "HAVE_SYS_DIR_H")
check_cxx_source_compiles(
    "
    #include <ndir.h>
    DIR *a = 0;
    int main(){return(0);}
    "
    HAVE_NDIR_H)
AddHeaderTestToConfig("ndir.h" "HAVE_NDIR_H")
check_cxx_source_compiles(
    "
    #include <sys/ndir.h>
    DIR *a = 0;
    int main(){return(0);}
    "
    HAVE_SYS_NDIR_H)
AddHeaderTestToConfig("sys/ndir.h" "HAVE_SYS_NDIR_H")

check_include_files("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)
check_cxx_source_compiles(
  "
    #include <signal.h>

    static void sg(int b) { return;}
    int main(int c, char *av[]) { signal(0, sg); return(0);}
    "
  RETSIGTYPE_VOID)

if(RETSIGTYPE_VOID)
  set(RETSIGTYPE void)
else(RETSIGTYPE_VOID)
  set(RETSIGTYPE int)
endif(RETSIGTYPE_VOID)

test_big_endian(WORDS_BIGENDIAN)
check_cxx_source_compiles(
  "
    #include <sys/types.h>
    off_t a = 0;
    int main(){return(0);}
    "
  HAVE_OFF_T)

check_cxx_source_compiles(
  "
    #include <sys/types.h>
    size_t a = 0;
    int main(){return(0);}
    "
  HAVE_SIZE_T)
check_cxx_source_compiles(
  "
  #include <string.h>
  int main(){char b;char *a = strerror_r(0, &b, 0); return(0);}
  "
  STRERROR_R_CHAR_P)

SET(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
SET(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
check_cxx_source_compiles(
  "
    #include <string>
    #include <type_traits>
    #include <utility>
    
    int main()
    {
      std::string s(\"foo\");
      const std::string& const_s = std::as_const(s);
      return 0;
    }
  "  
  
  HAVE_AS_CONST
)

check_cxx_source_compiles(
  "
    #include <filesystem>
    #include <string>
    
    int main()
    {
      std::tr2::sys::path filepath;
      std::string s=filepath.filename().string();
      return 0;    
    };
  "
  HAVE_STD_FILESYSTEM
)

SET(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
  
# workaround for older boost versions (<1.55)
set(CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIR})
check_cxx_source_compiles(
  "
    #include <boost/config/suffix.hpp>
    #ifndef BOOST_FALLTHROUGH
      #error \"BOOST_FALLTHROUGH not defined\"
    #endif
    int main(){return(0);}
    "
    HAVE_BOOST_FALLTHROUGH
)
