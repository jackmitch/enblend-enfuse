# Check if compiler supports C++11 features 
# and which compiler switches are necessary
# CXX11_FLAG : contains the necessary compiler flag

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(FindPackageHandleStandardArgs)

SET(CXX11_FLAG_CANDIDATES
  "--std=gnu++11"
  "--std=c++11"
  "--std=gnu++0x"
)

# sample openmp source code to test
SET(CXX11_TEST_SOURCE 
"
template <typename T>
struct check
{
    static_assert(sizeof(int) <= sizeof(T), \"not big enough\");
};

typedef check<check<bool>> right_angle_brackets;

class TestDeleted
{
public:
    TestDeleted() = delete;
};


int a;
decltype(a) b;

typedef check<int> check_type;
check_type c;
check_type&& cr = static_cast<check_type&&>(c);

auto d = a;

int main() {
  return 0;
};
")

# check c compiler
set(NUM 1)
FOREACH(FLAG ${CXX11_FLAG_CANDIDATES})
  SET(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  SET(CMAKE_REQUIRED_FLAGS "${FLAG}")
  CHECK_CXX_SOURCE_COMPILES("${CXX11_TEST_SOURCE}" CXX11_FLAG_DETECTED${NUM})
  SET(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
  IF(CXX11_FLAG_DETECTED${NUM})
    SET(CXX11_FLAG "${FLAG}")
    BREAK()
  ELSE()
    MATH(EXPR NUM "${NUM}+1")
  ENDIF(CXX11_FLAG_DETECTED${NUM})
ENDFOREACH()

# handle the standard arguments for find_package
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CXX11Compiler DEFAULT_MSG CXX11_FLAG)

MARK_AS_ADVANCED(CXX11_FLAG)
