set(DYNAMIC_LOADER_FILES
    dynamic_loader.h dynamic_loader.cc
    dynamic_loader_implementation.h dynamic_loader_implementation.cc
    win32_implementation.h sunny_implementation.h gmodule_implementation.h
)

if(WIN32)
  list(APPEND DYNAMIC_LOADER_FILES win32_implementation.cc)
else()
  if(HAVE_DL)
    list(APPEND DYNAMIC_LOADER_FILES sunny_implementation.cc)
  else()
    if(HAVE_GMODULE)
      list(APPEND DYNAMIC_LOADER_FILES gmodule_implementation.cc)
    else()
      list(APPEND DYNAMIC_LOADER_FILES null_implementation.h)
    endif()
  endif()
endif()
add_library(dynamic_loader STATIC ${DYNAMIC_LOADER_FILES})
