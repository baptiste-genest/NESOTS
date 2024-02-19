if (TARGET nanothread)
  return()
endif()

include(FetchContent)

message(STATUS "Fetching nanothread")

FetchContent_Declare(
    nanothread
    GIT_REPOSITORY https://github.com/mitsuba-renderer/nanothread
    GIT_SHALLOW    TRUE
    )
FetchContent_MakeAvailable(nanothread)
