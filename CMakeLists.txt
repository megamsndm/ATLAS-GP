cmake_minimum_required(VERSION 3.16)
project(atlas CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ── Release flags: -O3 + native CPU tuning ────────────────────────────────────
# -march=native tells GCC to use every instruction set available on this CPU
# (AVX2, FMA, etc.) — gives 10-20% speedup on matrix/Eigen code for free.
# -Wno-interference-size silences a noisy ABI warning from GCC 12+.
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG -march=native -funroll-loops -Wno-interference-size")

find_package(OpenMP REQUIRED)

# ── LBFGSpp ───────────────────────────────────────────────────────────────────
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/LBFGSpp/include/LBFGSB.h")
  message(FATAL_ERROR
    "LBFGSpp not found.\n"
    "Run ATLAS.bat — it clones LBFGSpp automatically.\n"
    "Or manually: git clone --depth=1 https://github.com/yixuan/LBFGSpp external/LBFGSpp")
endif()
set(LBFGS_INC "${CMAKE_SOURCE_DIR}/external/LBFGSpp/include")

# ── Eigen: search in order, download as last resort ───────────────────────────
set(EIGEN3_INCLUDE_DIR "")

# 1. Bundled inside LBFGSpp (most common — LBFGSpp ships Eigen headers)
if(EXISTS "${LBFGS_INC}/Eigen/Dense")
  set(EIGEN3_INCLUDE_DIR "${LBFGS_INC}")
  message(STATUS "Eigen: bundled inside LBFGSpp")
endif()

# 2. Previously downloaded to external/eigen
if(NOT EIGEN3_INCLUDE_DIR)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external/eigen/Eigen/Dense")
    set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external/eigen")
    message(STATUS "Eigen: ${EIGEN3_INCLUDE_DIR}")
  endif()
endif()

# 3. System-installed Eigen (Linux / MSYS2 pacman eigen3)
if(NOT EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 QUIET)
  if(Eigen3_FOUND)
    set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIRS}")
    message(STATUS "Eigen: system (${EIGEN3_INCLUDE_DIR})")
  endif()
endif()

# 4. Download from GitLab as last resort
if(NOT EIGEN3_INCLUDE_DIR)
  message(STATUS "Eigen not found locally — downloading Eigen 3.4.0 ...")
  set(EIGEN_TAR "${CMAKE_SOURCE_DIR}/external/eigen.tar.gz")
  set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/external/eigen")
  file(DOWNLOAD
    "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
    "${EIGEN_TAR}"
    SHOW_PROGRESS
    TLS_VERIFY ON
    STATUS DL_STATUS)
  list(GET DL_STATUS 0 DL_CODE)
  list(GET DL_STATUS 1 DL_MSG)
  if(DL_CODE EQUAL 0)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar xzf "${EIGEN_TAR}"
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/external")
    file(GLOB _ED "${CMAKE_SOURCE_DIR}/external/eigen-3*")
    if(_ED)
      file(RENAME "${_ED}" "${EIGEN_DIR}")
    endif()
  else()
    message(WARNING "Eigen download failed: ${DL_MSG}")
  endif()

  if(EXISTS "${EIGEN_DIR}/Eigen/Dense")
    set(EIGEN3_INCLUDE_DIR "${EIGEN_DIR}")
    message(STATUS "Eigen: downloaded to ${EIGEN_DIR}")
  else()
    message(FATAL_ERROR
      "Eigen not found and download failed.\n"
      "Manual fix: download\n"
      "  https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz\n"
      "and extract so that  external/eigen/Eigen/Dense  exists.\n"
      "Or run:  pacman -S mingw-w64-x86_64-eigen3  in MSYS2.")
  endif()
endif()

# ── Build ──────────────────────────────────────────────────────────────────────
add_executable(atlas "${CMAKE_SOURCE_DIR}/src/main.cpp")

target_include_directories(atlas PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/src/potential"
    "${LBFGS_INC}"
    "${EIGEN3_INCLUDE_DIR}"
)

target_link_libraries(atlas PRIVATE OpenMP::OpenMP_CXX)

if(MINGW)
  target_link_libraries(atlas PRIVATE -lmingwex)
endif()

target_compile_definitions(atlas PRIVATE
    NOMINMAX          # prevents Windows.h from defining min/max macros
    _USE_MATH_DEFINES # provides M_PI, M_E etc. on Windows
)

message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "Eigen      : ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Compiler   : ${CMAKE_CXX_COMPILER}")
