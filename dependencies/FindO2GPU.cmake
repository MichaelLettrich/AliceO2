# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

if(NOT DEFINED ENABLE_CUDA)
  set(ENABLE_CUDA "AUTO")
endif()
if(NOT DEFINED ENABLE_OPENCL)
  set(ENABLE_OPENCL "AUTO")
endif()
if(NOT DEFINED ENABLE_HIP)
  set(ENABLE_HIP "AUTO")
endif()
string(TOUPPER "${ENABLE_CUDA}" ENABLE_CUDA)
string(TOUPPER "${ENABLE_OPENCL}" ENABLE_OPENCL)
string(TOUPPER "${ENABLE_HIP}" ENABLE_HIP)
if(NOT DEFINED CMAKE_BUILD_TYPE_UPPER)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
endif()
if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
  set(GPUCA_BUILD_DEBUG 1)
endif()

if(CUDA_COMPUTETARGET AND CUDA_COMPUTETARGET STREQUAL "default")
  set(CUDA_COMPUTETARGET 86 89)
endif()

if(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET STREQUAL "default")
  set(HIP_AMDGPUTARGET gfx906;gfx908)
endif()

function(set_target_cuda_arch target)
  if(CUDA_COMPUTETARGET AND (CUDA_COMPUTETARGET MATCHES "86" OR CUDA_COMPUTETARGET MATCHES "89"))
    message(STATUS "Using optimized CUDA settings for Ampere GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_AMPERE)
  elseif(CUDA_COMPUTETARGET AND CUDA_COMPUTETARGET MATCHES "75")
    message(STATUS "Using optimized CUDA settings for Turing GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_TURING)
  else()
    message(STATUS "Defaulting optimized CUDA settings for Ampere GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_AMPERE)
  endif()
endfunction()

function(set_target_hip_arch target)
  if(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx906")
    message(STATUS "Using optimized HIP settings for MI50 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_VEGA)
  elseif(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx908")
    message(STATUS "Using optimized HIP settings for MI100 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_MI2xx)
  elseif(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx90a")
    message(STATUS "Using optimized HIP settings for MI210 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_MI2xx)
  else()
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_VEGA)
  endif()
endfunction()

# Detect and enable CUDA
STRING(REGEX REPLACE "\-std=[^ ]*" "" O2_GPU_CMAKE_CXX_FLAGS_NOSTD "${CMAKE_CXX_FLAGS}") # Need to strip c++17 imposed by alidist defaults

if(ENABLE_CUDA)
  set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
  set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
  include(CheckLanguage)
  check_language(CUDA)
  if (NOT ENABLE_CUDA STREQUAL "AUTO")
    if (NOT CMAKE_CUDA_COMPILER)
      set(CMAKE_CUDA_COMPILER "nvcc")
    endif()
    set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler")
  endif()
  if(CMAKE_CUDA_COMPILER)
    if(GPUCA_CUDA_GCCBIN)
      message(STATUS "Using as CUDA GCC version: ${GPUCA_CUDA_GCCBIN}")
      set(CMAKE_CUDA_HOST_COMPILER "${GPUCA_CUDA_GCCBIN}")
    endif()
    if(CUDA_COMPUTETARGET)
      set(CMAKE_CUDA_ARCHITECTURES ${CUDA_COMPUTETARGET} CACHE STRING "" FORCE)
    else()
      set(CMAKE_CUDA_ARCHITECTURES 61-virtual CACHE STRING "" FORCE)
    endif()
    enable_language(CUDA)
    get_property(LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
    if (ENABLE_CUDA STREQUAL "AUTO")
      set(FAILURE_SEVERITY STATUS)
    else()
      set(FAILURE_SEVERITY FATAL_ERROR)
    endif()
    if(NOT CUDA IN_LIST LANGUAGES)
      message(${FAILURE_SEVERITY} "CUDA was found but cannot be enabled")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
    find_path(THRUST_INCLUDE_DIR thrust/version.h PATHS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} NO_DEFAULT_PATH)
    if(THRUST_INCLUDE_DIR STREQUAL "THRUST_INCLUDE_DIR-NOTFOUND")
      message(${FAILURE_SEVERITY} "CUDA found but thrust not available")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
    if (NOT CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.6")
      message(${FAILURE_SEVERITY} "CUDA Version too old: ${CMAKE_CUDA_COMPILER_VERSION}, 12.6 required")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
  endif()
  if(CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_FLAGS "-Xcompiler \"${O2_GPU_CMAKE_CXX_FLAGS_NOSTD}\" ${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda -Xcompiler -Wno-attributes")
    if(GPUCA_KERNEL_RESOURCE_USAGE_VERBOSE)
      string(APPEND CMAKE_CUDA_FLAGS " -Xptxas -v")
    endif()
    string(APPEND CMAKE_CUDA_FLAGS " -Xcudafe --diag_suppress=114")
    if (NOT ENABLE_CUDA STREQUAL "AUTO")
      string(APPEND CMAKE_CUDA_FLAGS " --allow-unsupported-compiler")
    endif()
    set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "-Xcompiler \"${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}\" ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
    if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -lineinfo -Xptxas -O0 -Xcompiler -O0")
    else()
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -Xptxas -O4 -Xcompiler -O4")
    endif()
    set(GPUCA_CUDA_NO_FAST_MATH_FLAGS "--ftz=false --prec-div=true --prec-sqrt=true --fmad false")
    if(DEFINED GPUCA_NO_FAST_MATH AND "${GPUCA_NO_FAST_MATH}")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${GPUCA_CUDA_NO_FAST_MATH_FLAGS}")
    elseif(NOT CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -use_fast_math --ftz=true")#
    endif()
    if(CMAKE_CXX_FLAGS MATCHES "(^| )-Werror( |$)")
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Werror=cross-execution-space-call")
    endif()
    if(GPUCA_CUDA_GCCBIN)
      list(FILTER CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES EXCLUDE REGEX "^/usr/lib.*/gcc/") # Workaround, since CMake adds old GCC lib paths implicitly if we request that gcc for CUDA
    endif()

    set(CUDA_ENABLED ON)
    message(STATUS "CUDA found (Version ${CMAKE_CUDA_COMPILER_VERSION})")
  elseif(NOT ENABLE_CUDA STREQUAL "AUTO")
    message(FATAL_ERROR "CUDA not found (Compiler: ${CMAKE_CUDA_COMPILER})")
  else()
    set(CUDA_ENABLED OFF)
  endif()
endif()

# Detect and enable OpenCL 1.2 from AMD
if(ENABLE_OPENCL)
  find_package(OpenCL)
  if(ENABLE_OPENCL AND NOT ENABLE_OPENCL STREQUAL "AUTO")
    set_package_properties(OpenCL PROPERTIES TYPE REQUIRED)
  else()
    set_package_properties(OpenCL PROPERTIES TYPE OPTIONAL)
  endif()
endif()

# Detect and enable OpenCL 2.x
if(ENABLE_OPENCL)
  find_package(OpenCL)
  find_package(LLVM)
  if(LLVM_FOUND)
    find_package(Clang)
  endif()
  if (GPUCA_OPENCL_CLANGBIN)
    set(LLVM_CLANG ${GPUCA_OPENCL_CLANGBIN})
    execute_process(COMMAND "which" "${GPUCA_OPENCL_CLANGBIN}" OUTPUT_VARIABLE TMP_LLVM_SPIRV_PATH COMMAND_ERROR_IS_FATAL ANY)
    cmake_path(GET TMP_LLVM_SPIRV_PATH PARENT_PATH TMP_LLVM_SPIRV_PATH)
    find_program(LLVM_SPIRV llvm-spirv HINTS "${TMP_LLVM_SPIRV_PATH}")
  else()
    find_program(LLVM_CLANG clang HINTS "${Clang_DIR}/../../../bin-safe")
    find_program(LLVM_SPIRV llvm-spirv HINTS "${Clang_DIR}/../../../bin-safe")
  endif()
  if(Clang_FOUND
     AND LLVM_FOUND
     AND NOT LLVM_CLANG STREQUAL "LLVM_CLANG-NOTFOUND"
     AND LLVM_PACKAGE_VERSION VERSION_GREATER_EQUAL 13.0)
    set(OPENCL_COMPATIBLE_CLANG_FOUND ON)
  endif()
  if(OpenCL_VERSION_STRING VERSION_GREATER_EQUAL 2.2
     AND NOT LLVM_SPIRV STREQUAL "LLVM_SPIRV-NOTFOUND"
     AND OPENCL_COMPATIBLE_CLANG_FOUND)
    set(OPENCL_ENABLED_SPIRV ON)
    message(STATUS "Using CLANG ${LLVM_CLANG} and ${LLVM_SPIRV} for SPIR-V compilation")
  endif ()
  if(OPENCL_COMPATIBLE_CLANG_FOUND AND
     (OpenCL_VERSION_STRING VERSION_GREATER_EQUAL 2.2
     OR OPENCL_ENABLED_SPIRV))
    set(OPENCL_ENABLED ON)
    message(STATUS "Found OpenCL 2 (${OpenCL_VERSION_STRING} SPIR-V ${OPENCL_ENABLED_SPIRV} with CLANG ${LLVM_PACKAGE_VERSION})")
  elseif(NOT ENABLE_OPENCL STREQUAL "AUTO")
    message(FATAL_ERROR "OpenCL 2.x not available")
  else()
    set(OPENCL_ENABLED OFF)
  endif()
endif()

# Detect and enable HIP
if(ENABLE_HIP)
  if("$ENV{CMAKE_PREFIX_PATH}" MATCHES "rocm")
    set(CMAKE_HIP_STANDARD ${CMAKE_CXX_STANDARD})
    set(CMAKE_HIP_STANDARD_REQUIRED TRUE)
    if(HIP_AMDGPUTARGET)
      set(AMDGPU_TARGETS "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
      set(GPU_TARGETS "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
      set(CMAKE_HIP_ARCHITECTURES "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
    endif()
    set(TMP_ROCM_DIR_LIST $ENV{CMAKE_PREFIX_PATH})
    string(REPLACE ":" ";" TMP_ROCM_DIR_LIST "${TMP_ROCM_DIR_LIST}")
    list(FILTER TMP_ROCM_DIR_LIST INCLUDE REGEX rocm)
    list(POP_FRONT TMP_ROCM_DIR_LIST TMP_ROCM_DIR)
    get_filename_component(TMP_ROCM_DIR ${TMP_ROCM_DIR}/../../ ABSOLUTE)
    if (NOT DEFINED CMAKE_HIP_COMPILER)
      set(CMAKE_HIP_COMPILER "${TMP_ROCM_DIR}/llvm/bin/clang++")
      if(NOT EXISTS ${CMAKE_HIP_COMPILER})
        unset(CMAKE_HIP_COMPILER)
      endif()
    endif()
    include(CheckLanguage)
    check_language(HIP)
    find_package(hip)
    find_package(hipcub)
    find_package(rocprim)
    find_package(rocthrust)
    find_program(hip_HIPIFY_PERL_EXECUTABLE "hipify-perl")
    if(NOT hip_HIPIFY_PERL_EXECUTABLE)
      find_program(hip_HIPIFY_PERL_EXECUTABLE "hipify-perl" HINTS "${TMP_ROCM_DIR}/bin")
    endif()
    if(ENABLE_HIP STREQUAL "AUTO")
      set_package_properties(hip PROPERTIES TYPE OPTIONAL)
      set_package_properties(hipcub PROPERTIES TYPE OPTIONAL)
      set_package_properties(rocprim PROPERTIES TYPE OPTIONAL)
      set_package_properties(rocthrust PROPERTIES TYPE OPTIONAL)
    else()
      set_package_properties(hip PROPERTIES TYPE REQUIRED)
      set_package_properties(hipcub PROPERTIES TYPE REQUIRED)
      set_package_properties(rocprim PROPERTIES TYPE REQUIRED)
      set_package_properties(rocthrust PROPERTIES TYPE REQUIRED)
      if(NOT hip_HIPIFY_PERL_EXECUTABLE)
        message(FATAL_ERROR "Could not find hipify-perl command")
      endif()
    endif()
    if (CMAKE_HIP_COMPILER)
      enable_language(HIP)
      message(STATUS "HIP language enabled: ${CMAKE_HIP_COMPILER}")
    endif()
  elseif(NOT ENABLE_HIP STREQUAL "AUTO")
    message(FATAL_ERROR "HIP requested, but CMAKE_PREFIX_PATH env variable does not contain rocm folder!")
  endif()
  if(hip_FOUND AND NOT hip_VERSION VERSION_GREATER_EQUAL "5.5")
    set(hip_FOUND 0)
  endif()
  if(hip_FOUND AND hipcub_FOUND AND rocthrust_FOUND AND rocprim_FOUND AND hip_HIPCC_EXECUTABLE AND hip_HIPIFY_PERL_EXECUTABLE)
    set(HIP_ENABLED ON)
    set_target_properties(roc::rocthrust PROPERTIES IMPORTED_GLOBAL TRUE)
    message(STATUS "HIP Found (${hip_HIPCC_EXECUTABLE} version ${hip_VERSION})")
    set(O2_HIP_CMAKE_CXX_FLAGS "-fgpu-defer-diag -mllvm -amdgpu-enable-lower-module-lds=false -mllvm -amdgpu-function-calls=true -Wno-invalid-command-line-argument -Wno-unused-command-line-argument -Wno-invalid-constexpr -Wno-ignored-optimization-argument -Wno-unused-private-field -Wno-pass-failed")
    if(hip_VERSION VERSION_GREATER_EQUAL "6.0" AND NOT hip_VERSION VERSION_GREATER_EQUAL "6.2")
      string(APPEND O2_HIP_CMAKE_CXX_FLAGS " -mllvm -amdgpu-legacy-sgpr-spill-lowering=true") # TODO: Cleanup
    endif()
    if(GPUCA_KERNEL_RESOURCE_USAGE_VERBOSE)
      string(APPEND O2_HIP_CMAKE_CXX_FLAGS " -Rpass-analysis=kernel-resource-usage")
    endif()
    string(REGEX REPLACE "(gfx1[0-9]+;?)" "" CMAKE_HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}") # ROCm currently doesn’t support integrated graphics
    if(HIP_AMDGPUTARGET)
      set(CMAKE_HIP_ARCHITECTURES "${HIP_AMDGPUTARGET}") # If GPU build is enforced we override autodetection
    endif()
    if(NOT DEFINED GPUCA_NO_FAST_MATH OR NOT ${GPUCA_NO_FAST_MATH})
      string(APPEND O2_HIP_CMAKE_CXX_FLAGS " -fgpu-flush-denormals-to-zero -ffast-math")
    endif()
    set(CMAKE_HIP_FLAGS "${O2_GPU_CMAKE_CXX_FLAGS_NOSTD} ${CMAKE_HIP_FLAGS} ${O2_HIP_CMAKE_CXX_FLAGS}")
    set(CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
    if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
      set(CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -O0 -ggdb")
    else()
      set(CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER} "${CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -O3")
    endif()
  else()
    set(HIP_ENABLED OFF)
  endif()
  if(NOT HIP_ENABLED AND NOT ENABLE_HIP STREQUAL "AUTO")
    if (NOT hip_FOUND)
      message(WARNING "HIP not found")
    endif()
    if (NOT hipcub_FOUND)
      message(WARNING "hipcub not found")
    endif()
    if (NOT rocthrust_FOUND)
      message(WARNING "rocthrust not found")
    endif()
    if (NOT rocprim_FOUND)
      message(WARNING "rocprim not found")
    endif()
    if (NOT hip_HIPCC_EXECUTABLE)
      message(WARNING "hipcc executable not found")
    endif()
    if (NOT hip_HIPIFY_PERL_EXECUTABLE)
      message(WARNING "hipify-perl executable not found")
    endif()
    message(FATAL_ERROR "HIP requested but some of the above packages are not found")
  endif()

endif()

# if we end up here without a FATAL, it means we have found the "O2GPU" package
set(O2GPU_FOUND TRUE)
include("${CMAKE_CURRENT_LIST_DIR}/../GPU/GPUTracking/cmake/kernel_helpers.cmake")
