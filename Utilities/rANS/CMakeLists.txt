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

o2_add_header_only_library(rANS
               TARGETVARNAME LIBRANS
               INTERFACE_LINK_LIBRARIES FairLogger::FairLogger Microsoft.GSL::GSL)
target_compile_options(${LIBRANS} INTERFACE -mavx2 -mfma -march=native)
if(${ENABLE_VTUNE_PROFILER})
  target_compile_definitions(${LIBRANS} INTERFACE ENABLE_VTUNE_PROFILER)
  target_link_libraries(${LIBRANS} INTERFACE PkgConfig::Vtune)
endif()

o2_add_test(Utils
            NAME Utils
            SOURCES test/test_ransUtils.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(Iterators
            NAME Iterators
            SOURCES test/test_ransIterators.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(Pack
            NAME Pack
            SOURCES test/test_ransPack.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(SIMDEncoderKernels
            NAME SIMDEncoderKernels
            SOURCES test/test_ransSIMDEncoderKernels.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(SIMD
            NAME SIMD
            SOURCES test/test_ransSIMD.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(AlignedArray
            NAME AlignedArray
            SOURCES test/test_ransAlignedArray.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)
  
o2_add_test(HistogramView
            NAME HistogramView
            SOURCES test/test_ransHistogramView.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(Histograms
            NAME Histograms
            TARGETVARNAME testHistograms
            SOURCES test/test_ransHistograms.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(SymbolTables
            NAME SymbolTables
            SOURCES test/test_ransSymbolTables.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(ReverseSymbolLookupTable
            NAME ReverseSymbolLookupTable
            SOURCES test/test_ransReverseSymbolLookupTable.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)

o2_add_test(EncodeDecode
            NAME EncodeDecode
            SOURCES test/test_ransEncodeDecode.cxx
            PUBLIC_LINK_LIBRARIES O2::rANS
            COMPONENT_NAME rANS
            LABELS utils)
            
if (TARGET benchmark::benchmark)
  o2_add_executable(CombinedIterator
                    SOURCES benchmarks/bench_ransCombinedIterator.cxx
                    COMPONENT_NAME rANS
                    IS_BENCHMARK
                    PUBLIC_LINK_LIBRARIES O2::rANS benchmark::benchmark)

  o2_add_executable(Encode
                    SOURCES benchmarks/bench_ransEncode.cxx
                    COMPONENT_NAME rANS
                    IS_BENCHMARK
                    TARGETVARNAME ransEncoderBenchmark
                    PUBLIC_LINK_LIBRARIES O2::rANS benchmark::benchmark RapidJSON::RapidJSON)
  target_compile_options(${ransEncoderBenchmark} PRIVATE  -O3 -march=native -fopenmp)
  target_link_libraries(${ransEncoderBenchmark} PRIVATE -ltbb)

  o2_add_executable(EncodeImpl
                    SOURCES benchmarks/bench_ransEncodeImpl.cxx
                    COMPONENT_NAME rANS
                    IS_BENCHMARK
                    TARGETVARNAME ransEncodeBenchmark
                    PUBLIC_LINK_LIBRARIES O2::rANS benchmark::benchmark)
  target_compile_options(${ransEncodeBenchmark} PRIVATE  -O3 -march=native -fopenmp)
  target_link_libraries(${ransEncodeBenchmark} PRIVATE -ltbb)

  o2_add_executable(Streaming
                  SOURCES benchmarks/bench_ransStreaming.cxx
                  COMPONENT_NAME rANS
                  IS_BENCHMARK
                  TARGETVARNAME ransStreamingBenchmark
                  PUBLIC_LINK_LIBRARIES O2::rANS benchmark::benchmark)
    target_compile_options(${ransStreamingBenchmark} PRIVATE  -O2 -march=native -fopenmp)
    target_link_libraries(${ransStreamingBenchmark} PRIVATE -ltbb)
    if(${ENABLE_VTUNE_PROFILER})
      target_compile_definitions(${ransStreamingBenchmark} PRIVATE ENABLE_VTUNE_PROFILER)
      target_link_libraries(${ransStreamingBenchmark} PRIVATE PkgConfig::Vtune)
    endif()

  o2_add_executable(TPCEncodeDecode
                    SOURCES benchmarks/bench_ransTPC.cxx
                    COMPONENT_NAME rANS
                    IS_BENCHMARK
                    PUBLIC_LINK_LIBRARIES O2::rANS 
                    TARGETVARNAME ransTPCBenchmark
                    PUBLIC_LINK_LIBRARIES O2::rANS Boost::program_options FairLogger::FairLogger RapidJSON::RapidJSON)
  target_compile_options(${ransTPCBenchmark} PRIVATE  -O3 -march=native -fopenmp)
  target_link_libraries(${ransTPCBenchmark} PRIVATE -ltbb)
  if(${ENABLE_VTUNE_PROFILER})
    target_compile_definitions(${ransTPCBenchmark} PRIVATE ENABLE_VTUNE_PROFILER)
    target_link_libraries(${ransTPCBenchmark} PRIVATE PkgConfig::Vtune)
  endif()

endif()

o2_add_executable(rans-encode-decode-8
                  TARGETVARNAME targetName
                  SOURCES run/bin-encode-decode.cxx
                  PUBLIC_LINK_LIBRARIES O2::rANS Boost::program_options)
target_compile_definitions(${targetName} PRIVATE -DSOURCE_T=uint8_t)

o2_add_executable(rans-encode-decode-16
                  TARGETVARNAME targetName
                  SOURCES run/bin-encode-decode.cxx
                  PUBLIC_LINK_LIBRARIES O2::rANS Boost::program_options)
target_compile_definitions(${targetName} PRIVATE -DSOURCE_T=uint16_t)

o2_add_executable(rans-encode-decode-32
                  TARGETVARNAME targetName
                  SOURCES run/bin-encode-decode.cxx
                  PUBLIC_LINK_LIBRARIES O2::rANS Boost::program_options)
target_compile_definitions(${targetName} PRIVATE -DSOURCE_T=uint32_t)