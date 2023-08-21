// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   typetraits.h
/// @author michael.lettrich@cern.ch
/// @brief  manipulation of types at compile time

#ifndef RANS_INTERNAL_COMMON_CONTAINERTRAITS_H_
#define RANS_INTERNAL_COMMON_CONTAINERTRAITS_H_

#include <type_traits>
#include <utility>

#include "rANS/internal/common/defaults.h"

namespace o2::rans
{

namespace internal
{
template <typename source_T, typename value_T>
class VectorContainer;

template <typename source_T, typename value_T>
class SparseVectorContainer;

template <typename source_T, typename value_T>
class HashContainer;

template <typename source_T, typename value_T>
class HashTable;

} // namespace internal

// forward declarations
template <typename source_T, typename>
class Histogram;

template <typename source_T>
class SparseHistogram;

template <typename source_T>
class HashHistogram;

template <typename container_T>
class RenormedHistogramImpl;

template <typename source_T, typename value_T>
class SymbolTable;

template <typename source_T, typename value_T>
class SparseSymbolTable;

template <typename source_T, typename value_T>
class HashSymbolTable;

template <typename source_T>
using RenormedHistogram = RenormedHistogramImpl<internal::VectorContainer<source_T, uint32_t>>;

template <typename source_T>
using RenormedSparseHistogram = RenormedHistogramImpl<internal::SparseVectorContainer<source_T, uint32_t>>;

template <typename source_T>
using RenormedHashHistogram = RenormedHistogramImpl<internal::HashContainer<source_T, uint32_t>>;

namespace internal
{

template <typename T>
struct isSymbolTable : std::false_type {
};

template <typename source_T, typename value_T>
struct isSymbolTable<SymbolTable<source_T, value_T>> : std::true_type {
};

template <typename source_T, typename value_T>
struct isSymbolTable<SparseSymbolTable<source_T, value_T>> : std::true_type {
};

template <typename source_T, typename value_T>
struct isSymbolTable<HashSymbolTable<source_T, value_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isSymbolTable_v = isSymbolTable<T>::value;

template <typename T>
struct isHistogram : std::false_type {
};

template <typename source_T>
struct isHistogram<Histogram<source_T, void>> : std::true_type {
};

template <typename source_T>
struct isHistogram<SparseHistogram<source_T>> : std::true_type {
};

template <typename source_T>
struct isHistogram<HashHistogram<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isHistogram_v = isHistogram<T>::value;

template <typename T>
struct isRenormedHistogram : std::false_type {
};

template <typename source_T>
struct isRenormedHistogram<RenormedHistogram<source_T>> : std::true_type {
};

template <typename source_T>
struct isRenormedHistogram<RenormedSparseHistogram<source_T>> : std::true_type {
};

template <typename source_T>
struct isRenormedHistogram<RenormedHashHistogram<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isRenormedHistogram_v = isRenormedHistogram<T>::value;

template <typename T>
struct isDenseContainer : std::false_type {
};
template <typename source_T>
struct isDenseContainer<Histogram<source_T, void>> : std::true_type {
};

template <typename source_T>
struct isDenseContainer<RenormedHistogram<source_T>> : std::true_type {
};

template <typename source_T, typename value_T>
struct isDenseContainer<SymbolTable<source_T, value_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isDenseContainer_v = isDenseContainer<T>::value;

template <typename T>
struct isSparseContainer : std::false_type {
};
template <typename source_T>
struct isSparseContainer<SparseHistogram<source_T>> : std::true_type {
};

template <typename source_T>
struct isSparseContainer<RenormedSparseHistogram<source_T>> : std::true_type {
};

template <typename source_T, typename value_T>
struct isSparseContainer<SparseSymbolTable<source_T, value_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isSparseContainer_v = isSparseContainer<T>::value;

template <typename T>
struct isHashContainer : std::false_type {
};

template <typename key_T, typename value_T>
struct isHashContainer<HashTable<key_T, value_T>> : std::true_type {
};

template <typename source_T>
struct isHashContainer<HashHistogram<source_T>> : std::true_type {
};

template <typename source_T>
struct isHashContainer<RenormedHashHistogram<source_T>> : std::true_type {
};

template <typename source_T, typename value_T>
struct isHashContainer<HashSymbolTable<source_T, value_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isHashContainer_v = isHashContainer<T>::value;

template <typename T>
struct isPair : public std::false_type {
};

template <typename A, typename B>
struct isPair<std::pair<A, B>> : public std::true_type {
};

template <typename T>
inline constexpr bool isPair_v = isPair<T>::value;

} // namespace internal
} // namespace o2::rans

#endif /* RANS_INTERNAL_COMMON_CONTAINERTRAITS_H_ */