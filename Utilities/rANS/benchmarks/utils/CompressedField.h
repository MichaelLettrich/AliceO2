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

/// @file   CompressedField.h
/// @author Michael Lettrich
/// @brief  common functionality for rANS benchmarks.

#ifndef RANS_BENCHMARKS_UTILS_COMPRESSEDFIELD_H_
#define RANS_BENCHMARKS_UTILS_COMPRESSEDFIELD_H_

#include <vector>
#include <concepts>
#include <span>
#include <cstdint>

#include "DetectorsCommonDataFormats/Metadata.h"

namespace o2::rans::benchmarks
{

namespace detail
{
// Source - https://stackoverflow.com/a/21028912
// Posted by Casey, modified by community. See post 'Timeline' for change history
// Retrieved 2025-12-19, License - CC BY-SA 3.0

// Allocator adaptor that interposes construct() calls to
// convert value initialization into default initialization.
template <typename T, typename A = std::allocator<T>>
class default_init_allocator : public A
{
  using a_t = std::allocator_traits<A>;

 public:
  template <typename U>
  struct rebind {
    using other =
      default_init_allocator<
        U, typename a_t::template rebind_alloc<U>>;
  };

  using A::A;

  template <typename U>
  void construct(U* ptr) noexcept(std::is_nothrow_default_constructible_v<U>)
  {
    ::new (static_cast<void*>(ptr)) U;
  }
  template <typename U, typename... Args>
  void construct(U* ptr, Args&&... args)
  {
    a_t::construct(static_cast<A&>(*this),
                   ptr, std::forward<Args>(args)...);
  }
};
} // namespace detail

template <std::integral source_T,
          std::integral stream_T = uint32_t,
          std::integral buffer_T = std::byte,
          class Allocator = std::allocator<std::byte>>
class CompressedField
{

 public:
  using source_type = source_T;
  using stream_type = stream_T;
  using buffer_type = buffer_T;

  using allocator_type = Allocator;

 private:
  using PayloadAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<stream_type>;
  using BufferAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<buffer_type>;

 public:
  explicit CompressedField(const allocator_type& alloc = allocator_type())
    : mPayload(alloc),
      mLiterals(alloc),
      mDict(alloc)
  {
  }

  CompressedField(size_t EncodeBufferSizeBytes, size_t literalsSizeBytes, size_t dictSizeBytes, const allocator_type& alloc = allocator_type())
    : mPayload(EncodeBufferSizeBytes / sizeof(stream_type), alloc),
      mLiterals(literalsSizeBytes / sizeof(buffer_type), alloc),
      mDict(dictSizeBytes / sizeof(buffer_type), alloc)
  {
  }

  CompressedField(const CompressedField& other, const allocator_type& alloc)
    : mPayload(other.mPayload, alloc),
      mLiterals(other.mLiterals, alloc),
      mDict(other.mDict, alloc)
  {
  }

  CompressedField(CompressedField&& other, const allocator_type& alloc)
    : mPayload(std::move(other.mPayload), alloc),
      mLiterals(std::move(other.mLiterals), alloc),
      mDict(std::move(other.mDict), alloc)
  {
  }

  CompressedField(const CompressedField&) = default;
  CompressedField(CompressedField&&) = default;
  CompressedField& operator=(const CompressedField&) = default;
  CompressedField& operator=(CompressedField&&) = default;

  allocator_type get_allocator() const { return allocator_type{}; }

  std::span<stream_type> getPayload() { return {mPayload}; }

  void commitPayload(stream_type* end_ptr)
  {
    mPayload.resize(static_cast<size_t>(std::distance(mPayload.data(), end_ptr)));
  }

  std::span<buffer_type> getLiterals() { return {mLiterals}; }

  void commitLiterals(buffer_type* end_ptr)
  {
    mLiterals.resize(static_cast<size_t>(std::distance(mLiterals.data(), end_ptr)));
  }

  std::span<buffer_type> getDict() { return {mDict}; }

  void commitDict(buffer_type* end_ptr)
  {
    mDict.resize(static_cast<size_t>(std::distance(mDict.data(), end_ptr)));
  }

  std::span<const stream_type> getPayload() const { return {mPayload}; }
  std::span<const buffer_type> getLiterals() const { return {mLiterals}; }
  std::span<const buffer_type> getDict() const { return {mDict}; }

  const o2::ctf::Metadata& getMetadata() const { return mMetadata; }

  o2::ctf::Metadata& getMetadata()
  {
    return const_cast<o2::ctf::Metadata&>(static_cast<const CompressedField&>(*this).getMetadata());
  }

  void resize(size_t encodeBufferSizeBytes, size_t literalsSizeBytes, size_t dictSizeBytes)
  {
    mPayload.resize(encodeBufferSizeBytes / sizeof(stream_type));
    mLiterals.resize(literalsSizeBytes / sizeof(buffer_type));
    mDict.resize(dictSizeBytes / sizeof(buffer_type));
  }

  void reset()
  {
    mPayload.clear();
    mLiterals.clear();
    mDict.clear();
  }

  void shrink_to_fit()
  {
    mPayload.shrink_to_fit();
    mLiterals.shrink_to_fit();
    mDict.shrink_to_fit();
  }

  size_t capacity() const
  {
    return mPayload.capacity() * sizeof(stream_type) +
           mLiterals.capacity() * sizeof(buffer_type) +
           mDict.capacity() * sizeof(buffer_type);
  }

  size_t size() const
  {
    return mPayload.size() * sizeof(stream_type) +
           mLiterals.size() * sizeof(buffer_type) +
           mDict.size() * sizeof(buffer_type);
  }

 private:
  std::vector<stream_type, PayloadAllocator> mPayload;
  std::vector<buffer_type, BufferAllocator> mLiterals;
  std::vector<buffer_type, BufferAllocator> mDict;
  o2::ctf::Metadata mMetadata;
};

} // namespace o2::rans::benchmarks

#endif /* RANS_BENCHMARKS_UTILS_COMPRESSEDFIELD_H_ */