// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   AttachableBuffer.h
/// @author Michael Lettrich
/// @brief Data buffer that can self manage or attach external data source transparently.

#ifndef RANS_INTERNAL_CONTAINERS_ATTACHABLEBUFFER_H_
#define RANS_INTERNAL_CONTAINERS_ATTACHABLEBUFFER_H_

#include <cstdint>
#include <cstring>
#include <array>

#include <gsl/span>

#include "rANS/internal/common/utils.h"

namespace o2::rans::internal
{

template <typename T>
class AttachableBuffer
{
 public:
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  AttachableBuffer(pointer ptr, size_type size) : mSize{size}, mPtr{ptr} {};
  AttachableBuffer(size_type size) : mSize{size}, mBuffer{std::make_unique<T[]>(size)}, mPtr{mBuffer.get()} {};
  AttachableBuffer() = default;
  ~AttachableBuffer() = default;
  AttachableBuffer(const AttachableBuffer& other);
  AttachableBuffer(AttachableBuffer&& other);

  AttachableBuffer& operator=(const AttachableBuffer& other);
  AttachableBuffer& operator=(AttachableBuffer&& other);

  [[nodiscard]] inline pointer data() const noexcept { return mPtr; };
  [[nodiscard]] inline pointer data() noexcept { return const_cast<pointer>(static_cast<const AttachableBuffer&>(*this).data()); };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return mPtr; };
  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };
  [[nodiscard]] inline iterator begin() noexcept { return const_cast<iterator>(static_cast<const AttachableBuffer&>(*this).begin()); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return cbegin() + size(); };
  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };
  [[nodiscard]] inline iterator end() noexcept { return const_cast<iterator>(static_cast<const AttachableBuffer&>(*this).end()); };

  [[nodiscard]] inline const_reverse_iterator crbegin() const noexcept { return std::reverse_iterator(cend()); };
  [[nodiscard]] inline const_reverse_iterator rbegin() const noexcept { return crbegin(); };
  [[nodiscard]] inline reverse_iterator rbegin() noexcept { return std::reverse_iterator(end()); };

  [[nodiscard]] inline const_reverse_iterator crend() const noexcept { return std::reverse_iterator(cbegin()); };
  [[nodiscard]] inline const_reverse_iterator rend() const noexcept { return crend(); };
  [[nodiscard]] inline reverse_iterator rend() noexcept { return std::reverse_iterator(begin()); };

  [[nodiscard]] inline size_type size() const noexcept { return mSize; };
  [[nodiscard]] inline bool empty() const noexcept { return size() == 0; };
  inline void reset() noexcept
  {
    AttachableBuffer<T> tmp{};
    swap(*this, tmp);
  };

  [[nodiscard]] inline const T& operator[](size_t idx) const
  {
    assert(idx < size());
    return mPtr[idx];
  };

  [[nodiscard]] inline T& operator[](size_t idx) { return const_cast<T&>(static_cast<const AttachableBuffer&>(*this)[idx]); };

  inline void assign(const T& value) { std::fill(begin(), end(), value); };
  void resize(size_type newSize);
  void attach(pointer ptr, size_type newSize);
  [[nodiscard]] std::unique_ptr<T[]> release() &&;

  inline friend void swap(AttachableBuffer& first, AttachableBuffer& second) noexcept
  {
    using std::swap;

    swap(first.mBuffer, second.mBuffer);
    swap(first.mPtr, second.mPtr);
    swap(first.mSize, second.mSize);
  };

 private:
  inline void updatePtr() noexcept { mPtr = mBuffer.get(); };

  size_type mSize{};
  std::unique_ptr<T[]> mBuffer{};
  pointer mPtr{mBuffer.get()};
};

template <typename T>
AttachableBuffer<T>::AttachableBuffer(const AttachableBuffer& other) : mSize{other.mSize}
{
  if (other.mBuffer) {
    this->mBuffer = std::make_unique<T[]>(this->size());
    this->updatePtr();
    std::copy(other.begin(), other.end(), this->begin());
  } else {
    this->mPtr = other.mPtr;
    assert(this->mBuffer == nullptr);
  }
};

template <typename T>
AttachableBuffer<T>& AttachableBuffer<T>::operator=(const AttachableBuffer<T>& other)
{
  auto tmp{other};
  swap(*this, tmp);
  return *this;
};

template <typename T>
AttachableBuffer<T>::AttachableBuffer(AttachableBuffer<T>&& other)
{
  swap(*this, other);
};

template <typename T>
AttachableBuffer<T>& AttachableBuffer<T>::operator=(AttachableBuffer<T>&& other)
{
  swap(*this, other);
  return *this;
};

} // namespace o2::rans::internal

namespace gsl
{
template <typename T>
[[nodiscard]] span<T> make_span(o2::rans::internal::AttachableBuffer<T>& buffer)
{
  return gsl::span(buffer.data(), buffer.size());
};

template <typename T>
[[nodiscard]] span<const T> make_span(const o2::rans::internal::AttachableBuffer<T>& buffer)
{
  return gsl::span(buffer.data(), buffer.size());
};

} // namespace gsl

namespace o2::rans::internal
{

template <typename T>
void AttachableBuffer<T>::resize(size_type newSize)
{
  AttachableBuffer<T> tmp(newSize);
  size_t copySize = tmp.size() < this->size() ? tmp.size() : this->size();
  std::copy(this->begin(), this->begin() + copySize, tmp.begin());
  swap(*this, tmp);
};

template <typename T>
[[nodiscard]] inline std::unique_ptr<T[]> AttachableBuffer<T>::release() &&
{
  auto ret = std::move(mBuffer);
  reset();
  return ret;
};

template <typename T>
inline void AttachableBuffer<T>::attach(pointer ptr, size_type newSize)
{
  AttachableBuffer<T> tmp{ptr, newSize};
  swap(*this, tmp);
};

} // namespace o2::rans::internal

#endif /* RANS_INTERNAL_CONTAINERS_ATTACHABLEBUFFER_H_ */