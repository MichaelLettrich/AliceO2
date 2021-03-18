// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   kernels.h
/// @author Michael Lettrich
/// @since  2021-04-29
/// @brief

#ifndef RANS_INTERNAL_SIMD_KERNEL_H
#define RANS_INTERNAL_SIMD_KERNEL_H

#include "rANS/internal/backend/simd/sseKernel.h"
#include "rANS/internal/backend/simd/avxKernel.h"
#include "rANS/internal/backend/simd/avx512Kernel.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD_KERNEL_H */