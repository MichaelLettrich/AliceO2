// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef O2_FRAMEWORK_SERVICEREGISTRYREF_H_
#define O2_FRAMEWORK_SERVICEREGISTRYREF_H_

#include "Framework/ServiceRegistry.h"
#include "Framework/Logger.h"

namespace o2::framework
{

class ServiceRegistryRef
{
 public:
  struct thread_safe_zone {
    // Cannot be copied / stored
    thread_safe_zone(thread_safe_zone const&) = delete;
    thread_safe_zone& operator=(thread_safe_zone const&) = delete;

    // Creating a thread safe zone will unlock
    // the registry, and re-lock it when it goes out of scope.
    // Use this to mark a section of code where we are
    // guaranteed to not be accessing the registry the current thread.
    thread_safe_zone(ServiceRegistryRef& ref)
      : mRef{ref}
    {
      mRef.unlock();
    }

    ~thread_safe_zone()
    {
      mRef.lock();
    }

    ServiceRegistryRef& mRef;
  };
  // The streamId is used to identify the stream in case we have multiple
  // threads. We cannot merely used the thread id because that does not
  // work in case the thread is created ad-hoc, like it appears to happen
  // for both libuv and FairMQ. This behaviour, BTW, makes usage
  // of thread local storage basically impossible (i.e. you lose state).
  // We use the following convention:
  // - streamId == 0 means the main thread
  // - streamId > 0 means one of the libuv worker threads
  // - streamId == -1 means the region callback thread
  // - streamId < -1 means some other worker thread of FairMQ which
  //  we do not know about.
  //
  // The getter will also make sure that a service of kind Stream
  // cannot be accessed if the streamId is <= 0 and complain accordingly.
  // The dataProcessorId will be used to distinguish between different
  // data processors when
  ServiceRegistryRef(ServiceRegistry& registry, ServiceRegistry::Salt salt = ServiceRegistry::globalDeviceSalt())
    : mRegistry(registry),
      mSalt(salt)
  {
  }

  // Wether or not this is the main thread
  bool isMainThread()
  {
    return mSalt.streamId == 0;
  }

  /// Check if service of type T is currently active.
  template <typename T>
    requires(std::is_const_v<T> == false)
  [[nodiscard]] bool active() const
  {
    return mRegistry.active<T>(mSalt);
  }

  /// Get a service for the given interface T. The returned reference exposed to
  /// the user is actually of the last concrete type C registered, however this
  /// should not be a problem.
  template <typename T>
  T& get() const
  {
    return mRegistry.get<T>(mSalt);
  }

  void registerService(ServiceTypeHash typeHash, void* service, ServiceKind kind, char const* name = nullptr) const
  {
    mRegistry.registerService(typeHash, service, kind, mSalt, name);
  }

  /// Register a service given an handle, notice how
  /// the service will be created in the current salt,
  /// so that from a dataprocessor you cannot create a service
  /// globally, or in a stream you cannot create services for
  /// a dataprocessor.
  void registerService(ServiceHandle handle)
  {
    mRegistry.registerService({handle.hash}, handle.instance, handle.kind, mSalt, handle.name.c_str());
  }

  void lock()
  {
    mRegistry.lock(mSalt);
  }

  void unlock()
  {
    mRegistry.unlock(mSalt);
  }

 private:
  ServiceRegistry& mRegistry;
  ServiceRegistry::Salt mSalt;
};

} // namespace o2::framework

#endif // O2_FRAMEWORK_SERVICEREGISTRY_H_
