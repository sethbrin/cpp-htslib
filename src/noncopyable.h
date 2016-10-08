//
// Created by zp on 9/4/16.
//

#ifndef EASEHTSLIB_NONCOPYABLE_H
#define EASEHTSLIB_NONCOPYABLE_H


namespace ncic {

namespace easehts {

class NonCopyable {
 protected:
  NonCopyable() = default;
  ~NonCopyable() = default;
  NonCopyable( const NonCopyable& ) = delete;
  NonCopyable& operator=( const NonCopyable& ) = delete;
};

} // easehts
} // ncic
#endif //EASEHTSLIB_NONCOPYABLE_H
