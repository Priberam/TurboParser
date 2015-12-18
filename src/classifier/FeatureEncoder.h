// Copyright (c) 2012-2015 Andre Martins
// All Rights Reserved.
//
// This file is part of TurboParser 2.3.
//
// TurboParser 2.3 is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TurboParser 2.3 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TurboParser 2.3.  If not, see <http://www.gnu.org/licenses/>.

#ifndef FEATURE_ENCODER_H_
#define FEATURE_ENCODER_H_

// This class implements several methods to pack a conjunction of atomic
// features into a 64-bit word.
class FeatureEncoder {
public:
  FeatureEncoder() {};
  virtual ~FeatureEncoder() {};

  // Several methods for forming a 64-bit word of features.
  // Argument "type" denotes the feature type (8 most significant bits);
  // "flags" may contain additional flags (8 less significant bits);
  // the remaining 48 bits can be formed by 16-bit or 8-bit words.
  // 16 bits are convenient for representing lexical features,
  // while 8 bits are usually convenient for POS tags and other non-lexical
  // features.
  // Some features may just have the "type" (8 most significant bits) and use
  // the remaining bits to flag whether a element of such feature was
  // verified/occured or not
  // (cardinality of bit flags is equal to such feature element cardinality).

  uint64_t AddHashMapKeyBit(uint64_t key) {
    return key;
  }

  uint64_t CreateFKey_NONE(uint8_t type,
                           uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_W(uint8_t type,
                        uint16_t w,
                        uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w) << 40);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WP(uint8_t type,
                         uint16_t w,
                         uint8_t p,
                         uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w) << 40);
    fkey |= (((uint64_t)p) << 32);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WPP(uint8_t type,
                          uint16_t w,
                          uint8_t p1, uint8_t p2,
                          uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w) << 40);
    fkey |= (((uint64_t)p1) << 32);
    fkey |= (((uint64_t)p2) << 24);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WPPP(uint8_t type,
                           uint16_t w,
                           uint8_t p1, uint8_t p2, uint8_t p3,
                           uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w) << 40);
    fkey |= (((uint64_t)p1) << 32);
    fkey |= (((uint64_t)p2) << 24);
    fkey |= (((uint64_t)p3) << 16);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WPPPP(uint8_t type,
                            uint16_t w,
                            uint8_t p1, uint8_t p2, uint8_t p3, uint8_t p4,
                            uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w) << 40);
    fkey |= (((uint64_t)p1) << 32);
    fkey |= (((uint64_t)p2) << 24);
    fkey |= (((uint64_t)p3) << 16);
    fkey |= (((uint64_t)p4) << 8);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WW(uint8_t type,
                         uint16_t w1, uint16_t w2,
                         uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w1) << 40);
    fkey |= (((uint64_t)w2) << 24);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WWW(uint8_t type,
                          uint16_t w1, uint16_t w2, uint16_t w3,
                          uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w1) << 40);
    fkey |= (((uint64_t)w2) << 24);
    fkey |= (((uint64_t)w3) << 8);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WWPP(uint8_t type,
                           uint16_t w1, uint16_t w2,
                           uint8_t p1, uint8_t p2,
                           uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w1) << 40);
    fkey |= (((uint64_t)w2) << 24);
    fkey |= (((uint64_t)p1) << 16);
    fkey |= (((uint64_t)p2) << 8);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_WWP(uint8_t type,
                          uint16_t w1, uint16_t w2,
                          uint8_t p,
                          uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)w1) << 40);
    fkey |= (((uint64_t)w2) << 24);
    fkey |= (((uint64_t)p) << 16);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_P(uint8_t type,
                        uint8_t p,
                        uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p) << 48);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_PP(uint8_t type,
                         uint8_t p1, uint8_t p2,
                         uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p1) << 48);
    fkey |= (((uint64_t)p2) << 40);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_PPP(uint8_t type,
                          uint8_t p1, uint8_t p2, uint8_t p3,
                          uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p1) << 48);
    fkey |= (((uint64_t)p2) << 40);
    fkey |= (((uint64_t)p3) << 32);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }
  uint64_t CreateFKey_PPPP(uint8_t type,
                           uint8_t p1, uint8_t p2, uint8_t p3, uint8_t p4,
                           uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p1) << 48);
    fkey |= (((uint64_t)p2) << 40);
    fkey |= (((uint64_t)p3) << 32);
    fkey |= (((uint64_t)p4) << 24);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_PPPPP(uint8_t type,
                            uint8_t p1, uint8_t p2,
                            uint8_t p3, uint8_t p4, uint8_t p5,
                            uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p1) << 48);
    fkey |= (((uint64_t)p2) << 40);
    fkey |= (((uint64_t)p3) << 32);
    fkey |= (((uint64_t)p4) << 24);
    fkey |= (((uint64_t)p5) << 16);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_PPPPPP(uint8_t type,
                             uint8_t p1, uint8_t p2, uint8_t p3,
                             uint8_t p4, uint8_t p5, uint8_t p6,
                             uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)p1) << 48);
    fkey |= (((uint64_t)p2) << 40);
    fkey |= (((uint64_t)p3) << 32);
    fkey |= (((uint64_t)p4) << 24);
    fkey |= (((uint64_t)p5) << 16);
    fkey |= (((uint64_t)p6) << 8);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_S(uint8_t type,
                        uint32_t s,
                        uint8_t flags) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)s) << 24);
    fkey |= ((uint64_t)flags);
    return AddHashMapKeyBit(fkey);
  }

  uint64_t CreateFKey_48B(uint8_t type, uint8_t flags, uint64_t b) {
    uint64_t fkey = (((uint64_t)type) << 56);
    fkey |= (((uint64_t)flags) << 48);
    fkey |= (((uint64_t)b) << 40);
    return fkey;
  }
};
#endif /* ENTITY_FEATURE_ENCODER_H_ */
