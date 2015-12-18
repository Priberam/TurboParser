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

#ifndef ENTITYFEATURES_H_
#define ENTITYFEATURES_H_

#include "SequenceFeatures.h"
#include "FeatureEncoder.h"
#include "EntityFeatureTemplates.h"
#include "EntityInstanceNumeric.h"
#include <bitset>

typedef std::vector<BinaryFeatures*> MultiBinaryFeatures;

class EntityFeatures : public SequenceFeatures {
public:
  EntityFeatures(Pipe* pipe) { pipe_ = pipe; }
  virtual ~EntityFeatures() { Clear(); }

  void AddUnigramFeatures(SequenceInstanceNumeric *sentence,
                          int position);

  void AddBigramFeatures(SequenceInstanceNumeric *sentence,
                         int position);

  void AddTrigramFeatures(SequenceInstanceNumeric *sentence,
                          int position);
protected:

  void ProcessFeature(uint64_t key,
                      BinaryFeatures * features,
                      BinaryFeatures * temp_popular_features,
                      std::bitset<NUM_SPECIAL_KEYS> * bitset) {
    if (FLAGS_train) {
      if (pipe_->GetParameters()->find_popular_keys == 0 &&
          pipe_->GetParameters()->process_popular_keys == 0) {
        AddFeature(key, features);
      }

      if (pipe_->GetParameters()->find_popular_keys == 1) {
        IncrementFeatureCounter(key);
        AddFeature(key, features);
      }

      if (pipe_->GetParameters()->process_popular_keys == 1) {
        //std::unordered_map < uint64_t, index_t> & popular_keys =
        //  pipe_->GetParameters()->popular_keys_index_in_bitset;
        //std::unordered_map < uint64_t, index_t>::const_iterator iterator =
        //  popular_keys.find(key);

        std::array<std::pair< key_64bits_t, index_t>, NUM_SPECIAL_KEYS>& popular_keys =
          pipe_->GetParameters()->popular_keys;
        int i = 0;
        for (; popular_keys[i].first != key && i < popular_keys.size(); i++) {}

        if (popular_keys[i].first == key && i < popular_keys.size()) {
          AddFeature(key, temp_popular_features);
          bitset->set(popular_keys[i].second);
        } else {
          AddFeature(key, features);
        }
      }
    }

    if (FLAGS_test) {
      /*   std::unordered_map < uint64_t, index_t> & popular_keys =
           pipe_->GetParameters()->popular_keys_index_in_bitset;*/
      std::array<std::pair< key_64bits_t, index_t>, NUM_SPECIAL_KEYS>& popular_keys =
        pipe_->GetParameters()->popular_keys;

      if (key & pipe_->GetParameters()->or_key) {
        //std::unordered_map < uint64_t, index_t>::const_iterator iterator =
        //  popular_keys.find(key);
        int i = 0;
        for (; popular_keys[i].first != key && i < popular_keys.size(); i++) {}

        if (popular_keys[i].first == key && i < popular_keys.size()) {
          AddFeature(key, temp_popular_features);
          bitset->set(popular_keys[i].second);
        } else {
          AddFeature(key, features);
        }
      } else {
        AddFeature(key, features);
      }
    }
  }

  void PostProcessFeatures(BinaryFeatures * features,
                           BinaryFeatures * temp_popular_features,
                           std::bitset<NUM_SPECIAL_KEYS> * bitset) {
    bitmap_64bits_t bitset_value = bitset->to_ullong();
    if (FLAGS_train) {
      if ((pipe_->GetParameters()->process_popular_keys) == 1) {
        pipe_->GetParameters()->add_unigram_features_calls += features->size() + temp_popular_features->size();

        CHECK(temp_popular_features->size() == bitset->count());
        if (bitset->count() >= 1) {
          pipe_->GetParameters()->add_unigram_features_calls_w_special_key += features->size() + 1;
          IncrementPopularKeyBitmapCounter(bitset_value);
        }
      }
    }
    if (FLAGS_test) {
      std::unordered_map<bitmap_64bits_t, key_64bits_t> & generated_keys_storage =
        pipe_->GetParameters()->popular_keybitmap_generated_key;
      std::unordered_map<bitmap_64bits_t, key_64bits_t>::const_iterator key_iterator =
        generated_keys_storage.find(bitset_value);
      if (key_iterator == generated_keys_storage.end()) {
        for (auto binfeat : *temp_popular_features) {
          AddFeature(binfeat, features);
        }
      } else {
        AddFeature((*key_iterator).second, features);
      }
    }
  }

  void AddFeature(key_64bits_t key,
                  BinaryFeatures * features) {
    features->push_back(key);
  }

  void IncrementFeatureCounter(key_64bits_t key) {
    Pipe * pipe = (pipe_);
    std::unordered_map< key_64bits_t, occurrences_counter_t> & counter =
      pipe->GetParameters()->featurekeys_occurrences_counter;
    std::unordered_map< key_64bits_t, occurrences_counter_t>::iterator iterator =
      counter.find(key);
    if (iterator != counter.end())
      ((*iterator).second)++;
    else
      counter.insert({ key, 1 });
  }

  void IncrementPopularKeyBitmapCounter(bitmap_64bits_t bitset_value) {
    Pipe * pipe = (pipe_);
    std::unordered_map<bitmap_64bits_t, occurrences_counter_t> & counter =
      pipe->GetParameters()->popular_keybitmap_occurrences_counter;
    std::unordered_map<bitmap_64bits_t, occurrences_counter_t>::iterator counter_iterator =
      counter.find(bitset_value);
    if (counter_iterator != counter.end())
      ((*counter_iterator).second)++;
    else
      counter.insert({ bitset_value, 1 });

    std::unordered_map<bitmap_64bits_t, key_64bits_t> & generated_keys_storage =
      pipe->GetParameters()->popular_keybitmap_generated_key;
    std::unordered_map<bitmap_64bits_t, key_64bits_t>::const_iterator key_iterator =
      generated_keys_storage.find(bitset_value);
    if (key_iterator == generated_keys_storage.end())
      generated_keys_storage.insert({ bitset_value, encoder_.CreateFKey_48B(KeyFeatureCombination,
                                                                  0xff,
                                                                  (key_64bits_t)bitset_value) });
  }
protected:
  // Encoder that converts features into a codeword.
  FeatureEncoder encoder_;
};

#endif /* ENTITYFEATURES_H_ */
