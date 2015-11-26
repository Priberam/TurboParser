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

typedef std::unordered_map<int, vector<double>> WordFeatureLabelScoresHashMap;

class WordFeatureLabelScoresCache {
public:
  WordFeatureLabelScoresCache() {
    hits_ = 0;
    misses_ = 0;
  };
  virtual ~WordFeatureLabelScoresCache() {};

  int hits() const { return hits_; };
  int misses() const { return misses_; };
  int GetSize() const { return cache_.size(); };

  void IncrementHits() { hits_ += 1; };
  void IncrementMisses() { misses_ += 1; };

  // Insert a new pair {key, value} in the hash-table.
  void Insert(int key, vector<double> value) {
    cache_.insert({ key, value });
  };

  // Searches for a given key in the hash-table.
  // If found, value is returned in argument 'value'.
  // return: true if found, false otherwise.
  bool Find(int key, vector<double> &value) {
    WordFeatureLabelScoresHashMap::const_iterator caching_iterator;
    caching_iterator = cache_.find(key);
    if (caching_iterator != cache_.end()) {
      value = caching_iterator->second;
      return true;
    };
    return false;
  };

protected:
  WordFeatureLabelScoresHashMap cache_;
  uint64_t hits_;
  uint64_t misses_;
};

class EntityFeatures : public SequenceFeatures {
public:
  EntityFeatures(Pipe* pipe) : SequenceFeatures(pipe) {}
  virtual ~EntityFeatures() {};

public:
  void Clear() {
    (*this).SequenceFeatures::Clear();
    for (int i = 0; i < input_features_cacheable_unigrams_.size(); ++i) {
      if (!input_features_cacheable_unigrams_[i]) continue;
      input_features_cacheable_unigrams_[i]->clear();
      delete input_features_cacheable_unigrams_[i];
      input_features_cacheable_unigrams_[i] = NULL;
    }
    input_features_cacheable_unigrams_.clear();

    for (int i = 0; i < input_features_non_cacheable_unigrams_.size(); ++i) {
      if (!input_features_non_cacheable_unigrams_[i]) continue;
      input_features_non_cacheable_unigrams_[i]->clear();
      delete input_features_non_cacheable_unigrams_[i];
      input_features_non_cacheable_unigrams_[i] = NULL;
    }
    input_features_non_cacheable_unigrams_.clear();
  }

  void Initialize(Instance *instance, Parts *parts) {
    Clear();
    (*this).SequenceFeatures::Initialize(instance, parts);
    int length = static_cast<SequenceInstanceNumeric*>(instance)->size();
    input_features_cacheable_unigrams_.resize(length, static_cast<BinaryFeatures*>(NULL));
    input_features_non_cacheable_unigrams_.resize(length, static_cast<BinaryFeatures*>(NULL));
  }

  const BinaryFeatures &GetUnigramCacheableFeatures(int i) const {
    return *(input_features_cacheable_unigrams_[i]);
  };

  const BinaryFeatures &GetUnigramNonCacheableFeatures(int i) const {
    return *(input_features_non_cacheable_unigrams_[i]);
  };

  void AddUnigramFeatures(SequenceInstanceNumeric *sentence,
                          int position);

  void AddBigramFeatures(SequenceInstanceNumeric *sentence,
                         int position);

  void AddTrigramFeatures(SequenceInstanceNumeric *sentence,
                          int position);

protected:
  void AddFeature(uint64_t fkey, BinaryFeatures* features) {
    features->push_back(fkey);
  }

protected:
  FeatureEncoder encoder_; // Encoder that converts features into a codeword.
  // Vectors of input features.
  vector<BinaryFeatures*> input_features_cacheable_unigrams_;
  vector<BinaryFeatures*> input_features_non_cacheable_unigrams_;
public:
  WordFeatureLabelScoresCache caching_wordfeatures_labelscores_;
};

#endif /* ENTITYFEATURES_H_ */
