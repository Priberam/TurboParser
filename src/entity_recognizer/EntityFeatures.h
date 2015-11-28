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

#include "Features.h"
#include "FeatureEncoder.h"
#include "EntityFeatureEncoder.h"
#include "EntityFeatureTemplates.h"
#include "EntityInstanceNumeric.h"

typedef std::vector<BinaryFeatures*> MultiBinaryFeatures;

class EntityFeatures : public Features {
public:
  EntityFeatures(Pipe* pipe) { pipe_ = pipe; }
  virtual ~EntityFeatures() { Clear(); }

public:
  void Clear() {
    for (int i = 0; i < input_features_unigrams_.size(); ++i) {
      if (!input_features_unigrams_[i]) continue;
      for (int j = 0; j < input_features_unigrams_[i]->size(); ++j) {
        if (!(*input_features_unigrams_[i])[j]) continue;
        (*input_features_unigrams_[i])[j]->clear();
        delete (*input_features_unigrams_[i])[j];
        (*input_features_unigrams_[i])[j] = NULL;
      }
      input_features_unigrams_[i]->clear();
      delete input_features_unigrams_[i];
      input_features_unigrams_[i] = NULL;
    }
    input_features_unigrams_.clear();

    for (int i = 0; i < input_features_bigrams_.size(); ++i) {
      if (!input_features_bigrams_[i]) continue;
      for (int j = 0; j < input_features_bigrams_[i]->size(); ++j) {
        if (!(*input_features_bigrams_[i])[j]) continue;
        (*input_features_bigrams_[i])[j]->clear();
        delete (*input_features_bigrams_[i])[j];
        (*input_features_bigrams_[i])[j] = NULL;
      }
      input_features_bigrams_[i]->clear();
      delete input_features_bigrams_[i];
      input_features_bigrams_[i] = NULL;
    }
    input_features_bigrams_.clear();

    for (int i = 0; i < input_features_trigrams_.size(); ++i) {
      if (!input_features_trigrams_[i]) continue;
      for (int j = 0; j < input_features_trigrams_[i]->size(); ++j) {
        if (!(*input_features_trigrams_[i])[j]) continue;
        (*input_features_trigrams_[i])[j]->clear();
        delete (*input_features_trigrams_[i])[j];
        (*input_features_trigrams_[i])[j] = NULL;
      }
      input_features_trigrams_[i]->clear();
      delete input_features_trigrams_[i];
      input_features_trigrams_[i] = NULL;
    }
    input_features_trigrams_.clear();
  }

  void Initialize(Instance *instance, Parts *parts) {
    Clear();
    int length = static_cast<EntityInstanceNumeric*>(instance)->size();
    input_features_unigrams_.resize(length,
                                    static_cast<MultiBinaryFeatures*>(NULL));
    input_features_bigrams_.resize(length + 1,
                                   static_cast<MultiBinaryFeatures*>(NULL));
    // Make this optional?
    input_features_trigrams_.resize(length + 1,
                                    static_cast<MultiBinaryFeatures*>(NULL));


    for (int i = 0; i < length; i++) {
      input_features_unigrams_[i] = new MultiBinaryFeatures;
      input_features_unigrams_[i]->resize(EntityFeatureTemplateUnigram::size,
                                          static_cast<BinaryFeatures*>(NULL));
    }
    for (int i = 0; i < length + 1; i++) {
      input_features_bigrams_[i] = new MultiBinaryFeatures;
      input_features_bigrams_[i]->resize(EntityFeatureTemplateBigram::size,
                                         static_cast<BinaryFeatures*>(NULL));
      // Make this optional?
	  input_features_trigrams_[i] = new MultiBinaryFeatures;
      input_features_trigrams_[i]->resize(EntityFeatureTemplateTrigram::size,
                                          static_cast<BinaryFeatures*>(NULL));
    }

  }

  const BinaryFeatures &GetPartFeatures(int r) const {
    CHECK(false) << "All part features are specific to unigrams, bigrams, "
      "or trigrams.";
    // Do this to avoid compilation error.
    return *new BinaryFeatures;
  };

  BinaryFeatures *GetMutablePartFeatures(int r) const {
    CHECK(false) << "All part features are specific to unigrams, bigrams, "
      "or trigrams.";
    return NULL;
  };

  const MultiBinaryFeatures &GetUnigramMultiFeatures(int i) const {
    return *(input_features_unigrams_[i]);
  };

  const BinaryFeatures &GetUnigramFeatures(int i, int j) const {
    return *(GetUnigramMultiFeatures(i))[j];
  };
  const BinaryFeatures &GetUnigramFeatures(const MultiBinaryFeatures &UnigramMultiFeatures,
                                           int j) const {
    return *(UnigramMultiFeatures[j]);
  };

  const MultiBinaryFeatures &GetBigramMultiFeatures(int i) const {
    return *(input_features_bigrams_[i]);
  };

  const BinaryFeatures &GetBigramFeatures(int i, int j) const {
    return *(GetBigramMultiFeatures(i))[j];
  };
  const BinaryFeatures &GetBigramFeatures(const MultiBinaryFeatures &BigramMultiFeatures,
                                           int j) const {
    return *(BigramMultiFeatures[j]);
  };

  const MultiBinaryFeatures &GetTrigramMultiFeatures(int i) const {
    return *(input_features_trigrams_[i]);
  };

  const BinaryFeatures &GetTrigramFeatures(int i, int j) const {
    return *(GetTrigramMultiFeatures(i))[j];
  };
  const BinaryFeatures &GetTrigramFeatures(const MultiBinaryFeatures &TrigramMultiFeatures,
                                           int j) const {
    return *(TrigramMultiFeatures[j]);
  };



public:
  void AddUnigramFeatures(EntityInstanceNumeric *sentence,
                          int position);

  void AddBigramFeatures(EntityInstanceNumeric *sentence,
                         int position);

  void AddTrigramFeatures(EntityInstanceNumeric *sentence,
                          int position);

protected:
  void AddFeature(uint64_t fkey,
                  MultiBinaryFeatures * multifeatures,
                  uint8_t familyfeature) {
    BinaryFeatures* features = (*multifeatures)[familyfeature];
    features->push_back(fkey);
  }

protected:
  // Encoder that converts features into a codeword.
  EntityFeatureEncoder encoder_;
  // Vectors of input features.
  vector<MultiBinaryFeatures*> input_features_unigrams_;
  vector<MultiBinaryFeatures*> input_features_bigrams_;
  vector<MultiBinaryFeatures*> input_features_trigrams_;
};

#endif /* ENTITYFEATURES_H_ */
