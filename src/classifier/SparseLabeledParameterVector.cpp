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

#include "SparseLabeledParameterVector.h"
#include "SerializationUtils.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <glog/logging.h>
#include <bitset>

using namespace std;

bool SparseLabelWeights::IsSparse() const { return true; }
int SparseLabelWeights::Size() const { return label_weights_.size(); }

double SparseLabelWeights::GetWeight(int label) const {
  for (int k = 0; k < label_weights_.size(); ++k) {
    if (label == label_weights_[k].first) {
      return label_weights_[k].second;
    }
  }
  return 0.0;
}
void SparseLabelWeights::SetWeight(int label, double weight) {
  for (int k = 0; k < label_weights_.size(); ++k) {
    if (label == label_weights_[k].first) {
      label_weights_[k].second = weight;
      return;
    }
  }
  label_weights_.push_back(std::pair<int, double>(label, weight));
}
void SparseLabelWeights::AddWeight(int label, double weight) {
  for (int k = 0; k < label_weights_.size(); ++k) {
    if (label == label_weights_[k].first) {
      label_weights_[k].second += weight;
      return;
    }
  }
  label_weights_.push_back(std::pair<int, double>(label, weight));
}

double SparseLabelWeights::SetWeightAndNormalize(int label,
                                                 double value,
                                                 double scale_factor) {
  double weight = value / scale_factor;
  double previous_value;
  for (int k = 0; k < label_weights_.size(); ++k) {
    if (label == label_weights_[k].first) {
      previous_value = label_weights_[k].second * scale_factor;
      label_weights_[k].second = weight;
      return previous_value;
    }
  }
  label_weights_.push_back(std::pair<int, double>(label, weight));
  return 0.0;
}

double SparseLabelWeights::AddWeightAndNormalize(int label,
                                                 double value,
                                                 double scale_factor) {
  double weight = value / scale_factor;
  double previous_value;
  for (int k = 0; k < label_weights_.size(); ++k) {
    if (label == label_weights_[k].first) {
      previous_value = label_weights_[k].second * scale_factor;
      label_weights_[k].second += weight;
      return previous_value;
    }
  }
  label_weights_.push_back(std::pair<int, double>(label, weight));
  return 0.0;
}

void SparseLabelWeights::GetLabelWeightByPosition(int position, int *label,
                                                  double *weight) const {
  *label = label_weights_[position].first;
  *weight = label_weights_[position].second;
  CHECK_GE(*label, 0);
}

void SparseLabelWeights::SetWeightByPosition(int position, double weight) {
  label_weights_[position].second = weight;
}

DenseLabelWeights::DenseLabelWeights(LabelWeights *label_weights) {
  CHECK(label_weights->IsSparse());
  for (int k = 0; k < label_weights->Size(); ++k) {
    int label;
    double weight;
    label_weights->GetLabelWeightByPosition(k, &label, &weight);
    CHECK_GE(label, 0);
    SetWeight(label, weight);
  }
}

bool DenseLabelWeights::IsSparse() const { return false; }
int DenseLabelWeights::Size() const { return weights_.size(); }

double DenseLabelWeights::GetWeight(int label) const {
  if (label >= weights_.size()) return 0.0;
  return weights_[label];
}
void DenseLabelWeights::SetWeight(int label, double weight) {
  CHECK_GE(label, 0);
  if (label >= weights_.size()) {
    weights_.resize(label + 1, 0.0);
  }
  weights_[label] = weight;
}
void DenseLabelWeights::AddWeight(int label, double weight) {
  CHECK_GE(label, 0);
  if (label >= weights_.size()) {
    weights_.resize(label + 1, 0.0);
  }
  weights_[label] += weight;
}

double DenseLabelWeights::SetWeightAndNormalize(int label,
                                                double value,
                                                double scaling_factor) {
  CHECK_GE(label, 0);
  if (label >= weights_.size()) {
    weights_.resize(label + 1, 0.0);
  }
  double weight = value / scaling_factor;
  double previous_value = weights_[label] * scaling_factor;
  weights_[label] = weight;
  return previous_value;
}

double DenseLabelWeights::AddWeightAndNormalize(int label,
                                                double value,
                                                double scaling_factor) {
  CHECK_GE(label, 0);
  if (label >= weights_.size()) {
    weights_.resize(label + 1, 0.0);
  }
  double weight = value / scaling_factor;
  double previous_value = weights_[label] * scaling_factor;
  weights_[label] += weight;
  return previous_value;
}

void DenseLabelWeights::GetLabelWeightByPosition(int position, int *label,
                                                 double *weight) const {
  CHECK_GE(position, 0);
  *label = position;
  *weight = weights_[position];
}

void DenseLabelWeights::SetWeightByPosition(int position, double weight) {
  weights_[position] = weight;
}

// Lock/unlock the parameter vector. If the vector is locked, no new features
// can be inserted.
void SparseLabeledParameterVector::StopGrowth() { growth_stopped_ = true; }
void SparseLabeledParameterVector::AllowGrowth() { growth_stopped_ = false; }
bool SparseLabeledParameterVector::growth_stopped() const { return growth_stopped_; }

// Clear the parameter vector.
void SparseLabeledParameterVector::Clear() {
  for (LabeledParameterMap::iterator iterator = map_values_.begin();
  iterator != map_values_.end();
    ++iterator) {
    delete iterator->second;
  }
  map_values_.clear();

  for (auto &family_feature : matrix_values_) {
    for (auto &feature_key : family_feature) {
      delete feature_key;
    }
    family_feature.clear();
  }
  matrix_values_.clear();
}

// Save/load the parameters to/from a file.
void SparseLabeledParameterVector::Save(FILE *fs) const {
  bool success;
  success = WriteInteger(fs, GetMapSize());
  CHECK(success);
  for (LabeledParameterMap::const_iterator iterator = map_values_.begin();
  iterator != map_values_.end();
    ++iterator) {
    success = WriteUINT64(fs, iterator->first);
    CHECK(success);
    const LabelWeights *label_weights = iterator->second;
    int length = label_weights->Size();
    success = WriteInteger(fs, length);
    CHECK(success);
    int label;
    double value;
    for (int k = 0; k < length; ++k) {
      label_weights->GetLabelWeightByPosition(k, &label, &value);
      CHECK_GE(label, 0);
      success = WriteInteger(fs, label);
      CHECK(success);
      success = WriteDouble(fs, value);
      CHECK(success);
    }
  }

  success = WriteInteger(fs, GetMatrixFamilyFeatureVectorSize());
  CHECK(success);
  for (vector<vector<LabelWeights*>>::const_iterator
       ptr_to_feature_key_vector = matrix_values_.begin();
       ptr_to_feature_key_vector != matrix_values_.end();
       ptr_to_feature_key_vector++) {
    success = WriteInteger(fs,
                           GetMatrixFeatureKeyVectorSize(ptr_to_feature_key_vector));
    CHECK(success);
    for (vector<LabelWeights*>::const_iterator
         ptr_to_label_weights_vector = ptr_to_feature_key_vector->begin();
         ptr_to_label_weights_vector != ptr_to_feature_key_vector->end();
         ptr_to_label_weights_vector++) {
      const LabelWeights *label_weights = *ptr_to_label_weights_vector;
      int length = label_weights->Size();
      success = WriteInteger(fs, length);
      CHECK(success);
      int label;
      double value;
      for (int k = 0; k < length; ++k) {
        label_weights->GetLabelWeightByPosition(k, &label, &value);
        CHECK_GE(label, 0);
        success = WriteInteger(fs, label);
        CHECK(success);
        success = WriteDouble(fs, value);
        CHECK(success);
      }
    }
  }
}

void SparseLabeledParameterVector::Load(FILE *fs) {
  bool success;
  int num_features;

  Initialize();
  success = ReadInteger(fs, &num_features);
  CHECK(success);
  for (int i = 0; i < num_features; ++i) {
    uint64_t key;
    success = ReadUINT64(fs, &key);
    CHECK(success);
    int length;
    success = ReadInteger(fs, &length);
    CHECK(success);
    int label;
    double value;
    for (int k = 0; k < length; ++k) {
      success = ReadInteger(fs, &label);
      CHECK(success);
      success = ReadDouble(fs, &value);
      CHECK(success);
      Set(key, label, value);
    }
  }

  int matrix_size_x;
  success = ReadInteger(fs, &matrix_size_x);
  CHECK(success);
  matrix_values_.resize(matrix_size_x);
  for (vector<vector<LabelWeights*>>::iterator
       ptr_to_feature_key_vector = matrix_values_.begin();
       ptr_to_feature_key_vector != matrix_values_.end();
       ptr_to_feature_key_vector++) {
    int matrix_size_y;
    success = ReadInteger(fs, &matrix_size_y);
    CHECK(success);
    ResizeFeatureKeyVector(ptr_to_feature_key_vector, matrix_size_y);
    for (vector<LabelWeights*>::iterator
         ptr_to_label_weights_vector = ptr_to_feature_key_vector->begin();
         ptr_to_label_weights_vector != ptr_to_feature_key_vector->end();
         ptr_to_label_weights_vector++) {
      int length;
      success = ReadInteger(fs, &length);
      CHECK(success);
      int label;
      double value;
      for (int k = 0; k < length; ++k) {
        success = ReadInteger(fs, &label);
        CHECK(success);
        success = ReadDouble(fs, &value);
        CHECK(success);
        SetValue(ptr_to_label_weights_vector, label, value);
      }
    }
  }

#define PRINT_STATISTICS
#ifdef PRINT_STATISTICS
  // Print some statistics:
  int num_sparse = 0;
  int num_total = 0;
  int num_labels_sparse = 0;
  for (LabeledParameterMap::iterator iterator = map_values_.begin();
  iterator != map_values_.end();
    ++iterator) {
    LabelWeights *label_weights = iterator->second;
    if (label_weights->IsSparse()) {
      ++num_sparse;
      int length = label_weights->Size();
      num_labels_sparse += length;
    }
    ++num_total;
  }
  LOG(INFO) << "Statistics for labeled parameter vector:";
  LOG(INFO) << "Features with sparse labels: " << num_sparse
    << " Total: " << num_total
    << " Sparse labels: " << num_labels_sparse;

  LOG(INFO) << "ALERT: Statistics for matrix-based values not included.";
#endif
}

void SparseLabeledParameterVector::Initialize() {
  Clear();
  scale_factor_ = 1.0;
  squared_norm_ = 0.0;
}

uint64_t SparseLabeledParameterVector::Size() const {
  return GetMapSize() + GetMatrixSize();
}
uint64_t SparseLabeledParameterVector::GetMapSize() const {
  return map_values_.size();
}
uint64_t SparseLabeledParameterVector::GetMatrixSize() const {
  return 0;
  LOG(INFO) << "ALERT: Size statistics for matrix-based values are not implemented.";
} //NOT IMPLEMENTED

uint64_t SparseLabeledParameterVector::GetMatrixFamilyFeatureVectorSize() const {
  return matrix_values_.size();
}
uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(
  int x) const {
  return matrix_values_[x].size();
}
uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(
  const vector<LabelWeights*> & feature_key_vector) const {
  return feature_key_vector.size();
}
uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(
  std::vector<std::vector<LabelWeights*>>::const_iterator
  ptr_to_feature_key_vector) const {
  return ptr_to_feature_key_vector->size();
}

uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(
  int x,
  int y) const {
  return ((matrix_values_[x])[y])->Size();
}
uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(
  LabelWeights* label_weights_vector) const {
  return label_weights_vector->Size();
}
uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(
  std::vector<LabelWeights*>::const_iterator ptr_to_label_weights_vector) const {
  return (*ptr_to_label_weights_vector)->Size();
}

bool SparseLabeledParameterVector::isHashMapKey(uint64_t key) const {
  uint64_t mask_matrixmap = ((uint64_t)0x8000000000000000);
  if ((key & mask_matrixmap) == mask_matrixmap) {
    return true;
  } else {
    return false;
  }
}
bool SparseLabeledParameterVector::isMatrixMapKey(uint64_t key) const {
  if (isHashMapKey(key)) {
    return false;
  } else {
    return true;
  }
}
bool SparseLabeledParameterVector::isMulti64bitKey(uint64_t key) const {
  uint64_t mask_multi64bit = ((uint64_t)0x0000000000000003);
  if (isHashMapKey(key)) {
    return false;
  } else {
    if ((key & mask_multi64bit) == mask_multi64bit) {
      return true;
    } else {
      return false;
    }
  }
}
uint64_t SparseLabeledParameterVector::GetFamilyFeatureKey(uint64_t key) const {
  uint64_t mask = ((uint64_t)0xFF00000000000000);
  uint64_t  masked_key = (key & mask) >> 56;
  return    masked_key;
}
uint64_t SparseLabeledParameterVector::GetBitKeys(uint64_t key) const {
  uint64_t  mask_obtain_bitkeys = ((uint64_t)0x00ffffffffffff00);
  uint64_t  masked_key = (key & mask_obtain_bitkeys) >> 8;
  return    masked_key;
}
uint64_t SparseLabeledParameterVector::Get3_5BitKeys(uint64_t key) const {
  uint64_t  mask_obtain_bitkeys = ((uint64_t)0x000000000000001c);
  uint64_t  masked_key = (key & mask_obtain_bitkeys) >> 2;
  return    masked_key;
}
uint64_t SparseLabeledParameterVector::Get6_8BitKeys(uint64_t key) const {
  uint64_t  mask_obtain_bitkeys = ((uint64_t)0x00000000000000e0);
  uint64_t  masked_key = (key & mask_obtain_bitkeys) >> 5;
  return    masked_key;
}

bool SparseLabeledParameterVector::Exists(uint64_t key) const {
  if (isHashMapKey(key)) {
    LabeledParameterMap::const_iterator iterator = map_values_.find(key);
    if (iterator == map_values_.end()) return false;
    return true;
  } else {
    uint64_t feature_family_index = GetFamilyFeatureKey(key);
    if (GetMatrixFamilyFeatureVectorSize() < feature_family_index) {
      return false;
    } else {
      return true;
    }
  }
}

bool SparseLabeledParameterVector::Get(uint64_t key, const vector<int> &labels,
                                       vector<double> *weights) const {
  if (isHashMapKey(key)) {
    LabeledParameterMap::const_iterator iterator = map_values_.find(key);
    if (iterator == map_values_.end()) {
      weights->clear();
      return false;
    }
    GetValues(iterator, labels, weights);
    return true;
  } else {
    if (!Exists(key)) {
      weights->clear();
      return false;
    }

    weights->clear();

    std::vector<std::vector<LabelWeights*>>::const_iterator
      feature_key_vector = FindFamilyFeatureVector(key);
    if (feature_key_vector == matrix_values_.end()) {
      weights->clear();
      return false;
    }
    // uint64_t num_64bitword_keys = Get3_5BitKeys(key);
    uint64_t current_64bitword_key = Get6_8BitKeys(key);
    bitset<48> bitset_keys = bitset<48>(GetBitKeys(key)); //uint64_t bitset_keys = GetBitKeys(key);
    for (int i = 0; i < 48; i++) {
      if (bitset_keys.test(i)) { //if ((bitset_keys & (1i64 << i)) != 0){
        std::vector<LabelWeights*>::const_iterator label_weights =
          FindFeatureKeyVector(feature_key_vector, i + 48 * current_64bitword_key);
        if (label_weights == (*feature_key_vector).end()) {
          continue;
        }
        GetValuesAndUpdate(label_weights, labels, weights);
      }
    }
    return true;
  }
}

double SparseLabeledParameterVector::GetSquaredNorm() const {
  return squared_norm_;
}

void SparseLabeledParameterVector::Scale(double scale_factor) {
  scale_factor_ *= scale_factor;
  squared_norm_ *= scale_factor * scale_factor;
  RenormalizeIfNecessary();
}

bool SparseLabeledParameterVector::Set(uint64_t key,
                                       int label,
                                       double value) {
  CHECK_GE(label, 0);

  if (isHashMapKey(key)) {
    LabeledParameterMap::iterator iterator = FindOrInsert(key);
    if (iterator != map_values_.end()) {
      SetValue(iterator, label, value);
      return true;
    } else {
      return false;
    }
  } else {
    std::vector<std::vector<LabelWeights*>>::iterator
      feature_key_vector = FindOrResizeFamilyFeatureVector(key);
    uint64_t num_64bitword_keys = Get3_5BitKeys(key);
    uint64_t current_64bitword_key = Get6_8BitKeys(key);
    bitset<48> bitset_keys = bitset<48>(GetBitKeys(key));
    for (int i = 0; i < 48; i++) {
      if (bitset_keys.test(i)) {
        std::vector<LabelWeights*>::iterator label_weights =
          FindOrResizeFeatureKeyVector(feature_key_vector, i + 48 * current_64bitword_key);
        SetValue(label_weights, label, value);
      }
    }
    return true;
  }
}

bool SparseLabeledParameterVector::Add(uint64_t key,
                                       int label,
                                       double value) {
  if (isHashMapKey(key)) {
    LabeledParameterMap::iterator iterator = FindOrInsert(key);
    if (iterator != map_values_.end()) {
      AddValue(iterator, label, value);
      return true;
    } else {
      return false;
    }
  } else {
    std::vector<std::vector<LabelWeights*>>::iterator
      feature_key_vector = FindOrResizeFamilyFeatureVector(key);
    uint64_t num_64bitword_keys = Get3_5BitKeys(key);
    uint64_t current_64bitword_key = Get6_8BitKeys(key);
    bitset<48> bitset_keys = bitset<48>(GetBitKeys(key));
    for (uint64_t i = 0; i < 48; i++) {
      if (bitset_keys.test(i)) {
        std::vector<LabelWeights*>::iterator label_weights =
          FindOrResizeFeatureKeyVector(feature_key_vector, i + 48 * current_64bitword_key);
        AddValue(label_weights, label, value);
      }
    }
    return true;
  }
}

void SparseLabeledParameterVector::Add(uint64_t key,
                                       const vector<int> &labels,
                                       const vector<double> &values) {
  if (isHashMapKey(key)) {
    LabeledParameterMap::iterator iterator = FindOrInsert(key);
    for (int k = 0; k < labels.size(); ++k) {
      AddValue(iterator, labels[k], values[k]);
    }
  } else {
    std::vector<std::vector<LabelWeights*>>::iterator
      feature_key_vector = FindOrResizeFamilyFeatureVector(key);
    uint64_t num_64bitword_keys = Get3_5BitKeys(key);
    uint64_t current_64bitword_key = Get6_8BitKeys(key);
    bitset<48> bitset_keys = bitset<48>(GetBitKeys(key));
    for (int i = 0; i < 48; i++) {
      if (bitset_keys.test(i)) {
        std::vector<LabelWeights*>::iterator label_weights =
          FindOrResizeFeatureKeyVector(feature_key_vector, i + 48 * current_64bitword_key);
        for (int k = 0; k < labels.size(); ++k) {
          AddValue(label_weights, labels[k], values[k]);
        }
      }
    }
  }
}

void SparseLabeledParameterVector::Add(const vector<uint64_t> &keys,
                                       const vector<int> &labels,
                                       const vector<double> &values) {
  for (int i = 0; i < keys.size(); ++i) {
    Add(keys[i], labels[i], values[i]);
  }
}

void SparseLabeledParameterVector::Add(const SparseLabeledParameterVector &parameters) {
  for (LabeledParameterMap::const_iterator iterator =
       parameters.map_values_.begin();
       iterator != parameters.map_values_.end();
       ++iterator) {
    uint64_t key = iterator->first;
    LabelWeights *label_weights = iterator->second;
    int label;
    double value;
    for (int k = 0; k < label_weights->Size(); ++k) {
      label_weights->GetLabelWeightByPosition(k, &label, &value);
      value *= parameters.scale_factor_;
      CHECK_EQ(value, value);
      //LOG(INFO) << value;
      Add(key, label, value);
    }
  }
}

void SparseLabeledParameterVector::GetValues(LabeledParameterMap::const_iterator iterator,
                                             const vector<int> &labels,
                                             vector<double> *values) const {
  values->resize(labels.size());
  LabelWeights *label_weights = iterator->second;
  for (int i = 0; i < labels.size(); ++i) {
    (*values)[i] = label_weights->GetWeight(labels[i]) * scale_factor_;
  }
}

void SparseLabeledParameterVector::GetValues(std::vector<LabelWeights *>::const_iterator iterator,
                                             const vector<int> &labels,
                                             vector<double> *values) const {
  values->resize(labels.size());
  LabelWeights *label_weights = *iterator;
  for (int i = 0; i < labels.size(); ++i) {
    (*values)[i] = label_weights->GetWeight(labels[i]) * scale_factor_;
  }
}

void SparseLabeledParameterVector::GetValuesAndUpdate(std::vector<LabelWeights *>::const_iterator iterator,
                                                      const vector<int> &labels,
                                                      vector<double> *values) const {
  values->resize(labels.size());
  LabelWeights *label_weights = *iterator;
  for (int i = 0; i < labels.size(); ++i) {
    (*values)[i] += label_weights->GetWeight(labels[i]) * scale_factor_;
  }
}

double SparseLabeledParameterVector::GetValue(
  LabeledParameterMap::const_iterator iterator,
  int label) const {
  LabelWeights *label_weights = iterator->second;
  return label_weights->GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  LabeledParameterMap::iterator iterator,
  int label) const {
  LabelWeights *label_weights = iterator->second;
  return label_weights->GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  std::vector<LabelWeights *>::const_iterator iterator,
  int label) const {
  LabelWeights *label_weights = *iterator;
  return label_weights->GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  std::vector<LabelWeights *>::iterator iterator,
  int label) const {
  LabelWeights *label_weights = *iterator;
  return label_weights->GetWeight(label) * scale_factor_;
}

void SparseLabeledParameterVector::SetValue(
  LabeledParameterMap::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights *label_weights = iterator->second;
  label_weights->SetWeight(label, value / scale_factor_);
#else
  LabelWeights *label_weights = iterator->second;
  double previous_value = label_weights->SetWeightAndNormalize(label,
                                                               value,
                                                               scale_factor_);
  squared_norm_ += value * value - previous_value * previous_value;
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights->Size() > kNumMaxSparseLabels &&
      label_weights->IsSparse()) {
    DenseLabelWeights *dense_label_weights =
      new DenseLabelWeights(label_weights);
    delete label_weights;
    iterator->second = dense_label_weights;
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

void SparseLabeledParameterVector::SetValue(
  std::vector<LabelWeights *>::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights *label_weights = *iterator;
  label_weights->SetWeight(label, value / scale_factor_);
#else
  LabelWeights *label_weights = *iterator;
  double previous_value = label_weights->SetWeightAndNormalize(label,
                                                               value,
                                                               scale_factor_);
  squared_norm_ += value * value - previous_value * previous_value;
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights->Size() > kNumMaxSparseLabels &&
      label_weights->IsSparse()) {
    DenseLabelWeights *dense_label_weights =
      new DenseLabelWeights(label_weights);
    delete label_weights;
    *iterator = dense_label_weights;
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

void SparseLabeledParameterVector::AddValue(
  LabeledParameterMap::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  value += current_value;
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights *label_weights = iterator->second;
  if (!label_weights)
    label_weights = new SparseLabelWeights;
  label_weights->SetWeight(label, value / scale_factor_);
#else
  LabelWeights *label_weights = iterator->second;
  if (!label_weights)
    label_weights = new SparseLabelWeights;
  double previous_value = label_weights->AddWeightAndNormalize(label,
                                                               value,
                                                               scale_factor_);
  //step1
  //squared_norm_ += (value + previous_value) * (value + previous_value)
  // - previous_value * previous_value;
  //step2
  //squared_norm_ += value * value + 2 * value * previous_value
  //+ ( previous_value * previous_value - previous_value * previous_value );
  squared_norm_ += value * value + 2 * value * previous_value;    //step3
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights->Size() > kNumMaxSparseLabels &&
      label_weights->IsSparse()) {
    DenseLabelWeights *dense_label_weights =
      new DenseLabelWeights(label_weights);
    delete label_weights;
    iterator->second = dense_label_weights;
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

void  SparseLabeledParameterVector::AddValue(
  std::vector<LabelWeights *>::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  value += current_value;
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights *label_weights = *iterator;
  if (!label_weights)
    label_weights = new SparseLabelWeights;
  label_weights->SetWeight(label, value / scale_factor_);
#else
  LabelWeights *label_weights = *iterator;
  if (!label_weights)
    label_weights = new SparseLabelWeights;
  double previous_value = label_weights->AddWeightAndNormalize(label,
                                                               value,
                                                               scale_factor_);
  //step1
  //squared_norm_ += (value + previous_value) * (value + previous_value)
  // - previous_value * previous_value;
  //step2
  //squared_norm_ += value * value + 2 * value * previous_value
  //+ ( previous_value * previous_value - previous_value * previous_value );
  squared_norm_ += value * value + 2 * value * previous_value;    //step3
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights->Size() > kNumMaxSparseLabels &&
      label_weights->IsSparse()) {
    DenseLabelWeights *dense_label_weights =
      new DenseLabelWeights(label_weights);
    delete label_weights;
    *iterator = dense_label_weights;
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

LabeledParameterMap::iterator
SparseLabeledParameterVector::FindOrInsert(uint64_t key) {
  LabeledParameterMap::iterator iterator = map_values_.find(key);
  if (iterator != map_values_.end() || growth_stopped()) return iterator;
  LabelWeights *label_weights = new SparseLabelWeights;
  pair<LabeledParameterMap::iterator, bool> result =
    map_values_.insert(pair<uint64_t, LabelWeights*>(key, label_weights));
  CHECK(result.second);
  return result.first;
}

std::vector<std::vector<LabelWeights*>>::const_iterator
SparseLabeledParameterVector::FindFamilyFeatureVector(uint64_t key) const {
  uint64_t feature_family_index = GetFamilyFeatureKey(key);
  if (GetMatrixFamilyFeatureVectorSize() <= feature_family_index) {
    return matrix_values_.end();
  }
  return matrix_values_.begin() + feature_family_index;
}

std::vector<std::vector<LabelWeights*>>::iterator
SparseLabeledParameterVector::FindOrResizeFamilyFeatureVector(uint64_t key) {
  uint64_t feature_family_index = GetFamilyFeatureKey(key);
  if (GetMatrixFamilyFeatureVectorSize() <= feature_family_index) {
    matrix_values_.resize(feature_family_index + 1);
  }
  return matrix_values_.begin() + feature_family_index;
}

std::vector<LabelWeights*>::const_iterator
SparseLabeledParameterVector::FindFeatureKeyVector(
  std::vector<std::vector<LabelWeights*>>::const_iterator feature_key_vector,
  uint64_t index) const {
  unsigned int matrix_size_y = GetMatrixFeatureKeyVectorSize(*feature_key_vector);
  if (matrix_size_y <= index) {
    return (*feature_key_vector).end();
  }
  return (*feature_key_vector).begin() + index;
}

std::vector<LabelWeights*>::iterator
SparseLabeledParameterVector::FindOrResizeFeatureKeyVector(
  std::vector<std::vector<LabelWeights*>>::iterator ptr_to_feature_key_vector,
  uint64_t index) {
  unsigned int matrix_size_y = GetMatrixFeatureKeyVectorSize(*ptr_to_feature_key_vector);
  if (matrix_size_y <= index) {
    (*ptr_to_feature_key_vector).resize(index + 1);
    //need to create a LabelWeights object for each new entry
    for (std::vector<LabelWeights*>::iterator iterator =
         (*ptr_to_feature_key_vector).begin() + matrix_size_y;
         iterator != (*ptr_to_feature_key_vector).end();
         iterator++) {
      (*iterator) = new SparseLabelWeights();
    }
  }
  return (*ptr_to_feature_key_vector).begin() + index;
}

void SparseLabeledParameterVector::ResizeFeatureKeyVector(
  std::vector<std::vector<LabelWeights*>>::iterator
  ptr_to_feature_key_vector, uint64_t size) {
  (*ptr_to_feature_key_vector).resize(size);
  //need to create a LabelWeights object for each entry
  for (std::vector<LabelWeights*>::iterator iterator =
       (*ptr_to_feature_key_vector).begin();
       iterator != (*ptr_to_feature_key_vector).end();
       iterator++) {
    (*iterator) = new SparseLabelWeights();
  }
}

// If the scale factor is too small, renormalize the entire parameter map.
void SparseLabeledParameterVector::RenormalizeIfNecessary() {
  if (scale_factor_ > -kLabeledScaleFactorThreshold &&
      scale_factor_ < kLabeledScaleFactorThreshold) {
    Renormalize();
  }
}

// Renormalize the entire parameter map (an expensive operation).
void SparseLabeledParameterVector::Renormalize() {
  LOG(INFO) << "Renormalizing the parameter map...";
  for (LabeledParameterMap::iterator iterator = map_values_.begin();
  iterator != map_values_.end();
    ++iterator) {
    LabelWeights *label_weights = iterator->second;
    int label;
    double value;
    for (int k = 0; k < label_weights->Size(); ++k) {
      label_weights->GetLabelWeightByPosition(k, &label, &value);
      label_weights->SetWeightByPosition(k, value * scale_factor_);
    }
  }
  scale_factor_ = 1.0;



  LOG(INFO) << 
    " ALERT: Renormalization is not yet implemented for the matrix-based map.";
}
