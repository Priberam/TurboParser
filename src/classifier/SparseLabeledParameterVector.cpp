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
#include <iostream>
#include <sstream>
#include <vector>
#include <glog/logging.h>
#include <bitset>

bool LabelWeights::IsSparse() const {
  if (dense_mode) {
    return false;
  } else {
    return true;
  }
}
int LabelWeights::SparseSize() const {
  return sparse_label_weights_.size();
}
int LabelWeights::DenseSize() const {
  return dense_label_weights_.size();
}
int LabelWeights::Size() const {
  if (IsSparse()) {
    return SparseSize();
  } else {
    return DenseSize();
  }
}

double LabelWeights::SparseGetWeight(int label) const {
  for (int k = 0; k < sparse_label_weights_.size(); ++k) {
    if (label == sparse_label_weights_[k].first) {
      return sparse_label_weights_[k].second;
    }
  }
  return 0.0;
}

double LabelWeights::DenseGetWeight(int label) const {
  if (label >= dense_label_weights_.size()) return 0.0;
  return dense_label_weights_[label];
}

double LabelWeights::GetWeight(int label) const {
  if (IsSparse()) {
    return SparseGetWeight(label);
  } else {
    return DenseGetWeight(label);
  }
}

void LabelWeights::SparseSetWeight(int label, double weight) {
  for (int k = 0; k < sparse_label_weights_.size(); ++k) {
    if (label == sparse_label_weights_[k].first) {
      sparse_label_weights_[k].second = weight;
      return;
    }
  }
  sparse_label_weights_.push_back(std::pair<int, double>(label, weight));
}
void LabelWeights::DenseSetWeight(int label, double weight) {
  CHECK_GE(label, 0);
  if (label >= dense_label_weights_.size()) {
    dense_label_weights_.resize(label + 1, 0.0);
  }
  dense_label_weights_[label] = weight;
}

void LabelWeights::SetWeight(int label, double weight) {
  if (IsSparse()) {
    SparseSetWeight(label, weight);
  } else {
    DenseSetWeight(label, weight);
  }
}

void LabelWeights::SparseAddWeight(int label, double weight) {
  for (int k = 0; k < sparse_label_weights_.size(); ++k) {
    if (label == sparse_label_weights_[k].first) {
      sparse_label_weights_[k].second += weight;
      return;
    }
  }
  sparse_label_weights_.push_back(std::pair<int, double>(label, weight));
}
void LabelWeights::DenseAddWeight(int label, double weight) {
  CHECK_GE(label, 0);
  if (label >= dense_label_weights_.size()) {
    dense_label_weights_.resize(label + 1, 0.0);
  }
  dense_label_weights_[label] += weight;
}

void LabelWeights::AddWeight(int label, double weight) {
  if (weight == 0.0)
    return;
  if (IsSparse()) {
    SparseAddWeight(label, weight);
  } else {
    DenseAddWeight(label, weight);
  }
}

double LabelWeights::SparseSetWeightAndNormalize(int label,
                                                 double value,
                                                 double scale_factor) {
  double weight = value / scale_factor;
  double previous_value;
  for (int k = 0; k < sparse_label_weights_.size(); ++k) {
    if (label == sparse_label_weights_[k].first) {
      previous_value = sparse_label_weights_[k].second * scale_factor;
      sparse_label_weights_[k].second = weight;
      return previous_value;
    }
  }
  sparse_label_weights_.push_back(std::pair<int, double>(label, weight));
  return 0.0;
}

double LabelWeights::DenseSetWeightAndNormalize(int label,
                                                double value,
                                                double scale_factor) {
  CHECK_GE(label, 0);
  if (label >= dense_label_weights_.size()) {
    dense_label_weights_.resize(label + 1, 0.0);
  }
  double weight = value / scale_factor;
  double previous_value = dense_label_weights_[label] * scale_factor;
  dense_label_weights_[label] = weight;
  return previous_value;
}

double LabelWeights::SetWeightAndNormalize(int label,
                                           double value,
                                           double scale_factor) {
  if (IsSparse()) {
    return SparseSetWeightAndNormalize(label,
                                       value,
                                       scale_factor);
  } else {
    return DenseSetWeightAndNormalize(label,
                                      value,
                                      scale_factor);
  }
}

double LabelWeights::SparseAddWeightAndNormalize(int label,
                                                 double value,
                                                 double scale_factor) {
  double weight = value / scale_factor;
  double previous_value;
  for (int k = 0; k < sparse_label_weights_.size(); ++k) {
    if (label == sparse_label_weights_[k].first) {
      previous_value = sparse_label_weights_[k].second * scale_factor;
      sparse_label_weights_[k].second += weight;
      return previous_value;
    }
  }
  sparse_label_weights_.push_back(std::pair<int, double>(label, weight));
  return 0.0;
}
double LabelWeights::DenseAddWeightAndNormalize(int label,
                                                double value,
                                                double scale_factor) {
  CHECK_GE(label, 0);
  if (label >= dense_label_weights_.size()) {
    dense_label_weights_.resize(label + 1, 0.0);
  }
  double weight = value / scale_factor;
  double previous_value = dense_label_weights_[label] * scale_factor;
  dense_label_weights_[label] += weight;
  return previous_value;
}
double LabelWeights::AddWeightAndNormalize(int label,
                                           double value,
                                           double scale_factor) {
  if (IsSparse()) {
    return SparseAddWeightAndNormalize(label,
                                       value,
                                       scale_factor);
  } else {
    return DenseAddWeightAndNormalize(label,
                                      value,
                                      scale_factor);
  }
}

void LabelWeights::SparseGetLabelWeightByPosition(int position, int *label,
                                                  double *weight) const {
  *label = sparse_label_weights_[position].first;
  *weight = sparse_label_weights_[position].second;
  CHECK_GE(*label, 0);
}
void LabelWeights::DenseGetLabelWeightByPosition(int position, int *label,
                                                 double *weight) const {
  CHECK_GE(position, 0);
  *label = position;
  *weight = dense_label_weights_[position];
}
void LabelWeights::GetLabelWeightByPosition(int position, int *label,
                                            double *weight) const {
  if (IsSparse()) {
    return SparseGetLabelWeightByPosition(position, label, weight);
  } else {
    return DenseGetLabelWeightByPosition(position, label, weight);
  }
}

void LabelWeights::SparseSetWeightByPosition(int position, double weight) {
  sparse_label_weights_[position].second = weight;
}
void LabelWeights::DenseSetWeightByPosition(int position, double weight) {
  dense_label_weights_[position] = weight;
}
void LabelWeights::SetWeightByPosition(int position, double weight) {
  if (IsSparse()) {
    SparseSetWeightByPosition(position, weight);
  } else {
    DenseSetWeightByPosition(position, weight);
  }
}

void LabelWeights::ChangeToDenseLabelWeights() {
  CHECK(IsSparse());
  int size = SparseSize();
  for (int k = 0; k < size; ++k) {
    int label;
    double weight;
    SparseGetLabelWeightByPosition(k, &label, &weight);
    CHECK_GE(label, 0);
    DenseSetWeight(label, weight);
  }
  SetDenseMode();
  ClearSparseData();
}

void LabelWeights::SetDenseMode() {
  dense_mode = true;
}

void LabelWeights::ClearSparseData() {
  sparse_label_weights_.clear();
}

void LabelWeights::UpdateExternalPartialScore(std::vector<double>* to,
                                              double multiplier) const {
  if (IsSparse()) {
    int size = sparse_label_weights_.size();
    auto find = std::max_element(begin(sparse_label_weights_),
                                 end(sparse_label_weights_),
                                 [](const std::pair<int, double>& pairA,
                                    const std::pair<int, double>& pairB) {
      return pairA.first < pairB.first;
    });

    int new_size = (*find).first + 1;
    if (to->size() < new_size)
      to->resize(new_size, 0.0);
    for (int i = 0; i < size; i++) {
      const std::pair<int, double> & pair = sparse_label_weights_[i];
      (*to)[pair.first] += pair.second * multiplier;
    }
  } else {
    int size = dense_label_weights_.size();
    if (to->size() < size)
      to->resize(size, 0.0);
    for (int i = 0; i < size; i++)
      (*to)[i] += dense_label_weights_[i] * multiplier;
  }
}

// Lock/unlock the parameter vector. If the vector is locked, no new features
// can be inserted.
void SparseLabeledParameterVector::StopGrowth() { growth_stopped_ = true; }
void SparseLabeledParameterVector::AllowGrowth() { growth_stopped_ = false; }
bool SparseLabeledParameterVector::growth_stopped() const { return growth_stopped_; }

// Clear the parameter vector.
void SparseLabeledParameterVector::Clear() {
  map_values_.clear();
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
    const LabelWeights &label_weights = iterator->second;
    int length = label_weights.Size();
    success = WriteInteger(fs, length);
    CHECK(success);
    int label;
    double value;
    for (int k = 0; k < length; ++k) {
      label_weights.GetLabelWeightByPosition(k, &label, &value);
      CHECK_GE(label, 0);
      success = WriteInteger(fs, label);
      CHECK(success);
      success = WriteDouble(fs, value);
      CHECK(success);
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

#define PRINT_STATISTICS
#ifdef PRINT_STATISTICS
  // Print some statistics:
  int num_sparse = 0;
  int num_total = 0;
  int num_labels_sparse = 0;
  for (LabeledParameterMap::iterator iterator = map_values_.begin();
  iterator != map_values_.end();
    ++iterator) {
    LabelWeights &label_weights = iterator->second;
    if (label_weights.IsSparse()) {
      ++num_sparse;
      int length = label_weights.Size();
      num_labels_sparse += length;
    }
    ++num_total;
  }

  LOG(INFO) << "Statistics for labeled parameter vector:";
  LOG(INFO) << "Features with sparse labels: " << num_sparse
    << " Total: " << num_total
    << " Sparse labels: " << num_labels_sparse;
#endif
}

void SparseLabeledParameterVector::Initialize() {
  Clear();
  scale_factor_ = 1.0;
  squared_norm_ = 0.0;
}

uint64_t SparseLabeledParameterVector::Size() const {
  return GetMapSize();
}
uint64_t SparseLabeledParameterVector::GetMapSize() const {
  return map_values_.size();
}

bool SparseLabeledParameterVector::Exists(uint64_t key) const {
  LabeledParameterMap::const_iterator iterator = map_values_.find(key);
  if (iterator == map_values_.end()) return false;
  return true;
}

bool SparseLabeledParameterVector::Get(uint64_t key, const vector<int> &labels,
                                       vector<double> *weights) const {
  LabeledParameterMap::const_iterator iterator = map_values_.find(key);
  if (iterator == map_values_.end()) {
    weights->clear();
    return false;
  }
  GetValues(iterator, labels, weights);
  return true;
}

const LabelWeights* SparseLabeledParameterVector::GetLabelWeights(uint64_t key) const {
  LabeledParameterMap::const_iterator iterator = map_values_.find(key);
  if (iterator == map_values_.end()) {
    return NULL;
  }
  return &iterator->second;
}

double SparseLabeledParameterVector::GetSquaredNorm() const {
  return squared_norm_;
}

double SparseLabeledParameterVector::GetScaleFactor() const {
  return scale_factor_;
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

  LabeledParameterMap::iterator iterator = FindOrInsert(key);
  if (iterator != map_values_.end()) {
    SetValue(iterator, label, value);
    return true;
  } else {
    return false;
  }
}

bool SparseLabeledParameterVector::Add(uint64_t key,
                                       int label,
                                       double value) {
  LabeledParameterMap::iterator iterator = FindOrInsert(key);
  if (iterator != map_values_.end()) {
    AddValue(iterator, label, value);
    return true;
  } else {
    return false;
  }
}

void SparseLabeledParameterVector::Add(uint64_t key,
                                       const vector<int> &labels,
                                       const vector<double> &values) {
  LabeledParameterMap::iterator iterator = FindOrInsert(key);
  for (int k = 0; k < labels.size(); ++k) {
    AddValue(iterator, labels[k], values[k]);
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
    const LabelWeights & label_weights = iterator->second;
    int label;
    double value;
    for (int k = 0; k < label_weights.Size(); ++k) {
      label_weights.GetLabelWeightByPosition(k, &label, &value);
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
  const LabelWeights &label_weights = iterator->second;
  for (int i = 0; i < labels.size(); ++i) {
    (*values)[i] = label_weights.GetWeight(labels[i]) * scale_factor_;
  }
}

void SparseLabeledParameterVector::GetValues(std::vector<LabelWeights>::const_iterator iterator,
                                             const vector<int> &labels,
                                             vector<double> *values) const {
  values->resize(labels.size());
  const LabelWeights & label_weights = *iterator;
  for (int i = 0; i < labels.size(); ++i) {
    (*values)[i] = label_weights.GetWeight(labels[i]) * scale_factor_;
  }
}

double SparseLabeledParameterVector::GetValue(
  LabeledParameterMap::const_iterator iterator,
  int label) const {
  const LabelWeights &label_weights = iterator->second;
  return label_weights.GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  LabeledParameterMap::iterator iterator,
  int label) const {
  LabelWeights &label_weights = iterator->second;
  return label_weights.GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  std::vector<LabelWeights>::const_iterator iterator,
  int label) const {
  const LabelWeights &label_weights = *iterator;
  return label_weights.GetWeight(label) * scale_factor_;
}

double SparseLabeledParameterVector::GetValue(
  std::vector<LabelWeights>::iterator iterator,
  int label) const {
  LabelWeights &label_weights = *iterator;
  return label_weights.GetWeight(label) * scale_factor_;
}

void SparseLabeledParameterVector::SetValue(
  LabeledParameterMap::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights &label_weights = iterator->second;
  label_weights.SetWeight(label, value / scale_factor_);
#else
  LabelWeights &label_weights = iterator->second;
  double previous_value = label_weights.SetWeightAndNormalize(label,
                                                              value,
                                                              scale_factor_);
  squared_norm_ += value * value - previous_value * previous_value;
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights.Size() > kNumMaxSparseLabels &&
      label_weights.IsSparse()) {
    label_weights.ChangeToDenseLabelWeights();
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

void SparseLabeledParameterVector::SetValue(
  std::vector<LabelWeights>::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights &label_weights = *iterator;
  label_weights.SetWeight(label, value / scale_factor_);
#else
  LabelWeights &label_weights = *iterator;
  double previous_value = label_weights.SetWeightAndNormalize(label,
                                                              value,
                                                              scale_factor_);
  squared_norm_ += value * value - previous_value * previous_value;
#endif

  // If the number of labels is growing large, make this into dense
  // label weights.
  if (label_weights.Size() > kNumMaxSparseLabels &&
      label_weights.IsSparse()) {
    label_weights.ChangeToDenseLabelWeights();
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
  LabelWeights &label_weights = iterator->second;
  label_weights.SetWeight(label, value / scale_factor_);
#else
  LabelWeights &label_weights = iterator->second;
  double previous_value = label_weights.AddWeightAndNormalize(label,
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
  if (label_weights.Size() > kNumMaxSparseLabels &&
      label_weights.IsSparse()) {
    label_weights.ChangeToDenseLabelWeights();
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

void  SparseLabeledParameterVector::AddValue(
  std::vector<LabelWeights>::iterator iterator,
  int label,
  double value) {
#if USE_N_OPTIMIZATIONS==0
  // TODO: Make this more efficient, avoiding two lookups in LabelWeights.
  double current_value = GetValue(iterator, label);
  value += current_value;
  squared_norm_ += value * value - current_value * current_value;
  LabelWeights &label_weights = *iterator;
  label_weights.SetWeight(label, value / scale_factor_);
#else
  LabelWeights &label_weights = *iterator;
  double previous_value = label_weights.AddWeightAndNormalize(label,
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
  if (label_weights.Size() > kNumMaxSparseLabels &&
      label_weights.IsSparse()) {
    label_weights.ChangeToDenseLabelWeights();
  }

  // This prevents numerical issues:
  if (squared_norm_ < 0.0) squared_norm_ = 0.0;
}

LabeledParameterMap::iterator
SparseLabeledParameterVector::FindOrInsert(uint64_t key) {
  LabeledParameterMap::iterator iterator = map_values_.find(key);
  if (iterator != map_values_.end() || growth_stopped()) return iterator;
  pair<LabeledParameterMap::iterator, bool> result =
    map_values_.emplace(key, LabelWeights());
  CHECK(result.second);
  return result.first;
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
    LabelWeights &label_weights = iterator->second;
    int label;
    double value;
    for (int k = 0; k < label_weights.Size(); ++k) {
      label_weights.GetLabelWeightByPosition(k, &label, &value);
      label_weights.SetWeightByPosition(k, value * scale_factor_);
    }
  }

  scale_factor_ = 1.0;
}
