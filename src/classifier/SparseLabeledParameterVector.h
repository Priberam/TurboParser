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

#ifndef SPARSELABELEDPARAMETERVECTOR_H_
#define SPARSELABELEDPARAMETERVECTOR_H_

//#define USE_CUSTOMIZED_HASH_TABLE

#ifdef USE_CUSTOMIZED_HASH_TABLE
#include "HashTable.h"
#else
#ifdef _WIN32
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
#endif
#include "SerializationUtils.h"

using namespace std;

// Threshold for renormalizing the parameter vector.
const double kLabeledScaleFactorThreshold = 1e-9;
// After more than kNumMaxSparseLabels labels, use a dense representation.
const int kNumMaxSparseLabels = 5;

// This class contains the weights for every label conjoined with a single
// feature. This is a pure virtual class, so that we can derive from it a
// class for sparse label sets, and another for dense label sets.
// When training, features that get conjoined with more than kNumMaxSparseLabels
// use the dense variant, while the others use the sparse variant.
class LabelWeights {
public:
  LabelWeights() {};
  virtual ~LabelWeights() {};

  // True if set of labels is sparse.
  virtual bool IsSparse() const = 0;

  // Number of allocated labels.
  virtual int Size() const = 0;

  // Get/set/add weight to a labeled feature.
  virtual double GetWeight(int label) const = 0;
  virtual void SetWeight(int label, double weight) = 0;
  virtual void AddWeight(int label, double weight) = 0;
  virtual double SetWeightAndNormalize(int label, double value,
                                       double scaling_factor) = 0;
  virtual double AddWeightAndNormalize(int label, double value,
                                       double scaling_factor) = 0;

  // Get/set weight querying by position rather than the label.
  virtual void GetLabelWeightByPosition(int position, int *label,
                                        double *weight) const = 0;
  virtual void SetWeightByPosition(int position, double weight) = 0;
};

// Sparse implementation of LabelWeights.
class SparseLabelWeights : public LabelWeights {
public:
  SparseLabelWeights() {};
  virtual ~SparseLabelWeights() {};

  bool IsSparse() const;
  int Size() const;

  double GetWeight(int label) const;
  void SetWeight(int label, double weight);
  void AddWeight(int label, double weight);
  // Sets new weight, normalize it and returns previous value of
  // label_weights_[k].second, with k such that label == label_weights_[k].first.
  double SetWeightAndNormalize(int label, double value, double scale_factor);
  // Add weight value to current weight, normalize it and
  // returns previous value of label_weights_[k].second,
  // with k such that label == label_weights_[k].first.
  double AddWeightAndNormalize(int label, double value, double scale_factor);

  void GetLabelWeightByPosition(int position, int *label,
                                double *weight) const;

  void SetWeightByPosition(int position, double weight);

protected:
  std::vector<std::pair<int, double> > label_weights_;
};

// Dense implementation of LabelWeights.
class DenseLabelWeights : public LabelWeights {
public:
  DenseLabelWeights() {};
  DenseLabelWeights(LabelWeights *label_weights);
  virtual ~DenseLabelWeights() {};

  bool IsSparse() const;
  int Size() const;

  double GetWeight(int label) const;
  void SetWeight(int label, double weight);
  void AddWeight(int label, double weight);
  // Sets new weight, normalize it and returns previous value of weights_[label].
  double SetWeightAndNormalize(int label, double value, double scaling_factor);
  // Add weight value to current weight,
  // normalize it and returns previous value of weights_[label].
  double AddWeightAndNormalize(int label, double value, double scaling_factor);

  void GetLabelWeightByPosition(int position, int *label,
                                double *weight) const;

  void SetWeightByPosition(int position, double weight);

protected:
  vector<double> weights_;
};

// A labeled parameter map maps from feature keys ("labeled" features) to
// LabelWeights, which contain the weights of several labels conjoined with
// that feature.
#ifdef USE_CUSTOMIZED_HASH_TABLE
typedef HashTable<uint64_t, LabelWeights*> LabeledParameterMap;
#else
typedef std::tr1::unordered_map <uint64_t, LabelWeights*> LabeledParameterMap;
#endif

// This class implements a sparse parameter vector, which contains weights for
// the labels conjoined with each feature key. For fast lookup, this is
// implemented using an hash table.
// We represent a weight vector as a triple
// (map_values_, scale_factor_ , squared_norm_), where map_values_ contains the
// weights up to a scale, scale_factor_ is a factor such that
// weights[k] = scale_factor_ * map_values_[k], and the squared norm is cached.
// This way we can scale the weight vector in constant time (this operation is
// necessary in some training algorithms such as SGD), and manipulating a few
// elements is still fast. Plus, we can obtain the norm in constant time.
class SparseLabeledParameterVector {
public:
  SparseLabeledParameterVector() { growth_stopped_ = false; }
  virtual ~SparseLabeledParameterVector() { Clear(); }

  // Lock/unlock the parameter vector. If the vector is locked, no new features
  // can be inserted.
  void StopGrowth();
  void AllowGrowth();
  bool growth_stopped() const;

  // Clear the parameter vector.
  void Clear();

  // Save/load the parameters to/from a file.
  void Save(FILE *fs) const;

  void Load(FILE *fs);

  // Initialize to all-zeros.
  void Initialize();

  // Get the number of instantiated features.
     // This is the number of parameters up to different labels.
  uint64_t Size() const;

  // True if this feature key is already instantiated.
  bool Exists(uint64_t key) const;


  // Get the weights for the specified labels. Returns false if no key was
  // found, in which case weights becomes empty.
  bool Get(uint64_t key, const vector<int> &labels,
           vector<double> *weights) const;

  // Get squared norm of the parameter vector.
  double GetSquaredNorm() const;

  // Scale the weight vector by a factor scale_factor.
  // w_k' = w_k * c_k
  void Scale(double scale_factor);

  // Set the weight of this feature key conjoined with this label to "value".
  // Return false if the feature is not instantiated and cannot be inserted.
  // w'[id] = val
  bool Set(uint64_t key, int label, double value);

  // Increment the weight of this feature key conjoined with this label by an
  // amount of "value".
  // Return false if the feature is not instantiated and cannot be inserted.
  // w'[id] = w[id] + val
  bool Add(uint64_t key, int label, double value);

  // Increment the weights of this feature key conjoined with these labels by an
  // amount of "value".
  // Return false if the feature is not instantiated and cannot be inserted.
  // w'[id] = w[id] + val
  void Add(uint64_t key, const vector<int> &labels,
           const vector<double> &values);

  // Increment the weights of these feature keys paired with these labels by an
  // amount of "value".
  // NOTE: Silently bypasses the ones that could not be inserted, if any.
  // w'[id] = w[id] + val
  void Add(const vector<uint64_t> &keys, const vector<int> &labels,
           const vector<double> &values);

  // Adds two parameter vectors. This has the effect of incrementing the weights
  // of several features.
  // NOTE: Silently bypasses the ones that could not be inserted, if any.
  void Add(const SparseLabeledParameterVector &parameters);

protected:

  uint64_t SparseLabeledParameterVector::GetMapSize() const;
  uint64_t SparseLabeledParameterVector::GetMatrixSize() const;
  uint64_t SparseLabeledParameterVector::GetMatrixFamilyFeatureVectorSize() const;
  uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(int x) const;
  uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(const vector<LabelWeights*>&  feature_key_vector) const;
  uint64_t SparseLabeledParameterVector::GetMatrixFeatureKeyVectorSize(std::vector<std::vector<LabelWeights*>>::const_iterator ptr_to_feature_key_vector) const;
  uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(int x, int y) const;
  uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(LabelWeights* label_weights_vector) const;
  uint64_t SparseLabeledParameterVector::GetMatrixLabelWeightsVectorSize(std::vector<LabelWeights*>::const_iterator ptr_to_label_weights_vector) const;

  bool isHashMapKey(uint64_t key) const;
  bool isMatrixMapKey(uint64_t key) const;
  bool isMulti64bitKey(uint64_t key) const;
  uint64_t GetFamilyFeatureKey(uint64_t key) const;
  uint64_t GetBitKeys(uint64_t key) const;
  uint64_t Get3_5BitKeys(uint64_t key) const;
  uint64_t Get6_8BitKeys(uint64_t key) const;

  // Get the weights for the specified labels.
  void GetValues(LabeledParameterMap::const_iterator iterator,
                 const vector<int> &labels,
                 vector<double> *values) const;

  void GetValues(std::vector<LabelWeights *>::const_iterator iterator,
                 const vector<int> &labels,
                 vector<double> *values) const;
  void GetValuesAndUpdate(std::vector<LabelWeights *>::const_iterator iterator,
                          const vector<int> &labels,
                          vector<double> *values) const;

  // Get the weight for the specified label.
  // Two versions of this function: one using a const_iterator,
  // another using an iterator.
  double GetValue(LabeledParameterMap::const_iterator iterator,
                  int label) const;
  double GetValue(LabeledParameterMap::iterator iterator,
                  int label) const;
  double GetValue(std::vector<LabelWeights *>::const_iterator iterator,
                  int label) const;
  double GetValue(std::vector<LabelWeights *>::iterator iterator,
                  int label) const;

  // Set the weight for the specified label.
  void SetValue(LabeledParameterMap::iterator iterator,
                int label,
                double value);

  // Set the weight for the specified label.
  void SetValue(std::vector<LabelWeights *>::iterator iterator,
                int label,
                double value);


  // Add weight for the specified label.
  void AddValue(LabeledParameterMap::iterator iterator,
                int label,
                double value);
  // Add weight for the specified label.
  void AddValue(std::vector<LabelWeights *>::iterator iterator,
                int label,
                double value);

  // Find a key, or insert it in case it does not exist.
  LabeledParameterMap::iterator FindOrInsert(uint64_t key);

  std::vector<std::vector<LabelWeights*>>::const_iterator FindFamilyFeatureVector(uint64_t key) const;
  std::vector<std::vector<LabelWeights*>>::iterator FindOrResizeFamilyFeatureVector(uint64_t key);

  std::vector<LabelWeights*>::const_iterator FindFeatureKeyVector(
    std::vector<std::vector<LabelWeights*>>::const_iterator feature_key_vector,
    uint64_t index) const;
  std::vector<LabelWeights*>::iterator FindOrResizeFeatureKeyVector(
    std::vector<std::vector<LabelWeights*>>::iterator feature_key_vector,
    uint64_t index);
  void SparseLabeledParameterVector::ResizeFeatureKeyVector(
    std::vector<std::vector<LabelWeights*>>::iterator feature_key_vector, uint64_t size);


  // If the scale factor is too small, renormalize the entire parameter map.
  void RenormalizeIfNecessary();

  // Renormalize the entire parameter map (an expensive operation).
  void Renormalize();

protected:
  LabeledParameterMap map_values_; // Weight values, up to a scale. Hash-table-based storage.
  vector<vector<LabelWeights*>> matrix_values_; // Weight values, up to a scale. Vector-based storage.
  double scale_factor_; // The scale factor, such that w = values * scale.
  double squared_norm_; // The squared norm of the parameter vector.
  bool growth_stopped_; // True if parameters are locked.
};

#endif /*SPARSELABELEDPARAMETERVECTOR_H_*/
