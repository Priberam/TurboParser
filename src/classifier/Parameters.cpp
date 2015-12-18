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

#include "Parameters.h"
#include <iostream>
#include <math.h>

void Parameters::Save(FILE *fs) {
  int total_occurences_of_combinations = 0;
  int counter1 = 0;
  int counter1_total_occurences_of_combinations = 0;
  int counter10 = 0;
  int counter10_total_occurences_of_combinations = 0;
  int counter100 = 0;
  int counter100_total_occurences_of_combinations = 0;
  int counter1000 = 0;
  int counter1000_total_occurences_of_combinations = 0;
  for (std::unordered_map<bitmap_64bits_t, occurrences_counter_t>::iterator iterator =
       popular_keybitmap_occurrences_counter.begin();
       iterator != popular_keybitmap_occurrences_counter.end();
       iterator++) {
    total_occurences_of_combinations += (*iterator).second;
    if ((*iterator).second == 1) {
      counter1++;
      counter1_total_occurences_of_combinations += (*iterator).second;
    }
    if ((*iterator).second > 10) {
      counter10++;
      counter10_total_occurences_of_combinations += (*iterator).second;
    }
    if ((*iterator).second > 100) {
      counter100++;
      counter100_total_occurences_of_combinations += (*iterator).second;
      //LOG(INFO) << "Found a combination with over 100 occurrences";
    }
    if ((*iterator).second > 1000) {
      counter1000++;
      counter1000_total_occurences_of_combinations += (*iterator).second;
      LOG(INFO) << "Found a combination with over 1000 occurrences "
        << "(" << (*iterator).second << ")";
    }
  }
  LOG(INFO) << std::endl;
  LOG(INFO) << "Found " << counter1
    << " combinations with only 1 occurrence out of "
    << popular_keybitmap_occurrences_counter.size()
    << " combinations." << std::endl << "\t\t\t\t"
    << "which account for "
    << 100.0 * counter1_total_occurences_of_combinations / total_occurences_of_combinations
    << "% of total occurrences of every processed combination";

  LOG(INFO) << "Found " << counter10
    << " combinations with over 10 occurrences out of "
    << popular_keybitmap_occurrences_counter.size()
    << " combinations." << std::endl << "\t\t\t\t"
    << "which account for "
    << 100.0 * counter10_total_occurences_of_combinations / total_occurences_of_combinations
    << "% of total occurrences of every processed combination";

  LOG(INFO) << "Found " << counter100
    << " combinations with over 100 occurrences out of "
    << popular_keybitmap_occurrences_counter.size()
    << " combinations." << std::endl << "\t\t\t\t"
    << "which account for "
    << 100.0 * counter100_total_occurences_of_combinations / total_occurences_of_combinations
    << "% of total occurrences of every processed combination";

  LOG(INFO) << "Found " << counter1000
    << " combinations with over 1000 occurrences out of "
    << popular_keybitmap_occurrences_counter.size()
    << " combinations." << std::endl << "\t\t\t\t"
    << "which account for "
    << 100.0 * counter1000_total_occurences_of_combinations / total_occurences_of_combinations
    << "% of total occurrences of every processed combination";

  LOG(INFO) << "'Add unigram features' calls with popular keys: " << add_unigram_features_calls_w_special_key;
  LOG(INFO) << "Total 'Add unigram features' calls: " << add_unigram_features_calls;

  //Save feature key combinations -> add them to label_weights_
  for (std::unordered_map<bitmap_64bits_t, occurrences_counter_t>::const_iterator iterator =
       popular_keybitmap_occurrences_counter.begin();
       iterator != popular_keybitmap_occurrences_counter.end();
       iterator++) {
    if ((*iterator).second >= 1) { //save them all
      LabelWeights combined_scores;
      std::bitset<NUM_SPECIAL_KEYS> bit_set = std::bitset<NUM_SPECIAL_KEYS>((*iterator).first);
      for (int i = 0; i < NUM_SPECIAL_KEYS; i++) {
        if (bit_set.test(i)) {
          const LabelWeights* temp = labeled_weights_.GetLabelWeights(popular_keys[i].first);
          temp->CopyToExternalLabelWeights(&combined_scores);
        }
      }
      uint64_t key = popular_keybitmap_generated_key.find((*iterator).first)->second;
      labeled_weights_.Set(key, combined_scores);
    }
  }

  weights_.Save(fs);
  labeled_weights_.Save(fs);

  bool success;
  //success = WriteInteger(fs, popular_keys.size());
  //CHECK(success);
  //for (std::array<key_64bits_t, NUM_SPECIAL_KEYS>::const_iterator
  //     iterator = popular_keys.begin();
  //     iterator != popular_keys.end();
  //     ++iterator) {
  //  success = WriteUINT64(fs, *iterator);
  //  CHECK(success);
  //}
  success = WriteUINT64(fs, or_key);
  CHECK(success);

  //success = WriteInteger(fs, popular_keys_index_in_bitset.size());
  success = WriteInteger(fs, popular_keys.size());
  CHECK(success);
  //for (std::unordered_map< key_64bits_t, index_t>::const_iterator
  //     iterator = popular_keys_index_in_bitset.begin();
  //     iterator != popular_keys_index_in_bitset.end();
  //     ++iterator) {
  for (int i = 0; i < popular_keys.size(); i++) {
    success = WriteUINT64(fs, popular_keys[i].first);
    CHECK(success);
    success = WriteUINT64(fs, popular_keys[i].second);
    CHECK(success);
  }

  success = WriteInteger(fs, popular_keybitmap_generated_key.size());
  CHECK(success);
  for (std::unordered_map<bitmap_64bits_t, key_64bits_t>::const_iterator
       iterator = popular_keybitmap_generated_key.begin();
       iterator != popular_keybitmap_generated_key.end();
       ++iterator) {
    success = WriteUINT64(fs, iterator->first);
    CHECK(success);
    success = WriteUINT64(fs, iterator->second);
    CHECK(success);
  }
}

void Parameters::Load(FILE *fs) {
  weights_.Load(fs);
  labeled_weights_.Load(fs);

  LOG(INFO) << "Squared norm of the weight vector = " << GetSquaredNorm();
  LOG(INFO) << "Number of features = " << Size();

  bool success;
  //int num_popular_keys;

  //success = ReadInteger(fs, &num_popular_keys);
  //CHECK(success);
  //for (int i = 0; i < num_popular_keys; ++i) {
  //  uint64_t key;
  //  success = ReadUINT64(fs, &key);
  //  CHECK(success);
  //  popular_keys[i] = key;
  //}
  success = ReadUINT64(fs, &or_key);
  CHECK(success);

  int num_popular_keys_index_in_bitset;

  success = ReadInteger(fs, &num_popular_keys_index_in_bitset);
  CHECK(success);
  for (int i = 0; i < num_popular_keys_index_in_bitset; ++i) {
    uint64_t key;
    uint64_t index;
    success = ReadUINT64(fs, &key);
    CHECK(success);
    success = ReadUINT64(fs, &index);
    CHECK(success);
    // popular_keys_index_in_bitset.insert({ key, index });
    popular_keys[i] = { key, index };
  }

  int num_popular_keybitmap_generated_key;

  success = ReadInteger(fs, &num_popular_keybitmap_generated_key);
  CHECK(success);
  for (int i = 0; i < num_popular_keybitmap_generated_key; ++i) {
    uint64_t bitset_value;
    uint64_t key;
    success = ReadUINT64(fs, &bitset_value);
    CHECK(success);
    success = ReadUINT64(fs, &key);
    CHECK(success);
    popular_keybitmap_generated_key.insert({ bitset_value, key });
  }
}
