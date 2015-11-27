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

#include "EntityPipe.h"
#include <iostream>
#include <sstream>
#include <vector>
#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

void EntityPipe::PreprocessData() {
  delete token_dictionary_;
  CreateTokenDictionary();
  static_cast<SequenceDictionary*>(dictionary_)->
    SetTokenDictionary(token_dictionary_);
  // To get the right reader (instead of the default sequence reader).
  static_cast<EntityTokenDictionary*>(token_dictionary_)->
    Initialize(GetEntityReader());
  static_cast<SequenceDictionary*>(dictionary_)->
    CreateTagDictionary(GetSequenceReader());
}


void EntityPipe::ComputeScores(Instance *instance,
                                 Parts *parts,
                                 Features *features,
                                 vector<double> *scores) {
  EntityInstanceNumeric *sentence =
    static_cast<EntityInstanceNumeric*>(instance);
  SequenceParts *sequence_parts = static_cast<SequenceParts*>(parts);
  EntityFeatures *entity_features =
    static_cast<EntityFeatures*>(features);
  SequenceDictionary *sequence_dictionary = GetSequenceDictionary();
  scores->resize(parts->size());

  // Compute scores for the unigram parts.
  for (int i = 0; i < sentence->size(); ++i) {
    const vector<int> &index_unigram_parts =
      sequence_parts->FindUnigramParts(i);
    vector<int> allowed_tags(index_unigram_parts.size());
    for (int k = 0; k < index_unigram_parts.size(); ++k) {
      SequencePartUnigram *unigram =
        static_cast<SequencePartUnigram*>((*parts)[index_unigram_parts[k]]);
      allowed_tags[k] = unigram->tag();
    }


    const MultiBinaryFeatures &UnigramMultiFeatures = 
      entity_features->GetUnigramMultiFeatures(i);
    for (int j; j < UnigramMultiFeatures.size(); j++) {
      // Conjoin unigram features with the tag.
      const BinaryFeatures &unigram_features =
        entity_features->GetUnigramFeatures(UnigramMultiFeatures, j);

      vector<double> tag_scores;
      parameters_->ComputeLabelScores(unigram_features,
                                      allowed_tags,
                                      &tag_scores);
      for (int k = 0; k < index_unigram_parts.size(); ++k) {
        (*scores)[index_unigram_parts[k]] = tag_scores[k];
      }
    }
  }

  // Compute scores for the bigram parts.
  if (GetSequenceOptions()->markov_order() >= 1) {
    for (int i = 0; i < sentence->size() + 1; ++i) {
      const vector<int> &index_bigram_parts = sequence_parts->FindBigramParts(i);
      vector<int> bigram_tags(index_bigram_parts.size());
      for (int k = 0; k < index_bigram_parts.size(); ++k) {
        SequencePartBigram *bigram =
          static_cast<SequencePartBigram*>((*parts)[index_bigram_parts[k]]);
        bigram_tags[k] = sequence_dictionary->GetBigramLabel(bigram->tag_left(),
                                                             bigram->tag());
      }

      const MultiBinaryFeatures &BigramMultiFeatures =
        entity_features->GetBigramMultiFeatures(i);
      for (int j; j < BigramMultiFeatures.size(); j++) {
        // Conjoin bigram features with the pair of tags.
        const BinaryFeatures &bigram_features =
          entity_features->GetBigramFeatures(BigramMultiFeatures, i);

        vector<double> tag_scores;
        parameters_->ComputeLabelScores(bigram_features,
                                        bigram_tags,
                                        &tag_scores);
        for (int k = 0; k < index_bigram_parts.size(); ++k) {
          (*scores)[index_bigram_parts[k]] = tag_scores[k];
        }
      }
    }
  }

  // Compute scores for the trigram parts.
  if (GetSequenceOptions()->markov_order() >= 2) {
    for (int i = 1; i < sentence->size() + 1; ++i) {
      const vector<int> &index_trigram_parts = sequence_parts->FindTrigramParts(i);
      vector<int> trigram_tags(index_trigram_parts.size());
      for (int k = 0; k < index_trigram_parts.size(); ++k) {
        SequencePartTrigram *trigram =
          static_cast<SequencePartTrigram*>((*parts)[index_trigram_parts[k]]);
        trigram_tags[k] = sequence_dictionary->GetTrigramLabel(
          trigram->tag_left_left(),
          trigram->tag_left(),
          trigram->tag());
      }

      const MultiBinaryFeatures &TrigramMultiFeatures =
        entity_features->GetTrigramMultiFeatures(i);
      for (int j; j < TrigramMultiFeatures.size(); j++) {
        // Conjoin trigram features with the triple of tags.
        const BinaryFeatures &trigram_features =
          entity_features->GetTrigramFeatures(TrigramMultiFeatures, i);

        vector<double> tag_scores;
        parameters_->ComputeLabelScores(trigram_features,
                                        trigram_tags,
                                        &tag_scores);
        for (int k = 0; k < index_trigram_parts.size(); ++k) {
          (*scores)[index_trigram_parts[k]] = tag_scores[k];
        }
      }
    }
  }
}

void EntityPipe::MakeGradientStep(Parts *parts,
                                    Features *features,
                                    double eta,
                                    int iteration,
                                    const vector<double> &gold_output,
                                    const vector<double> &predicted_output) {
  EntityFeatures *entity_features =
    static_cast<EntityFeatures*>(features);
  SequenceDictionary *sequence_dictionary = GetSequenceDictionary();

  for (int r = 0; r < parts->size(); ++r) {
    //LOG(INFO) << predicted_output[r] << " " << gold_output[r];
    if (predicted_output[r] == gold_output[r]) continue;

    if ((*parts)[r]->type() == SEQUENCEPART_UNIGRAM) {
      SequencePartUnigram *unigram =
        static_cast<SequencePartUnigram*>((*parts)[r]);


      const MultiBinaryFeatures &UnigramMultiFeatures =
        entity_features->GetUnigramMultiFeatures(unigram->position());
      for (int j; j < UnigramMultiFeatures.size(); j++) {
        const BinaryFeatures &unigram_features =
          entity_features->GetUnigramFeatures(UnigramMultiFeatures, j);

        parameters_->MakeLabelGradientStep(unigram_features, eta, iteration,
                                           unigram->tag(),
                                           predicted_output[r] - gold_output[r]);
      }
    } else if ((*parts)[r]->type() == SEQUENCEPART_BIGRAM) {
      SequencePartBigram *bigram =
        static_cast<SequencePartBigram*>((*parts)[r]);


      const MultiBinaryFeatures &BigramMultiFeatures =
        entity_features->GetBigramMultiFeatures(bigram->position());
      for (int j; j < BigramMultiFeatures.size(); j++) {
        const BinaryFeatures &bigram_features =
          entity_features->GetBigramFeatures(BigramMultiFeatures, j);
        int bigram_tag = sequence_dictionary->GetBigramLabel(bigram->tag_left(),
                                                             bigram->tag());

        parameters_->MakeLabelGradientStep(bigram_features, eta, iteration,
                                           bigram_tag,
                                           predicted_output[r] - gold_output[r]);
      }
    } else if ((*parts)[r]->type() == SEQUENCEPART_TRIGRAM) {
      SequencePartTrigram *trigram =
        static_cast<SequencePartTrigram*>((*parts)[r]);


      const MultiBinaryFeatures &TrigramMultiFeatures =
        entity_features->GetTrigramMultiFeatures(trigram->position());
      for (int j; j < TrigramMultiFeatures.size(); j++) {
        const BinaryFeatures &trigram_features =
          entity_features->GetTrigramFeatures(TrigramMultiFeatures, j);
        int trigram_tag =
          sequence_dictionary->GetTrigramLabel(trigram->tag_left_left(),
                                               trigram->tag_left(),
                                               trigram->tag());

        parameters_->MakeLabelGradientStep(trigram_features, eta, iteration,
                                           trigram_tag,
                                           predicted_output[r] - gold_output[r]);
      }
    } else {
      CHECK(false);
    }
  }
}


void EntityPipe::MakeFeatureDifference(Parts *parts,
                                         Features *features,
                                         const vector<double> &gold_output,
                                         const vector<double> &predicted_output,
                                         FeatureVector *difference) {
  EntityFeatures *entity_features =
    static_cast<EntityFeatures*>(features);
  SequenceDictionary *sequence_dictionary = GetSequenceDictionary();

  for (int r = 0; r < parts->size(); ++r) {
    if (predicted_output[r] == gold_output[r]) continue;

    if ((*parts)[r]->type() == SEQUENCEPART_UNIGRAM) {
      SequencePartUnigram *unigram =
        static_cast<SequencePartUnigram*>((*parts)[r]);


      const MultiBinaryFeatures &UnigramMultiFeatures =
        entity_features->GetUnigramMultiFeatures(unigram->position());
      for (int j; j < UnigramMultiFeatures.size(); j++) {
        const BinaryFeatures &unigram_features =
          entity_features->GetUnigramFeatures(UnigramMultiFeatures, j);
        for (int f = 0; f < unigram_features.size(); ++f) {
          difference->mutable_labeled_weights()->Add(unigram_features[f],
                                                     unigram->tag(),
                                                     predicted_output[r] - gold_output[r]);
        }
      }
    } else if ((*parts)[r]->type() == SEQUENCEPART_BIGRAM) {
      SequencePartBigram *bigram =
        static_cast<SequencePartBigram*>((*parts)[r]);


      const MultiBinaryFeatures &BigramMultiFeatures =
        entity_features->GetBigramMultiFeatures(bigram->position());
      for (int j; j < BigramMultiFeatures.size(); j++) {
        const BinaryFeatures &bigram_features =
          entity_features->GetBigramFeatures(BigramMultiFeatures, j);
        int bigram_tag = sequence_dictionary->GetBigramLabel(bigram->tag_left(),
                                                             bigram->tag());
        for (int f = 0; f < bigram_features.size(); ++f) {
          difference->mutable_labeled_weights()->Add(bigram_features[f],
                                                     bigram_tag,
                                                     predicted_output[r] - gold_output[r]);
        }
      }
    } else if ((*parts)[r]->type() == SEQUENCEPART_TRIGRAM) {
      SequencePartTrigram *trigram =
        static_cast<SequencePartTrigram*>((*parts)[r]);


      const MultiBinaryFeatures &TrigramMultiFeatures =
        entity_features->GetTrigramMultiFeatures(trigram->position());
      for (int j; j < TrigramMultiFeatures.size(); j++) {
        const BinaryFeatures &trigram_features =
          entity_features->GetTrigramFeatures(TrigramMultiFeatures, j);
        int trigram_tag =
          sequence_dictionary->GetTrigramLabel(trigram->tag_left_left(),
                                               trigram->tag_left(),
                                               trigram->tag());
        for (int f = 0; f < trigram_features.size(); ++f) {
          difference->mutable_labeled_weights()->Add(trigram_features[f],
                                                     trigram_tag,
                                                     predicted_output[r] - gold_output[r]);
        }
      }
    } else {
      CHECK(false);
    }
  }
}

void EntityPipe::MakeSelectedFeatures(Instance *instance,
                                        Parts *parts,
                                        const vector<bool> &selected_parts,
                                        Features *features) {
  EntityInstanceNumeric *sentence =
    static_cast<EntityInstanceNumeric*>(instance);
  EntityFeatures *entity_features =
    static_cast<EntityFeatures*>(features);

  int sentence_length = sentence->size();

  entity_features->Initialize(instance, parts);

  // Build features for words only. They will later be conjoined with the tags.
  for (int i = 0; i < sentence_length; ++i) {
    entity_features->AddUnigramFeatures(sentence, i);
  }

  if (GetSequenceOptions()->markov_order() >= 1) {
    for (int i = 0; i < sentence_length + 1; ++i) {
      entity_features->AddBigramFeatures(sentence, i);
    }
  }

  if (GetSequenceOptions()->markov_order() >= 2) {
    for (int i = 1; i < sentence_length + 1; ++i) {
      entity_features->AddTrigramFeatures(sentence, i);
    }
  }
}