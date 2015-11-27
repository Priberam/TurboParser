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
#include "EntityFeatures.h"
#include "SequencePart.h"
#include <bitset>

void EntityFeatures::AddUnigramFeatures(EntityInstanceNumeric *sentence,
                                        int position) {
  CHECK(!input_features_unigrams_[position]);

  MultiBinaryFeatures *features = new MultiBinaryFeatures;
  input_features_unigrams_[position] = features;

  int sentence_length = sentence->size();

  EntityInstanceNumeric *entity_sentence =
    static_cast<EntityInstanceNumeric*>(sentence);

  EntityOptions *options = static_cast<class EntityPipe*>(pipe_)->
    GetEntityOptions();

  // Array of form IDs.
  const vector<int>* word_ids = &entity_sentence->GetFormIds();

  // Array of POS IDs.
  const vector<int>* pos_ids = &entity_sentence->GetPosIds();

  // Words.
  uint16_t WID = (*word_ids)[position]; // Current word.
  // Word on the left.
  uint16_t pWID = (position > 0) ? (*word_ids)[position - 1] : TOKEN_START;
  // Word on the right.
  uint16_t nWID = (position < sentence_length - 1) ?
    (*word_ids)[position + 1] : TOKEN_STOP;
  // Word two positions on the left.
  uint16_t ppWID = (position > 1) ? (*word_ids)[position - 2] : TOKEN_START;
  // Word two positions on the right.
  uint16_t nnWID = (position < sentence_length - 2) ?
    (*word_ids)[position + 2] : TOKEN_STOP;

  // Gazetteer tags.
  std::vector<int> empty_GIDs;
  // Current gazetter tag.
  const std::vector<int> &GIDs = entity_sentence->GetGazetteerIds(position);
  // Gazetteer tag on the left.
  const std::vector<int> &pGIDs = (position > 0) ?
    entity_sentence->GetGazetteerIds(position - 1) : empty_GIDs;
  // Gazetteer tag on the right.
  const std::vector<int> &nGIDs = (position < sentence_length - 1) ?
    entity_sentence->GetGazetteerIds(position + 1) : empty_GIDs;
  // Gazetteer tag two positions on the left.
  const std::vector<int> &ppGIDs = (position > 1) ?
    entity_sentence->GetGazetteerIds(position - 2) : empty_GIDs;
  // Gazetteer tag two positions on the right.
  const std::vector<int> &nnGIDs = (position < sentence_length - 2) ?
    entity_sentence->GetGazetteerIds(position + 2) : empty_GIDs;

  // POS tags.
  uint8_t PID = (*pos_ids)[position]; // Current POS.
  // POS on the left.
  uint8_t pPID = (position > 0) ?
    (*pos_ids)[position - 1] : TOKEN_START;
  // POS on the right.
  uint8_t nPID = (position < sentence_length - 1) ?
    (*pos_ids)[position + 1] : TOKEN_STOP;
  // POS two positions on the left.
  uint8_t ppPID = (position > 1) ?
    (*pos_ids)[position - 2] : TOKEN_START;
  // POS two positions on the right.
  uint8_t nnPID = (position < sentence_length - 2) ?
    (*pos_ids)[position + 2] : TOKEN_STOP;

  // Word shapes.
  uint16_t SID = sentence->GetShapeId(position); // Current shape.
  // Shape on the left.
  uint16_t pSID = (position > 0) ?
    sentence->GetShapeId(position - 1) : TOKEN_START;
  // Shape on the right.
  uint16_t nSID = (position < sentence_length - 1) ?
    sentence->GetShapeId(position + 1) : TOKEN_STOP;
  // Shape two positions on the left.
  uint16_t ppSID = (position > 1) ?
    sentence->GetShapeId(position - 2) : TOKEN_START;
  // Shape two positions on the right.
  uint16_t nnSID = (position < sentence_length - 2) ?
    sentence->GetShapeId(position + 2) : TOKEN_STOP;

  // Prefixes/Suffixes.
  vector<uint16_t> AID(sentence->GetMaxPrefixLength(position), 0xffff);
  vector<uint16_t> ZID(sentence->GetMaxSuffixLength(position), 0xffff);
  for (int l = 0; l < AID.size(); ++l) {
    AID[l] = sentence->GetPrefixId(position, l + 1);
  }
  for (int l = 0; l < ZID.size(); ++l) {
    ZID[l] = sentence->GetSuffixId(position, l + 1);
  }

  // Several flags.
  uint8_t flag_all_digits = sentence->AllDigits(position) ? 0x1 : 0x0;
  uint8_t flag_all_digits_with_punctuation =
    sentence->AllDigitsWithPunctuation(position) ? 0x1 : 0x0;
  uint8_t flag_all_upper = sentence->AllUpper(position) ? 0x1 : 0x0;
  uint8_t flag_first_upper = position > 0 && sentence->FirstUpper(position) ?
    0x1 : 0x0;

  flag_all_digits = 0x0 | (flag_all_digits << 4);
  flag_all_digits_with_punctuation =
    0x1 | (flag_all_digits_with_punctuation << 4);
  flag_all_upper = 0x2 | (flag_all_upper << 4);
  flag_first_upper = 0x3 | (flag_first_upper << 4);

  uint64_t fkey;
  uint8_t flags = 0x0;

  // Maximum is 255 feature templates.
  CHECK_LT(EntityFeatureTemplateUnigram::COUNT, 256);

  // Bias feature.
  fkey = encoder_.CreateFKey_NONE(EntityFeatureTemplateUnigram::BIAS,
                                  flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::BIAS);

  // Lexical features.
  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::W,
                               WID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::W);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::pW,
                               pWID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::pW);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nW,
                               nWID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::nW);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::ppW,
                               ppWID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::ppW);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nnW,
                               nnWID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::nnW);

  // Gazetteer features.
  for (int k = 0; k < GIDs.size(); ++k) {
    uint16_t GID = GIDs[k];
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::G,
                                 GID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::G);
  }
  for (int k = 0; k < pGIDs.size(); ++k) {
    uint16_t pGID = pGIDs[k];
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::pG,
                                 pGID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::pG);
  }
  for (int k = 0; k < nGIDs.size(); ++k) {
    uint16_t nGID = nGIDs[k];
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nG,
                                 nGID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::nG);
  }
  for (int k = 0; k < ppGIDs.size(); ++k) {
    uint16_t ppGID = ppGIDs[k];
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::ppG,
                                 ppGID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::ppG);
  }
  for (int k = 0; k < nnGIDs.size(); ++k) {
    uint16_t nnGID = nnGIDs[k];
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nnG,
                                 nnGID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::nnG);
  }

  // POS features.
  fkey = encoder_.CreateFKey_P(EntityFeatureTemplateUnigram::P,
                               PID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::P);

  fkey = encoder_.CreateFKey_PP(EntityFeatureTemplateUnigram::PpP,
                                PID, pPID,
                                flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::PpP);

  fkey = encoder_.CreateFKey_PP(EntityFeatureTemplateUnigram::PnP,
                                PID, nPID,
                                flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::PnP);

  fkey = encoder_.CreateFKey_PPP(EntityFeatureTemplateUnigram::PpPppP,
                                 PID, pPID, ppPID,
                                 flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::PpPppP);

  fkey = encoder_.CreateFKey_PPP(EntityFeatureTemplateUnigram::PnPnnP,
                                 PID, nPID, nnPID,
                                 flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::PnPnnP);

  // Shape features.
  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::S,
                               SID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::S);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::pS,
                               pSID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::pS);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nS,
                               nSID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::nS);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::ppS,
                               ppSID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::ppS);

  fkey = encoder_.CreateFKey_W(EntityFeatureTemplateUnigram::nnS,
                               nnSID,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::nnS);

  // Prefix/Suffix features.
  for (int l = 0; l < AID.size(); ++l) {
    uint8_t flag_prefix_length = l;
    fkey = encoder_.CreateFKey_WP(EntityFeatureTemplateUnigram::A,
                                  AID[l], flag_prefix_length,
                                  flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::A);
  }
  for (int l = 0; l < ZID.size(); ++l) {
    uint8_t flag_suffix_length = l;
    fkey = encoder_.CreateFKey_WP(EntityFeatureTemplateUnigram::Z,
                                  ZID[l], flag_suffix_length,
                                  flags);
    AddFeature(fkey, features, EntityFeatureTemplateUnigram::Z);
  }

  // Several flags.
  fkey = encoder_.CreateFKey_P(EntityFeatureTemplateUnigram::FLAG,
                               flag_all_digits,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::FLAG);

  fkey = encoder_.CreateFKey_P(EntityFeatureTemplateUnigram::FLAG,
                               flag_all_digits_with_punctuation,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::FLAG);

  fkey = encoder_.CreateFKey_P(EntityFeatureTemplateUnigram::FLAG,
                               flag_all_upper,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::FLAG);

  fkey = encoder_.CreateFKey_P(EntityFeatureTemplateUnigram::FLAG,
                               flag_first_upper,
                               flags);
  AddFeature(fkey, features, EntityFeatureTemplateUnigram::FLAG);
}

void EntityFeatures::AddBigramFeatures(EntityInstanceNumeric *sentence,
                                       int position) {
  CHECK(!input_features_bigrams_[position]) << position
    << " " << sentence->size();

  MultiBinaryFeatures *features = new MultiBinaryFeatures;
  input_features_bigrams_[position] = features;

  uint64_t fkey;
  uint8_t flags = 0x0;
  flags |= EntityFeatureTemplateParts::BIGRAM;

  // Note: position ranges between 0 and N (inclusive), where N is the number
  // of words. If position = N, we need to be careful not to access invalid
  // memory in arrays.

  // Add other bigram features.
  int sentence_length = sentence->size();

  EntityInstanceNumeric *entity_sentence =
    static_cast<EntityInstanceNumeric*>(sentence);

  EntityOptions *options = static_cast<class EntityPipe*>(pipe_)->
    GetEntityOptions();
  std::bitset<32> *feature_set_bitmap;
  feature_set_bitmap = new bitset<32>(options->large_feature_set());

  // Array of form IDs.
  const vector<int>* word_ids = &entity_sentence->GetFormIds();

  // Array of POS IDs.
  const vector<int>* pos_ids = &entity_sentence->GetPosIds();

  // Words.
  uint16_t WID, pWID, nWID, ppWID, nnWID;
  if (feature_set_bitmap->test(0)) {
    WID = (position < sentence_length) ?
      (*word_ids)[position] : TOKEN_STOP; // Current word.
    // Word on the left.
    pWID = (position > 0) ?
      (*word_ids)[position - 1] : TOKEN_START;
    // Word on the right.
    nWID = (position < sentence_length - 1) ?
      (*word_ids)[position + 1] : TOKEN_STOP;
    // Word two positions on the left.
    ppWID = (position > 1) ?
      (*word_ids)[position - 2] : TOKEN_START;
    // Word two positions on the right.
    nnWID = (position < sentence_length - 2) ?
      (*word_ids)[position + 2] : TOKEN_STOP;
  }

  // POS tags.
  uint8_t PID, pPID, nPID, ppPID, nnPID;
  if (feature_set_bitmap->test(1)) {
    PID = (position < sentence_length) ?
      (*pos_ids)[position] : TOKEN_STOP; // Current POS.
    // POS on the left.
    pPID = (position > 0) ?
      (*pos_ids)[position - 1] : TOKEN_START;
    // POS on the right.
    nPID = (position < sentence_length - 1) ?
      (*pos_ids)[position + 1] : TOKEN_STOP;
    // POS two positions on the left.
    ppPID = (position > 1) ?
      (*pos_ids)[position - 2] : TOKEN_START;
    // POS two positions on the right.
    nnPID = (position < sentence_length - 2) ?
      (*pos_ids)[position + 2] : TOKEN_STOP;
  }

  // Word shapes.
  uint16_t SID, pSID, nSID, ppSID, nnSID;
  if (feature_set_bitmap->test(2)) {
    SID = (position < sentence_length) ?
      sentence->GetShapeId(position) : TOKEN_STOP; // Current shape.
     // Shape on the left.
    pSID = (position > 0) ?
      sentence->GetShapeId(position - 1) : TOKEN_START;
    // Shape on the right.
    nSID = (position < sentence_length - 1) ?
      sentence->GetShapeId(position + 1) : TOKEN_STOP;
    // Shape two positions on the left.
    ppSID = (position > 1) ?
      sentence->GetShapeId(position - 2) : TOKEN_START;
    // Shape two positions on the right.
    nnSID = (position < sentence_length - 2) ?
      sentence->GetShapeId(position + 2) : TOKEN_STOP;
  }

  // Maximum is 255 feature templates.
  CHECK_LT(EntityFeatureTemplateBigram::COUNT, 256);

  // Bias feature.
  fkey = encoder_.CreateFKey_NONE(EntityFeatureTemplateBigram::BIAS,
                                  flags);
  AddFeature(fkey, features, EntityFeatureTemplateBigram::BIAS);

  // Lexical features.
  if (feature_set_bitmap->test(0)) {
    //Word
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::W,
                                 WID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::W);
    //pContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::pW,
                                 pWID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::pW);
    //nContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::nW,
                                 nWID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::nW);
    //ppContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::ppW,
                                 ppWID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::ppW);
    //nnContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::nnW,
                                 nnWID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::nnW);
  }

  // POS features.
  if (feature_set_bitmap->test(1)) {
    //Word
    fkey = encoder_.CreateFKey_P(EntityFeatureTemplateBigram::P,
                                 PID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::P);
    //pContext
    fkey = encoder_.CreateFKey_PP(EntityFeatureTemplateBigram::PpP,
                                  PID, pPID,
                                  flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::PpP);
    //nContext
    fkey = encoder_.CreateFKey_PP(EntityFeatureTemplateBigram::PnP,
                                  PID, nPID,
                                  flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::PnP);
    //ppContext
    fkey = encoder_.CreateFKey_PPP(EntityFeatureTemplateBigram::PpPppP,
                                   PID, pPID, ppPID,
                                   flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::PpPppP);
    //nnContext
    fkey = encoder_.CreateFKey_PPP(EntityFeatureTemplateBigram::PnPnnP,
                                   PID, nPID, nnPID,
                                   flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::PnPnnP);
  }

  // Shape features.
  if (feature_set_bitmap->test(2)) {
    //Word
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::S,
                                 SID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::S);
    //pContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::pS,
                                 pSID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::pS);
    //nContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::nS,
                                 nSID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::nS);
    //ppContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::ppS,
                                 ppSID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::ppS);
    //nnContext
    fkey = encoder_.CreateFKey_W(EntityFeatureTemplateBigram::nnS,
                                 nnSID,
                                 flags);
    AddFeature(fkey, features, EntityFeatureTemplateBigram::nnS);
  }
}

void EntityFeatures::AddTrigramFeatures(EntityInstanceNumeric *sentence,
                                        int position) {
  CHECK(!input_features_trigrams_[position]) << position
    << " " << sentence->size();
  MultiBinaryFeatures *features = new MultiBinaryFeatures;
  input_features_trigrams_[position] = features;

  uint64_t fkey;
  uint8_t flags = 0x0;
  flags |= EntityFeatureTemplateParts::TRIGRAM;


  // Maximum is 255 feature templates.
  CHECK_LT(EntityFeatureTemplateTrigram::COUNT, 256);

  // Bias feature.
  fkey = encoder_.CreateFKey_NONE(EntityFeatureTemplateTrigram::BIAS,
                                  flags);
  AddFeature(fkey, features, EntityFeatureTemplateTrigram::BIAS);
}
