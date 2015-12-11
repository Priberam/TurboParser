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

#ifndef SEQUENCEPART_H_
#define SEQUENCEPART_H_

#include <stdio.h>
#include <vector>
#include "Part.h"

#define USE_MEMORY_POOl_FOR_SEQUENCE_PARTS 1

#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
#include "MemPool.h"
#endif

using namespace std;

enum {
  SEQUENCEPART_UNIGRAM = 0,
  SEQUENCEPART_BIGRAM,
  SEQUENCEPART_TRIGRAM,
  NUM_SEQUENCEPARTS
};

class SequencePartUnigram : public Part {
public:
  SequencePartUnigram() { position_ = tag_ = -1; };
  SequencePartUnigram(int position, int tag) :
    position_(position), tag_(tag) {};
  virtual ~SequencePartUnigram() {};

public:
  int position() { return position_; };
  int tag() { return tag_; };

public:
  int type() { return SEQUENCEPART_UNIGRAM; };

private:
  int position_; // Word position.
  int tag_; // Tag ID.
};

class SequencePartBigram : public Part {
public:
  SequencePartBigram() { position_ = tag_ = tag_left_ = -1; };
  SequencePartBigram(int position, int tag, int tag_left) :
    position_(position), tag_(tag), tag_left_(tag_left) {};
  virtual ~SequencePartBigram() {};

public:
  int position() { return position_; };
  int tag() { return tag_; };
  int tag_left() { return tag_left_; };

public:
  int type() { return SEQUENCEPART_BIGRAM; };

private:
  int position_; // Word position.
  int tag_; // Tag ID.
  int tag_left_; // Left tag ID
};

class SequencePartTrigram : public Part {
public:
  SequencePartTrigram() {
    position_ = tag_ = tag_left_ = tag_left_left_ = -1;
  };
  SequencePartTrigram(int position, int tag, int tag_left, int tag_left_left) :
    position_(position), tag_(tag), tag_left_(tag_left),
    tag_left_left_(tag_left_left) {};
  virtual ~SequencePartTrigram() {};

public:
  int position() { return position_; };
  int tag() { return tag_; };
  int tag_left() { return tag_left_; };
  int tag_left_left() { return tag_left_left_; };

public:
  int type() { return SEQUENCEPART_TRIGRAM; };

private:
  int position_; // Word position.
  int tag_; // Tag ID.
  int tag_left_; // Left tag ID
  int tag_left_left_; // Tag ID two words on the left.
};

class SequenceParts : public Parts {
public:
  SequenceParts()
#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
    : unigram_pool_(), bigram_pool_(), trigram_pool_()
#endif
  {};
  virtual ~SequenceParts() {
    DeleteAll();
  };
  void Initialize() {
    DeleteAll();
    for (int i = 0; i < NUM_SEQUENCEPARTS; ++i) {
      offsets_[i] = -1;
    }
  };

  Part *CreatePartUnigram(int position, int tag) {
#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
    // First, raw memory is requested from the memory pool
    SequencePartUnigram * get_allocated_part = (SequencePartUnigram *)unigram_pool_.GetNextBuffer();
    //Once that memory is granted, the new object is constructed in it, with a "placement new"
    //return new (get_allocated_part) SequencePartUnigram(position, tag);
    get_allocated_part->SequencePartUnigram::SequencePartUnigram(position, tag);
    return get_allocated_part;
#else
    return new SequencePartUnigram(position, tag);
#endif
  }
  Part *CreatePartBigram(int position, int tag, int tag_left) {
#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
    // First, raw memory is requested from the memory pool
    SequencePartBigram * get_allocated_part = (SequencePartBigram *)bigram_pool_.GetNextBuffer();
    //Once that memory is granted, the new object is constructed in it, with a "placement new"
    //return new (get_allocated_part) SequencePartBigram(position, tag, tag_left);
    get_allocated_part->SequencePartBigram::SequencePartBigram(position, tag, tag_left);
    return get_allocated_part;
#else
    return new SequencePartBigram(position, tag, tag_left);
#endif
  }
  Part *CreatePartTrigram(int position, int tag, int tag_left,
                          int tag_left_left) {
#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
    // First, raw memory is requested from the memory pool
    SequencePartTrigram * get_allocated_part = (SequencePartTrigram *)trigram_pool_.GetNextBuffer();
    //Once that memory is granted, the new object is constructed in it, with a "placement new"
    //return new (get_allocated_part) SequencePartTrigram(position, tag, tag_left, tag_left_left);
    get_allocated_part->SequencePartTrigram::SequencePartTrigram(position, tag, tag_left, tag_left_left);
    return get_allocated_part;
#else
    return new SequencePartTrigram(position, tag, tag_left, tag_left_left);
#endif
}

public:
  void DeleteAll();

public:
  void BuildUnigramIndices(int sentence_length);
  void BuildBigramIndices(int sentence_length);
  void BuildTrigramIndices(int sentence_length);
  void BuildIndices(int sentence_length);
  void DeleteUnigramIndices();
  void DeleteBigramIndices();
  void DeleteTrigramIndices();
  void DeleteIndices();
  const vector<int> &FindUnigramParts(int position) {
    return index_[position];
  }
  const vector<int> &FindBigramParts(int position) {
    return index_bigrams_[position];
  }
  const vector<int> &FindTrigramParts(int position) {
    return index_trigrams_[position];
  }

  // Set/Get offsets:
  void BuildOffsets() {
    for (int i = NUM_SEQUENCEPARTS - 1; i >= 0; --i) {
      if (offsets_[i] < 0) {
        offsets_[i] = (i == NUM_SEQUENCEPARTS - 1) ? size() : offsets_[i + 1];
      }
    }
  };
  void SetOffsetUnigram(int offset, int size) {
    SetOffset(SEQUENCEPART_UNIGRAM, offset, size);
  };
  void SetOffsetBigram(int offset, int size) {
    SetOffset(SEQUENCEPART_BIGRAM, offset, size);
  };
  void SetOffsetTrigram(int offset, int size) {
    SetOffset(SEQUENCEPART_TRIGRAM, offset, size);
  };
  void GetOffsetUnigram(int *offset, int *size) const {
    GetOffset(SEQUENCEPART_UNIGRAM, offset, size);
  };
  void GetOffsetBigram(int *offset, int *size) const {
    GetOffset(SEQUENCEPART_BIGRAM, offset, size);
  };
  void GetOffsetTrigram(int *offset, int *size) const {
    GetOffset(SEQUENCEPART_TRIGRAM, offset, size);
  };

private:
  // Get offset from part index.
  void GetOffset(int i, int *offset, int *size) const {
    *offset = offsets_[i];
    *size = (i < NUM_SEQUENCEPARTS - 1) ? offsets_[i + 1] - (*offset) :
      SequenceParts::size() - (*offset);
  }

  // Set offset from part index.
  void SetOffset(int i, int offset, int size) {
    offsets_[i] = offset;
    if (i < NUM_SEQUENCEPARTS - 1) offsets_[i + 1] = offset + size;
  }

private:
  vector<vector<int> >  index_;
  vector<vector<int> >  index_bigrams_;
  vector<vector<int> >  index_trigrams_;
  int offsets_[NUM_SEQUENCEPARTS];

#if USE_MEMORY_POOl_FOR_SEQUENCE_PARTS == 1
  MemPool<SequencePartUnigram> unigram_pool_;
  MemPool<SequencePartBigram> bigram_pool_;
  MemPool<SequencePartTrigram> trigram_pool_;
#endif
};

#endif /* SEQUENCEPART_H_ */
