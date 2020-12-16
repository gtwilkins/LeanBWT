/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the LeanBWT software package <https://github.com/gtwilkins/LeanBWT>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INDEX_READER_H
#define INDEX_READER_H

#include "types.h"
#include "filenames.h"
#include "index_structs.h"

class IndexReader
{
public:
    IndexReader( Filenames* fns );
    ~IndexReader();
    
    void countRange( uint8_t i, CharId rank, CharId count, CharCount &ranks, CharCount &counts );
    void countRange( uint8_t i, CharId rank, CharId edge, CharId count, CharCount &ranks, CharCount &edges, CharCount &counts );
    void createSeeds( string &fn, int mer );
    int primeOverlap( uint8_t* q, CharId &rank, CharId &count );
    void primeOverlap( string &seq, vector<uint8_t> &q, CharId &rank, CharId &count, int &ol, bool drxn );
    int setBaseAll( vector<uint8_t> &q, CharId &rank, CharId &count );
    void setBaseAll( uint8_t i, uint8_t j, CharId &rank, CharId &count );
    void setBaseAll( uint8_t i, uint8_t j, CharId &rank, CharId &edge, CharId &count );
    ReadId setBaseMap( uint8_t i, uint8_t j, CharId &rank, CharId &count );
    void setBaseOverlap( uint8_t i, uint8_t j, CharId &rank, CharId &count );
    
private:
    void advance( CharCount &ranks, CharId &bwtIndex, CharId &toCount );
    void createSeeds( FILE* fp, int i, int it, int limit, CharId rank, CharId edge, CharId count );
    void setRank( uint8_t i, CharId rank, CharCount &ranks );
    CharId setRankIndex( CharId rankIndex, CharCount &ranks );
    
    
    FILE* bwt;
    uint8_t* buff,* mers;
    
    CharId bwtSize, indexSize, markSize;
    ReadId bwtPerIndex, indexPerMark;
    int kmerLen;
    uint8_t beginBwt, beginIdx;
    uint8_t sizePerIndex;
    
    // Index data
    uint8_t* index_;
    ReadId* marks_;
    CharId charRanks[4], charCounts[5];
    CharId baseCounts[5][4], midRanks[4][4];
    
    bool isBaseRun[256];
    uint8_t decodeBaseChar[256], decodeBaseRun[256];
    uint8_t runFlag, runMask;
};

#endif /* INDEX_READER_H */

