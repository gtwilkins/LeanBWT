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

#ifndef TRANSFORM_BINARY_H
#define TRANSFORM_BINARY_H

#include "types.h"
#include "constants.h"
#include "transform_structs.h"
#include "transform_functions.h"
#include "transform_bwt.h"

struct BinaryReader
{
    BinaryReader( PreprocessFiles* filenames );
    ~BinaryReader();
    
    void finish();
    void read();
    void update();
    
    PreprocessFiles* fns;
    FILE* bin;
    
    uint8_t* chars,* buff,* ends;
    
    CharId id;
    uint8_t endBitArray[8];
    
    uint8_t seqsBegin, lineLen, cycle, readLen;
    CharId buffSize, fileSize;
    CharId endCount, prevEndCount;
    ReadId seqCount;
    bool anyEnds;
};

struct BinaryWriter
{
    BinaryWriter( PreprocessFiles* filenames, uint8_t inLibCount, uint8_t inReadLen );
    ~BinaryWriter();
    
    void close();
    void dumpBin();
    void dumpIds( uint8_t i, uint8_t j );
    void setNextLibrary();
    void write( string &read );
    void writeBwt();
    void writeEnd();
    void writeIds();
    void writeIns();
    
    // File pointers
    PreprocessFiles* fns;
    FILE* bin,* bwt,* ends,* ins[4],* ids[4][4];
    
    // Buffers
    uint8_t* binBuff;
    ReadId* idsBuff[4][4];
    
    // Buffer pointers
    CharId pBin, pIds[4][4];
    uint8_t iBwt[4];
    
    CharId id;
    CharId buffSize;
    CharId charCounts[5];
    ReadId idsCounts[4][4];
    ReadId seqCount,* libCounts;
    uint8_t lineLen, readLen, currLib, libCount, seqsBegin, cycle;
};


#endif /* TRANSFORM_BINARY_H */

