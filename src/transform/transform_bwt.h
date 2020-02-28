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

#ifndef TRANFORM_BWT_H
#define TRANFORM_BWT_H

#include "types.h"
#include "filenames.h"
#include "transform_constants.h"
#include "transform_functions.h"

struct BwtCycler
{
public:
    BwtCycler( PreprocessFiles* filenames );
    ~BwtCycler();
    
    void run( uint8_t* inChars, uint8_t* inEnds, uint8_t cycle );
    void finish( uint8_t cycle );
    
private:
    void finishIter( uint8_t i );
    void flush( uint8_t cycle );
    void prepIn();
    void prepIter();
    void prepOut();
    void prepOutFinal();
    void readIds();
    void readBwtIn();
    void readInsertBuff();
    void readNextId();
    void readNextPos();
    void readNextSap();
    void rewriteEnd( ReadId runLen );
    void runIter( uint8_t i );
    void setReadEnds();
    void setWriteEnds();
    void writeIdsToFile( uint8_t i, uint8_t j );
    void writeBwt();
    void writeBwtByte( uint8_t c );
    void writeEnd();
    void writeInsBuff( uint8_t i );
    void writeInsBytes();
    void writeLast();
    void writeNext();
    void writeNextId();
    void writeRun( uint8_t c, ReadId runLen );
    void writeSame();
    void writeSplit();
    
    // Files
    PreprocessFiles* fns;
    FILE* inBwt,* outBwt;
    FILE* inIns,* outIns[4];
    FILE* inIds[5],* outIds[4][5];
    FILE* inEnd,* outEnd;
//    FILE* dupes;
    CharId id;
    
    // Buffers
    uint8_t* chars,* ends;
    uint8_t* inBwtBuff,* outBwtBuff;
    uint8_t* inInsBuff,* outInsBuff[4];
    ReadId* inIdsBuff[5],* outIdsBuff[4][5];
    ReadId* inEndBuff,* outEndBuff;
    
    // Buffer pointers
    ReadId pInBwt, pOutBwt;
    ReadId pInIns, pOutIns[4];
    ReadId pInIds[5], pOutIds[4][5];
    ReadId pInEnd, pOutEnd;
    
    // Counts
    CharId bwtCount, bwtLeft;
    CharId charCounts[5], charSizes[5];
    CharId insLeft, lastIns[4], insCounts[4];
    ReadId endLeft, idsLeft[5];
    ReadId inSapCount[5], outSapCount[5];
    ReadId endCount;
    CharId idsCounts[4][5];
    
    uint64_t inSap8[5], outSap8[5];
    uint32_t inSap4[5], outSap4[5];
    uint16_t inSap2[5], outSap2[5];
    uint8_t inSap1[5], outSap1[5];
    
    uint8_t sapSize;
    
    bool nextSame, bwtFirst, currSplit;
    uint8_t thisChar, nextChar, lastChar, splitChar;
    ReadId nextId, lastRun, splitRun;
    CharId currPos, nextPos;
    ReadId basePos[4];
    
    uint8_t sameByteFlag, sameByteMask;
    CharId samePosFlag, samePosMask;
    CharId insMax1, insMax2, insMax4;
    ReadId sapMax1, sapMax2, sapMax3;
    
    // Variables defining whether ends are currently being input/output
    bool anyEnds, isPenultimate, isFinal;
    bool readEndBwt, writeEndBwt;
    bool readEndIds, writeEndIds;
    
    // Tables for byte encoding
    bool isRunArray[256];
    uint8_t readArray[256], runLenArray[256], endBitArray[8];
    uint8_t writeFullByte[5], writeMaxBase[5], writeBaseBit[5];
};

#endif /* TRANFORM_BWT_H */

