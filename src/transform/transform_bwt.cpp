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

#include "transform_bwt.h"
#include "filenames.h"
#include <cassert>
#include <string.h>
#include <iostream>
//#include <chrono>
//#include <iomanip>

BwtCycler::BwtCycler( PreprocessFiles* filenames )
: fns( filenames )
{
    // Create buffers
    inBwtBuff = new uint8_t[BWT_BUFFER];
    outBwtBuff = new uint8_t[BWT_BUFFER];
    inInsBuff = new uint8_t[BWT_BUFFER];
    inEndBuff = new ReadId[IDS_BUFFER];
    outEndBuff = new ReadId[IDS_BUFFER];
    for ( int i ( 0 ); i < 4; i++ )
    {
        inIdsBuff[i] = new ReadId[IDS_BUFFER];
        outInsBuff[i] = new uint8_t[BWT_BUFFER];
        for ( int j ( 0 ); j < 5; j++ )
        {
            outIdsBuff[i][j] = new ReadId[IDS_BUFFER];
        }
    }
    inIdsBuff[4] = new ReadId[IDS_BUFFER];
    
    samePosFlag = (CharId)1 << 63;
    samePosMask = ~samePosFlag;
    sameByteFlag = (uint8_t)1 << 7;
    sameByteMask = ~sameByteFlag;
    
    isFinal = isPenultimate = false;
    
    FILE* bin = fns->getBinary( true, false );
    readEndBwt = writeEndBwt = readEndIds = writeEndIds = false;
    fseek( bin, 1, SEEK_SET );
    fread( &id, 8, 1, bin );
    fclose( bin );
    
    memset( &isRunArray, false, 256 );
    memset( &readArray[0], 0, 64 );
    memset( &readArray[64], 1, 64 );
    memset( &readArray[128], 2, 64 );
    memset( &readArray[192], 3, 64 );
    insMax1 = 256;
    insMax2 = insMax1 * insMax1;
    insMax4 = insMax2 * insMax2;
    sapMax1 = 256;
    sapMax2 = sapMax1 * 256;
    sapMax3 = sapMax2 * 256;
    for ( int i ( 0 ); i < 4; i++ )
    {
        for ( int j ( 0 ); j < 64; j++ )
        {
            runLenArray[ i * 64 + j] = j + 1;
        }
        isRunArray[ i * 64 + 63 ] = true;
        writeBaseBit[i] = 64 * i;
        writeMaxBase[i] = 63;
        writeFullByte[i] = 64 * i + 63;
    }
    
    for ( int i = 0; i < 8; i++ )
    {
        endBitArray[i] = 1 << ( 7 - i );
    }
}

BwtCycler::~BwtCycler()
{
    if ( inBwtBuff ) delete[] inBwtBuff;
    if ( outBwtBuff ) delete[] outBwtBuff;
    if ( inInsBuff ) delete[] inInsBuff;
    if ( inEndBuff ) delete[] inEndBuff;
    if ( outEndBuff ) delete[] outEndBuff;
    for ( int i ( 0 ); i < 4; i++ )
    {
        if ( inIdsBuff[i] ) delete[] inIdsBuff[i];
        if ( outInsBuff[i] ) delete[] outInsBuff[i];
        for ( int j ( 0 ); j < 5; j++ )
        {
            if ( outIdsBuff[i][j] ) delete[] outIdsBuff[i][j];
        }
    }
    if ( inIdsBuff[4] ) delete[] inIdsBuff[4];
}

void BwtCycler::finish( uint8_t cycle )
{
    fns->setCyclerFinal( inBwt, outBwt, inEnd, outEnd, cycle );
    prepIn();
    prepOutFinal();
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fns->setCyclerFinalIter( inIns, inIds[4], cycle, i );
        prepIter();
        finishIter( i );
    }
    
    // Flush buffers and close write files
    writeLast();
    fclose( inBwt );
    fwrite( outBwtBuff, 1, pOutBwt, outBwt );
    fclose( outBwt );
    fwrite( outEndBuff, 4, pOutEnd, outEnd );
    fclose( outEnd );
    outBwt = fns->getReadPointer( fns->bwt, true );
    uint8_t bwtBegin = 57;
    fwrite( &bwtBegin, 1, 1, outBwt );
    fwrite( &id, 8, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts[4], 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    fclose( outBwt );
}
//
void BwtCycler::finishIter( uint8_t i )
{
    if ( insLeft ) readNextPos();
    
    while( insLeft )
    {
        if ( currPos == nextPos )
        {
            writeEnd();
            if ( --insLeft ) readNextPos();
        }
        else
        {
            writeBwt();
        }
    }
    
    nextPos = -1;
    if ( currSplit )
    {
        writeRun( splitChar, splitRun );
        currSplit = false;
        if ( splitChar == 4 ) rewriteEnd( splitRun );
    }
    
    while ( currPos < charSizes[i] )
    {
        writeBwt();
    }
    assert( currPos >= charSizes[i] );
    currPos -= charSizes[i];
    
    fclose( inIns );
    fclose( inIds[4] );
}

void BwtCycler::flush( uint8_t cycle )
{
    // Flush buffers and close write files
    writeLast();
    fwrite( outBwtBuff, 1, pOutBwt, outBwt );
    fclose( outBwt );
    fwrite( outEndBuff, 4, pOutEnd, outEnd );
    fclose( outEnd );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        writeInsBuff( i );
        fclose( outIns[i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            fwrite( outIdsBuff[i][j], 4, pOutIds[i][j], outIds[i][j] );
            fclose( outIds[i][j] );
        }
    }
    
    // Edit in counts
    fns->setCyclerUpdate( outBwt, outEnd, outIns, outIds, cycle );
    fseek( outBwt, 9, SEEK_SET );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    fclose( outBwt );
    fwrite ( &endCount, 4, 1, outEnd );
    fclose( outEnd );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fwrite( &insCounts[i], 8, 1, outIns[i] );
        fclose( outIns[i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            fwrite( &idsCounts[i][j], 4, 1, outIds[i][j] );
            fclose( outIds[i][j] );
        }
    }
}

void BwtCycler::prepIn()
{
    // Read sizes for this cycle
    bool doReadBwtEnds;
    CharId thisId;
    fread( &thisId, 8, 1, inBwt );
    if ( thisId != id )
    {
        cerr << "Error: temporary files do not appear to be from the same session." << endl;
        exit( EXIT_FAILURE );
    }
    fread( &doReadBwtEnds, 1, 1, inBwt );
    fread( &bwtLeft, 8, 1, inBwt );
    fread( &charSizes, 8, 5, inBwt );
    fread( &basePos, 4, 4, inBwt );
    fread( &endLeft, 4, 1, inEnd );
    if ( doReadBwtEnds )
    {
        if ( !readEndBwt ) setReadEnds();
        if ( !writeEndBwt ) setWriteEnds();
    }
    
    // Reset pointers and counts for cycle
    bwtFirst = true;
    currSplit = false;
    lastChar = -1;
    bwtCount = currPos = 0;
    pInBwt = BWT_BUFFER;
    pInEnd = IDS_BUFFER;
    pOutBwt = pOutEnd = 0;
    for ( int i ( 0 ); i < 4; i++ )
    {
        charCounts[i] = basePos[i];
    }
    charCounts[4] = 0;
    memset( lastIns, 0, 32 );
    memset( inSapCount, 0, 20 );
    memset( outSapCount, 0, 20 );
}

void BwtCycler::prepIter()
{
    fread( &insLeft, 8, 1, inIns );
    pInIns = BWT_BUFFER;
    nextPos = insLeft ? 0 : -1;
    
    for ( int j ( 0 ); j < 5; j++ )
    {
        if ( isFinal && j < 4 ) continue;
        fread( &idsLeft[j], 4, 1, inIds[j] );
        pInIds[j] = IDS_BUFFER;
    }
}

void BwtCycler::prepOut()
{
    anyEnds = false;
    if ( !chars )
    {
        isPenultimate = true;
        nextChar = 4;
    }
    else if ( ends )
    {
        anyEnds = true;
        writeEndIds = true;
        if ( !writeEndBwt ) setWriteEnds();
    }
    endCount = 0;
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        insCounts[i] = 0;
        pOutIns[i] = 0;
        fwrite( &insCounts[i], 8, 1, outIns[i] );
        
        for ( int j ( 0 ); j < 5; j++ )
        {
            idsCounts[i][j] = 0;
            pOutIds[i][j] = 0;
            fwrite( &idsCounts[i][j], 4, 1, outIds[i][j] );
        }
    }
    
    // Set place holders for this cycle size
    fwrite( &id, 8, 1, outBwt );
    fwrite( &writeEndBwt, 1, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 5, outBwt );
    fwrite( &basePos, 4, 4, outBwt );
    fwrite( &endCount, 4, 1, outEnd );
}

void BwtCycler::prepOutFinal()
{
    isFinal = true;
    if ( !writeEndBwt ) setWriteEnds();
    
    endCount = 0;
    uint8_t bwtBegin = 57, idsBegin = 9;
    CharId finalEndCount = basePos[0] + basePos[1] + basePos[2] + basePos[3];
    fwrite( &bwtBegin, 1, 1, outBwt );
    fwrite( &id, 8, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &finalEndCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    
    fwrite( &idsBegin, 1, 1, outEnd );
    fwrite( &id, 8, 1, outEnd );
    
    memset( &charCounts, 0, 40 );
    for ( int i ( 0 ); i < 4; i++ )
    {
        if ( basePos[i] )
        {
            writeRun( i, basePos[i] );
        }
    }
    currPos = 0;
}

void BwtCycler::readBwtIn()
{
    fread( inBwtBuff, 1, min( BWT_BUFFER, bwtLeft), inBwt );
    pInBwt = 0;
}

void BwtCycler::readInsertBuff()
{
    if ( pInIns < BWT_BUFFER )
    {
        CharId diff = BWT_BUFFER - pInIns;
        memcpy( inInsBuff, &inInsBuff[pInIns], diff );
        fread( &inInsBuff[diff], 1, min( insLeft - 1, BWT_BUFFER ) - diff, inIns );
    }
    else
    {
        fread( inInsBuff, 1, min( insLeft - 1, BWT_BUFFER ), inIns );
    }
    pInIns = 0;
}

void BwtCycler::readNextId()
{
    // Refresh IDs buffer if necessary
    if ( pInIds[thisChar] == IDS_BUFFER )
    {
        fread( inIdsBuff[thisChar], 4, min( idsLeft[thisChar], IDS_BUFFER ), inIds[thisChar] );
        pInIds[thisChar] = 0;
    }

    // Set id and next character
    nextId = inIdsBuff[thisChar][ pInIds[thisChar]++ ];
    if ( anyEnds && ( ends[ (nextId / 2) / 8 ] & endBitArray[ (nextId / 2) % 8 ] ) )
    {
        nextChar = 4;
    }
    else if ( !isPenultimate )
    {
        nextChar = byteToInt[ nextId & 0x3 ][ chars[ nextId / 4 ] ];
    }
    --idsLeft[thisChar];
}

void BwtCycler::readNextPos()
{
    // Refresh inserts buffer if necessary
    if ( pInIns == BWT_BUFFER )
    {
        fread( inInsBuff, 1, min( insLeft, BWT_BUFFER ), inIns );
        pInIns = 0;
    }
    
    thisChar = inInsBuff[pInIns++];
    nextSame = thisChar & uint8_t(128);
    uint8_t posBytes = ( thisChar >> 3 ) & 0x7;
    if ( pInIns == BWT_BUFFER ) readInsertBuff();
    --insLeft;
    CharId thisPos = inInsBuff[pInIns++];
    while ( posBytes-- )
    {
        if ( pInIns == BWT_BUFFER ) readInsertBuff();
        --insLeft;
        thisPos = ( thisPos << 8 ) ^ inInsBuff[pInIns++];
    }
    nextPos += thisPos;
    
    if ( !nextSame )
    {
        thisChar &= 0x7;
        if ( !isFinal ) readNextId();
    }
    
    if ( currSplit && currPos < nextPos )
    {
        if ( currPos + splitRun > nextPos )
        {
            ReadId thisRun = nextPos - currPos;
            splitRun -= thisRun;
            if ( splitChar == 4 ) rewriteEnd( thisRun );
            writeRun( splitChar, thisRun );
        }
        else
        {
            if ( splitChar == 4 ) rewriteEnd( splitRun );
            writeRun( splitChar, splitRun );
            currSplit = false;
        }
    }
}

void BwtCycler::run( uint8_t* inChars, uint8_t* inEnds, uint8_t cycle )
{
    fns->setCycler( inBwt, outBwt, inEnd, outEnd, outIns, outIds, cycle );
    chars = inChars;
    ends = inEnds;
    
    double cycleStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();

    prepIn();
    prepOut();
    
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fns->setCyclerIter( inIns, inIds, cycle, i );
        prepIter();
        runIter( i );
    }
    
    assert( !bwtLeft );
    flush( cycle );
    fclose( inBwt );
    fclose( inEnd );
    
//    cout << std::fixed << std::setprecision(2) << " cycle: " << ( clock() - cycleStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << std::flush;
}

void BwtCycler::rewriteEnd( ReadId runLen )
{
    assert( runLen );
    while ( runLen-- )
    {
        if ( pInEnd == IDS_BUFFER )
        {
            fread( inEndBuff, 4, min( endLeft, IDS_BUFFER ), inEnd );
            pInEnd = 0;
        }
        if ( pOutEnd == IDS_BUFFER )
        {
            fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
            pOutEnd = 0;
        }
        outEndBuff[ pOutEnd++ ] = inEndBuff[ pInEnd++ ];
        ++endCount;
        --endLeft;
    }
}

void BwtCycler::runIter( uint8_t i )
{
    if ( insLeft ) readNextPos();
    
    while ( insLeft )
    {
        if ( currPos == nextPos )
        {
            if ( nextSame )
            {
                writeSame();
            }
            else
            {
                writeNext();
            }
            
            if ( --insLeft ) readNextPos();
        }
        else
        {
            writeBwt();
        }
    }
    
    nextPos = -1;
    if ( currSplit )
    {
        writeRun( splitChar, splitRun );
        currSplit = false;
        if ( splitChar == 4 ) rewriteEnd( splitRun );
    }
    
    while ( currPos < charSizes[i] )
    {
        writeBwt();
    }
    assert( currPos >= charSizes[i] );
    currPos -= charSizes[i];
    
    fclose( inIns );
    for ( int j ( 0 ); j < 4; j++ )
    {
        fclose( inIds[j] );
    }
}

void BwtCycler::setReadEnds()
{
    memset( &isRunArray, false, 256 );
    memset( &readArray[0], 0, 63 );
    memset( &readArray[63], 1, 63 );
    memset( &readArray[126], 2, 63 );
    memset( &readArray[189], 3, 63 );
    memset( &readArray[252], 4, 4 );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        runLenArray[ 252 + i ] = i + 1;
        for ( int j ( 0 ); j < 63; j++ )
        {
            runLenArray[ i * 63 + j] = j + 1;
        }
        isRunArray[ i * 63 + 62 ] = true;
    }
    isRunArray[255] = true;
    
    readEndBwt = true;
}

void BwtCycler::setWriteEnds()
{
    for ( int i ( 0 ); i < 4; i++ )
    {
        writeBaseBit[i] = 63 * i;
        writeMaxBase[i] = 62;
        writeFullByte[i] = 63 * i + 62;
    }
    writeBaseBit[4] = 252;
    writeMaxBase[4] = 3;
    writeFullByte[4] = 255;
    
    writeEndBwt = true;
}

void BwtCycler::writeIdsToFile( uint8_t i, uint8_t j )
{
    fwrite( outIdsBuff[i][j], 4, pOutIds[i][j], outIds[i][j] );
    pOutIds[i][j] = 0;
}

void BwtCycler::writeBwt()
{
    if ( pInBwt == BWT_BUFFER ) readBwtIn();
    uint8_t currChar = inBwtBuff[pInBwt];
    uint8_t c = readArray[ currChar ];
    ReadId runLen = runLenArray[ currChar ];
    --bwtLeft;
    ++pInBwt;
    
    // Count run length if greater than 63
    if ( isRunArray[currChar] )
    {
        if ( pInBwt == BWT_BUFFER ) readBwtIn();
        ReadId thisRun = 0;
        uint8_t shiftCount = 0;
        do {
            if ( pInBwt == BWT_BUFFER ) readBwtIn();
            thisRun ^= ( inBwtBuff[pInBwt] & sameByteMask ) << ( shiftCount++ * 7 );
        } while ( inBwtBuff[pInBwt++] & sameByteFlag );
        
        runLen += thisRun;
        bwtLeft -= shiftCount;
    }
    
    if ( currPos + runLen > nextPos )
    {
        splitChar = c;
        splitRun = currPos + runLen - nextPos;
        runLen -= splitRun;
        currSplit = true;
    }
    
    writeRun( c, runLen );
    if ( c == 4 ) rewriteEnd( runLen );
}

void BwtCycler::writeBwtByte( uint8_t c )
{
    if ( pOutBwt == BWT_BUFFER )
    {
        fwrite( outBwtBuff, 1, BWT_BUFFER, outBwt );
        pOutBwt = 0;
    }
    
    outBwtBuff[ pOutBwt++ ] = c;
    ++bwtCount;
}

void BwtCycler::writeEnd()
{
    ReadId thisSap = 1;
//    bool isDupe = nextSame;
    if ( nextSame )
    {
        if ( pInIns == BWT_BUFFER ) readInsertBuff();
        --insLeft;
        thisSap = inInsBuff[pInIns++];
        for ( uint8_t j = thisChar & 0x3; j--; )
        {
            if ( pInIns == BWT_BUFFER ) readInsertBuff();
            --insLeft;
            thisSap = ( thisSap << 8 ) ^ inInsBuff[pInIns++];
        }
    }
    
    writeRun( 4, thisSap );
    
//    ReadId dupeCount = 1;
//    ReadId dupeArray[thisSap + 1];
//    if ( isDupe )
//    {
//        dupeArray[0] = thisSap;
//    }
    
    while ( thisSap-- )
    {
        if ( pInIds[4] == IDS_BUFFER )
        {
            fread( inIdsBuff[4], 4, min( IDS_BUFFER, idsLeft[4] ), inIds[4] );
            pInIds[4] = 0;
        }
        if ( pOutEnd == IDS_BUFFER )
        {
            fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
            pOutEnd = 0;
        }
        
//        dupeArray[dupeCount++] = inIdsBuff[4][ pInIds[4] ];
        outEndBuff[ pOutEnd++ ] = inIdsBuff[4][ pInIds[4]++ ];
    }
    
//    if ( isDupe )
//    {
//        fwrite( &dupeArray, 4, dupeCount, dupes );
//    }
}

void BwtCycler::writeInsBuff( uint8_t i )
{
    fwrite( outInsBuff[i], 1, pOutIns[i], outIns[i] );
    insCounts[i] += pOutIns[i];
    pOutIns[i] = 0;
}

void BwtCycler::writeInsBytes()
{
    CharId ins = charCounts[thisChar] - lastIns[thisChar];
    lastIns[thisChar] = charCounts[thisChar];
    
    uint8_t posBytes = ins > 255;
    if ( posBytes )
    {
        ReadId insRemain = ins >> 16;
        while ( insRemain )
        {
            posBytes++;
            insRemain >>= 8;
        }
        
        outInsBuff[thisChar][ pOutIns[thisChar]++ ] ^= ( posBytes << 3 );
        for ( int i = posBytes + 1; i--; )
        {
            if ( pOutIns[thisChar] == BWT_BUFFER ) writeInsBuff( thisChar );
            outInsBuff[thisChar][ pOutIns[thisChar]++ ] = ins >> ( 8 * i );
        }
    }
    else
    {
        if ( ++pOutIns[thisChar] == BWT_BUFFER ) writeInsBuff( thisChar );
        outInsBuff[thisChar][ pOutIns[thisChar]++ ] = ins;
    }
}

void BwtCycler::writeLast()
{
    if ( bwtFirst )
    {
        bwtFirst = false;
        return;
    }
    
    uint8_t maxBase = writeMaxBase[lastChar];
    if ( lastRun < maxBase )
    {
        writeBwtByte( writeBaseBit[lastChar] + lastRun );
    }
    // Encode a new byte for every power of 128
    else
    {
        writeBwtByte( writeFullByte[lastChar] );
        lastRun -= maxBase;
        while ( lastRun >= 128 )
        {
            writeBwtByte( sameByteFlag ^ ( lastRun & sameByteMask ) );
            lastRun >>= 7;
        }
        writeBwtByte( lastRun );
    }
}

void BwtCycler::writeNext()
{
    if ( thisChar != 4 )
    {
        // Write Insert buffer to file if full
        if ( pOutIns[thisChar] == BWT_BUFFER ) writeInsBuff( thisChar );
        outInsBuff[thisChar][ pOutIns[thisChar] ] = nextChar;
        writeInsBytes();
    }
    
    // Write id to IDs stack
    writeNextId();
    writeRun( thisChar, 1 );
}

void BwtCycler::writeNextId()
{
    if ( thisChar == 4 )
    {
        if ( pOutEnd == IDS_BUFFER )
        {
            fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
            pOutEnd = 0;
        }
        
        outEndBuff[ pOutEnd++ ] = nextId;
        ++endCount;
    }
    else
    {
        // Write IDs buffer to file if full
        if ( pOutIds[thisChar][nextChar] == IDS_BUFFER )
        {
            fwrite( outIdsBuff[thisChar][nextChar], 4, IDS_BUFFER, outIds[thisChar][nextChar] );
            pOutIds[thisChar][nextChar] = 0;
        }

        // Write id to buffer
        outIdsBuff[thisChar][nextChar][ pOutIds[thisChar][nextChar]++ ] = nextId;
        idsCounts[thisChar][nextChar]++;
    }
}

void BwtCycler::writeRun( uint8_t c, ReadId runLen )
{
    if ( c == lastChar )
    {
        lastRun += runLen;
    }
    else
    {
        writeLast();
        lastChar = c;
        lastRun = runLen - 1;
    }
    
    charCounts[c] += runLen;
    currPos += runLen;
}

void BwtCycler::writeSame()
{
    uint8_t baseBytes = ( thisChar & 0x3 );
    if ( thisChar & 0x4 )
    {
        if ( pInIns + 1 + baseBytes > BWT_BUFFER ) readInsertBuff();
        insLeft -= 1 + baseBytes;
        thisChar = 4;
        ReadId thisSap = inInsBuff[pInIns++];
        for ( uint8_t j = baseBytes; j--; )
        {
            thisSap = ( thisSap << 8 ) ^ inInsBuff[pInIns++];
        }
        writeRun( thisChar, thisSap );
        while ( thisSap-- )
        {
            readNextId();
            writeNextId();
        }
    }
    
    if ( pInIns + ( 1 + baseBytes ) * 4 > BWT_BUFFER ) readInsertBuff();
    insLeft -= ( 1 + baseBytes ) * 4;
    for ( int i = 0; i < 4; i++ )
    {
        thisChar = i;
        ReadId thisSap = inInsBuff[pInIns++];
        for ( uint8_t j = baseBytes; j--; )
        {
            thisSap = ( thisSap << 8 ) ^ inInsBuff[pInIns++];
        }
        
        if ( thisSap > 1 )
        {
            // Write IDs and count same runs
            memset( &outSapCount, 0, 20 );
            for ( ReadId j = 0; j < thisSap; j++ )
            {
                readNextId();
                writeNextId();
                outSapCount[nextChar]++;
            }
            ReadId maxCount = max( outSapCount[0], max( outSapCount[1], max( outSapCount[2], max( outSapCount[3], outSapCount[4] ) ) ) );
            uint8_t thisBytes = maxCount < sapMax1 ? 0 : ( maxCount < sapMax2 ? 1 : ( maxCount < sapMax3 ? 2 : 3 ) );
            
            if ( pOutIns[i] == BWT_BUFFER ) writeInsBuff( i );
            outInsBuff[i][ pOutIns[i] ] = 128 ^ ( outSapCount[4] ? 4 : 0 ) ^ thisBytes;
            writeInsBytes();
            writeRun( thisChar, thisSap );
            
            // Write end runs to insert buffer
            if ( isPenultimate || outSapCount[4] )
            {
                for ( uint8_t j = thisBytes + 1; j--; )
                {
                    if ( pOutIns[i] == BWT_BUFFER ) writeInsBuff( i );
                    outInsBuff[i][ pOutIns[i]++ ] = outSapCount[4] >> ( 8 * j );
                }
            }
            
            // Write same runs to insert buffer
            if ( !isPenultimate )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    for ( uint8_t k = thisBytes + 1; k--; )
                    {
                        if ( pOutIns[i] == BWT_BUFFER ) writeInsBuff( i );
                        outInsBuff[i][ pOutIns[i]++ ] = outSapCount[j] >> ( 8 * k );
                    }
                }
            }
        }
        else if ( thisSap )
        {
            readNextId();
            writeNext();
        }
    }
}
