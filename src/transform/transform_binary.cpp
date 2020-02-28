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

#include "transform_binary.h"
#include "filenames.h"
#include <cassert>
#include <iostream>
#include <string.h>

BinaryReader::BinaryReader( PreprocessFiles* filenames )
: fns( filenames )
{
    endCount = 0;
    bin = fns->getReadPointer( fns->bin, false );
    fread( &seqsBegin, 1, 1, bin );
    fread( &id, 8, 1, bin );
    fread( &readLen, 1, 1, bin );
    fread( &cycle, 1, 1, bin );
    fseek( bin, 16, SEEK_SET );
    fread( &seqCount, 4, 1, bin );
    
    lineLen = 1 + ( readLen + 3 ) / 4;
    fileSize = (CharId)seqCount * (CharId)lineLen;
    buffSize = 16777216 - ( 16777216 % lineLen );
    buff = new uint8_t[buffSize];
    chars = new uint8_t[ ( seqCount + 3 ) / 4 ];
    ends = new uint8_t[ ( seqCount + 15 ) / 16 ];
    
    for ( int i = 0; i < 8; i++ )
    {
        endBitArray[i] = 1 << ( 7 - i );
    }
    
    if ( !cycle )
    {
        cout << "Error: cannot resume as no cycles have been previously completed" << endl;
        exit( EXIT_FAILURE );
    }
    if ( cycle >= readLen + 1 )
    {
        cout << "Error: transformation has already been completed." << endl;
        exit( EXIT_FAILURE );
    }
}

BinaryReader::~BinaryReader()
{
    if ( chars ) delete chars;
    if ( buff ) delete buff;
    if ( ends ) delete ends;
    chars = buff = ends = NULL;
}

void BinaryReader::finish()
{
    ++cycle;
    assert( cycle == readLen + 1 );
    update();
}

void BinaryReader::read()
{
    prevEndCount = endCount;
    if ( ++cycle == readLen )
    {
        delete[] chars;
        chars = NULL;
        return;
    }
    
    anyEnds = false;
    memset( chars, 0, ( seqCount + 3 ) / 4 );
    memset( ends, 0, ( seqCount + 15 ) / 16 );
    
    fseek( bin, seqsBegin, SEEK_SET );
    CharId readsLeft = seqCount / 2, fileLeft = fileSize;
    CharId pBuff = buffSize, pChar = 0, pFwd, pRev = 1 + cycle / 4;
    ReadId readCount = 0;
    uint8_t iChar = 0, iFwd, iRev = cycle & 0x3;
    
    while ( readsLeft )
    {
//        if ( pBuff == buffSize ) rebuffBin( bin, buff, buffSize, fileLeft, pBuff );
        if ( pBuff == buffSize )
        {
            CharId thisBuff = min( buffSize, fileLeft );
            fread( buff, 1, thisBuff, bin );
            fileLeft -= thisBuff;
            pBuff = 0;
        }
//        if ( !iChar ) chars[pChar] = 0;
        if ( buff[pBuff] == cycle )
        {
            ends[readCount/8] ^= endBitArray[readCount % 8];
            anyEnds = true;
            endCount++;
        }
        pFwd = 1 + ( buff[pBuff] - 1 - cycle ) / 4;
        iFwd = ( buff[pBuff] - 1 - cycle ) & 0x3;
        
        chars[pChar] += intToByte[iChar++][ byteToInt[iFwd][ buff[ pBuff + pFwd ] ] ];
        chars[pChar] += intToByte[iChar][ complement[ byteToInt[iRev][ buff[ pBuff + pRev ] ] ] ];
        iter( pChar, iChar );
        pBuff += lineLen;
        readsLeft--;
        readCount++;
    }
}

void BinaryReader::update()
{
    FILE* fp = fns->getReadPointer( fns->bin, true );
    fseek( fp, 10, SEEK_SET );
    fwrite( &cycle, 1, 1, fp );
    fclose( fp );
}

BinaryWriter::BinaryWriter( PreprocessFiles* filenames, uint8_t inLibCount, uint8_t inReadLen )
: fns( filenames ), libCount( inLibCount ), readLen( inReadLen ), libCounts( NULL )
{
    pBin = 0;
    seqCount = 0;
    cycle = currLib = 0;
    
    srand( time(NULL) );
    CharId randMask = 65535;
    id = ( (CharId)( rand() & randMask ) << 48 ) 
            ^ ( (CharId)( rand() & randMask ) << 32 ) 
            ^ ( (CharId)( rand() & randMask ) << 16 ) 
            ^ (CharId)( rand() & randMask );
    
    memset( &charCounts, 0, 40 );
    fns->setBinaryWrite( bin, bwt, ends, ins, ids );
    lineLen = 1 + ( readLen + 3 ) / 4;
    buffSize = 16777216 - ( 16777216 % lineLen );
    binBuff = new uint8_t[buffSize];
    if ( libCount ) libCounts = new ReadId[libCount]{0};
    seqsBegin = 21 + ( libCount * 12 );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        for ( int j( 0 ); j < 4; j++ )
        {
            idsCounts[i][j] = pIds[i][j] = 0;
            idsBuff[i][j] = new ReadId[IDS_BUFFER];
            fwrite( &idsCounts[i][j], 4, 1, ids[i][j] );
        }
    }
    
    uint8_t dummy8 = 0;
    uint16_t dummy16 = 0;
    uint32_t dummy32 = 0;
    
    fwrite( &seqsBegin, 1, 1, bin );             // Byte offset of first sequence
    fwrite( &id, 8, 1, bin );                    // ID number for this transform session
    fwrite( &readLen, 1, 1, bin );               // Read length
    fwrite( &cycle, 1, 1, bin );                 // Current cycles complete
    fwrite( &dummy8, 1, 1, bin );                // Is calibrated
    fwrite( &dummy32, 4, 1, bin );               // Estimated coverage
    fwrite( &seqCount, 4, 1, bin );              // Sequence count for each library, total sequence count
    fwrite( &libCount, 1, 1, bin );              // Number of libraries
    for ( int i ( 0 ); i < libCount; i++ )
    {
        fwrite( &dummy32, 4, 1, bin );           // Library sequence count
        fwrite( &dummy16, 2, 3, bin );           // Library insert size estimates
        fwrite( &dummy8, 1, 2, bin );            // Library type details
    }
    
}

BinaryWriter::~BinaryWriter()
{
    delete binBuff;
    if ( libCount ) delete[] libCounts;
    for ( int i ( 0 ); i < 4; i++ )
    {
        for ( int j ( 0 ); j < 4; j++ )
        {
            delete[] idsBuff[i][j];
        }
    }
}

void BinaryWriter::close()
{
    assert( libCount == currLib && cycle == 0 );
    
    // Clear remaining binary buffer
    dumpBin();
    fclose( bin );
    
    // Fill in missing binary variables
    bin = fns->getBinary( true, true );
    cycle = 1;
    fseek( bin, 10, SEEK_SET );
    fwrite( &cycle, 1, 1, bin );
    fseek( bin, 16, SEEK_SET );
    fwrite( &seqCount, 4, 1, bin );
    fseek( bin, 21, SEEK_SET );
    for ( int i ( 0 ); i < libCount; i++ )
    {
        fwrite( &libCounts[i], 4, 1, bin );
        fseek( bin, 8, SEEK_CUR );
    }
    fclose( bin );
    
    // Write BWT, POS and SAP, fill in counts for BWT, POS, SAP and IDS
    writeBwt();
    writeEnd();
    writeIns();
    writeIds();
}

void BinaryWriter::dumpBin()
{
    fwrite( binBuff, 1, pBin, bin );
    pBin = 0;
}

void BinaryWriter::dumpIds( uint8_t i, uint8_t j )
{
    fwrite( idsBuff[i][j], 4, pIds[i][j], ids[i][j] );
    pIds[i][j] = 0;
}

void BinaryWriter::setNextLibrary()
{
    assert( currLib < libCount );
    ReadId thisCount = seqCount;
    for ( int i ( 0 ); i < currLib; i++ )
    {
        assert( libCount <= thisCount );
        thisCount -= libCounts[i];
    }
    libCounts[currLib] = thisCount;
    currLib++;
}

void BinaryWriter::write( string &read )
{
    if ( pBin == buffSize ) dumpBin();
    
    // Check and write sequence length into one byte
    uint8_t seqLen = read.length();
    binBuff[pBin] = seqLen;
    if ( seqLen > readLen )
    {
        read = read.substr( 0, readLen );
        seqLen = readLen;
//        cerr << "Error: Unexpectedly long read of length " << to_string( seqLen ) << " given set length of " << to_string( readLen ) << "." << endl;
//        exit( EXIT_FAILURE );
    }
    
    // Encode characters into 2 bits per byte
    CharId p = pBin + 1;
    uint8_t i = 0;
    for ( uint8_t j ( 0 ); j < seqLen; j++ )
    {
        assert( read[j] != 'N' );
        write2Bit( binBuff, p, i, charToInt[ read[j] ] );
    }
    pBin += lineLen;
    
    // Get first and second characters from each end
    uint8_t c[4];
    c[0] = charToInt[ read[seqLen - 1] ];
    c[1] = charToInt[ read[seqLen - 2] ];
    c[2] = complement[ charToInt[ read[0] ] ];
    c[3] = complement[ charToInt[ read[1] ] ];
    
    // Write first cycle of BWT
    for ( int k : { 0, 2 } )
    {
        if ( pIds[ c[k] ][ c[k+1] ] == IDS_BUFFER ) dumpIds( c[k], c[k+1] );
        idsBuff[ c[k] ][ c[k+1] ][ pIds[ c[k] ][ c[k+1] ]++ ] = seqCount++;
        idsCounts[ c[k] ][ c[k+1] ]++;
        charCounts[ c[k] ]++;
    }
}

void BinaryWriter::writeBwt()
{
    ReadId basePos[4];
    for ( int i ( 0 ); i < 4; i++ )
    {
        basePos[i] = charCounts[i];
    }
    bool writeEndBwt = false;
    CharId bwtCount = 0;
    fwrite( &id, 8, 1, bwt );
    fwrite( &writeEndBwt, 1, 1, bwt );
    fwrite( &bwtCount, 8, 1, bwt );
    fwrite( &charCounts, 8, 5, bwt );
    fwrite( &basePos, 4, 4, bwt );
    fclose( bwt );
}

void BinaryWriter::writeEnd()
{
    ReadId endcount = 0;
    fwrite( &endcount, 4, 1, ends );
    fclose( ends );
}

void BinaryWriter::writeIds()
{
    for ( int i( 0 ); i < 4; i++ )
    {
        for ( int j( 0 ); j < 4; j++ )
        {
            dumpIds( i, j );
            fclose( ids[i][j] );
            ids[i][j] = fns->getReadPointer( fns->tmpIds[0][i][j], true );
            fwrite( &idsCounts[i][j], 4, 1, ids[i][j] );
            fclose( ids[i][j] );
        }
    }
}

void BinaryWriter::writeIns()
{
    for ( int i ( 0 ); i < 4; i++ )
    {
        CharId maxCount = max( idsCounts[i][0], max( idsCounts[i][1], max( idsCounts[i][2], idsCounts[i][3] ) ) );
        CharId thisMax = 255;
        uint8_t sapByte = 0;
        while ( maxCount > thisMax )
        {
            thisMax <<= 8;
            sapByte++;
        }
        
        CharId insCount = 2 + ( ( sapByte + 1 ) * 4 );
        fwrite( &insCount, 8, 1, ins[i] );
        uint8_t writeBuff[ insCount ];
        writeBuff[0] = 128 ^ sapByte;
        writeBuff[1] = 0;
        uint8_t p = 2;
        for ( int j = 0; j < 4; j++ )
        {
            for ( uint8_t k = sapByte + 1; k--; )
            {
                writeBuff[p++] = ( idsCounts[i][j] >> ( 8 * k ) ) & uint8_t(255);
            }
        }
        fwrite( &writeBuff, 1, insCount, ins[i] );
        fclose( ins[i] );
    }
}

