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
//#include <chrono>
//#include <iomanip>

BinaryReader::BinaryReader( PreprocessFiles* filenames )
: fns( filenames )
{
    endCount = 0;
    bin = fns->getReadPointer( fns->bin, false );
    fread( &seqsBegin, 1, 1, bin );
    fread( &id, 8, 1, bin );
    fread( &readLen, 1, 1, bin );
    fread( &cycle, 1, 1, bin );
    fread( &revComp, 1, 1, bin );
    fseek( bin, 16, SEEK_SET );
    fread( &seqCount, 4, 1, bin );
    
    lineLen = 1 + ( readLen + 3 ) / 4;
    fileSize = (CharId)seqCount * (CharId)lineLen;
    buffSize = 16777216 - ( 16777216 % lineLen );
    charSize = ( seqCount + 3 ) / 4;
    buff = new uint8_t[buffSize];
    chars = new uint8_t[ ( seqCount + 3 ) / 4 ];
    ends = NULL;
    
    for ( int i = 0; i < 8; i++ )
    {
        endBitArray[i] = 1 << ( 7 - i );
    }
    if ( !cycle && !seqCount )
    {
        cerr << "Error: cannot resume as no cycles have been previously completed" << endl;
        exit( EXIT_FAILURE );
    }
    
    prep();
    if ( cycle < 2 )
    {
        init();
    }
    if ( cycle >= readLen + 1 )
    {
        cerr << "Error: transformation has already been completed." << endl;
        exit( EXIT_FAILURE );
    }
    chr = fns->getReadPointer( fns->tmpChr, false );
    trm = fns->getReadPointer( fns->tmpTrm, false );
    
    CharId trimSkip = trmBegin;
    if ( cycle >= minTrim ) for ( int i = 0; i < ( cycle+1-minTrim ); i++ ) trimSkip += trimCounts[i]*4;
    assert( !trimCounts.empty() || minTrim == readLen );
    if ( !trimCounts.empty() ) fseek( trm, trimSkip, SEEK_SET );
    if ( cycle > 2 ) fseek( chr, ( cycle-2 ) * charSize, SEEK_SET );
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

void BinaryReader::init()
{
    fseek( bin, seqsBegin, SEEK_SET );
    cycle = 1;
    FILE* ids[4][4];
    for ( int i( 0 ); i < 4; i++ ) for ( int j( 0 ); j < 4; j++ )
    {
        ids[i][j] = fns->getReadPointer( fns->tmpIds[0][i][j], true );
        fseek( ids[i][j], 4, SEEK_SET );
    }
    
    uint8_t line[lineLen];
    ReadId p = 0, i = 0;
    ReadId idsCounts[4][4]{0};
    
    {
        ReadId blockSize = 8*1000, pChar = 0;
        CharId seeks[readLen];
        uint8_t seq[readLen],* outs[readLen];
        for ( uint8_t j = 0; j < readLen; j++ ) outs[j] = new uint8_t[blockSize]{0};
        for ( uint8_t j = 3; j < readLen; j++ ) seeks[j] = ( j-3 ) * charSize;
        
        chr = fns->getWritePointer( fns->tmpChr );
        fseek( chr, CharId( readLen-3 ) * charSize - 1, SEEK_SET );
        fwrite( line, 1, 1, chr );
        fclose( chr );
        chr = fns->getReadPointer( fns->tmpChr, true );
        
        // Prep trim file
        int pTrim[readLen-minTrim]{0};
        ReadId* bufTrim[readLen-minTrim], totalTrims = 0;
        CharId fpTrims[readLen-minTrim];
        for ( uint8_t j = 0; j+minTrim < readLen; j++ )
        {
            fpTrims[j] = j ? fpTrims[j-1] + ( trimCounts[j-1] * 4 ) : trmBegin;
            bufTrim[j] = new ReadId[1000];
            totalTrims += trimCounts[j];
        }
        
        trm = fns->getReadPointer( fns->tmpTrm, true );
        fseek( trm, CharId( totalTrims ) * 4 - 1 + trmBegin, SEEK_SET );
        fwrite( line, 1, 1, trm );

        for ( ReadId id = 0; id < seqCount; id++ )
        {
            if ( i == 4 )
            {
                ++pChar;
                if ( ++p == blockSize )
                {
                    for ( uint8_t j = 3; j < readLen; j++ )
                    {
                        fseek( chr, seeks[j], SEEK_SET );
                        fwrite( outs[j], 1, p, chr );
                        memset( outs[j], 0, p );
                        seeks[j] += p;
                    }
                    p = 0;
                }
                i = 0;
            }
            fread( line, 1, lineLen, bin );
            
            if ( line[0] < readLen )
            {
                assert( line[0] >= minTrim );
                uint8_t j = line[0] - minTrim;
                
                if ( pTrim[j] == 1000 )
                {
                    fseek( trm, fpTrims[j], SEEK_SET );
                    fwrite( bufTrim[j], 4, pTrim[j], trm );
                    fpTrims[j] += pTrim[j]*4;
                    pTrim[j] = 0;
                }
                
                bufTrim[j][ pTrim[j]++ ] = revComp ? ( id / 2 ) : id;
            }
            
            uint8_t base = line[0]-1, nxt = line[0]-3;
            for ( uint8_t j = nxt; j < line[0]; j++ ) seq[j] = byteToInt[ j&0x3 ][ line[1+j/4] ];
            for ( uint8_t j = 0; j < nxt; j++ )
            {
                seq[j] = byteToInt[ j&0x3 ][ line[1+j/4] ];
                if ( !i ) outs[base-j][p] = intToByte[0][ seq[j] ];
                else outs[base-j][p] |= intToByte[i][ seq[j] ];
            }
            if ( !i ) chars[pChar] = intToByte[0][ seq[nxt] ];
            else chars[pChar] |= intToByte[i][ seq[nxt] ];
            i++;
            
            idsCounts[ seq[base] ][ seq[base-1] ]++;
            fwrite( &id, 4, 1, ids[ seq[base] ][ seq[base-1] ] );

            if ( !revComp ) continue;
            
            ++id;
            for ( uint8_t j = 3; j < line[0]; j++ ) outs[j][p] |= intToByte[i][ 3-seq[j] ];
            chars[pChar] |= intToByte[i++][ 3-seq[2] ];
            idsCounts[ 3-seq[0] ][ 3-seq[1] ]++;
            fwrite( &id, 4, 1, ids[ 3-seq[0] ][ 3-seq[1] ] );
        }
        
        if ( i ? ++p : p ) for ( uint8_t j = 3; j < readLen; j++ )
        {
            fseek( chr, seeks[j], SEEK_SET );
            fwrite( outs[j], 1, p, chr );
            seeks[j] += p;
        }
        
        for ( uint8_t j = 0; j+minTrim < readLen; j++ ) if ( pTrim[j] )
        {
            fseek( trm, fpTrims[j], SEEK_SET );
            fwrite( bufTrim[j], 4, pTrim[j], trm );
            fpTrims[j] += pTrim[j]*4;
        }
        
        for ( uint8_t j = 0; j < readLen; j++ ) delete outs[j];
        for ( uint8_t j = 0; j+minTrim < readLen; j++ ) delete bufTrim[j];
        
        fclose( chr );
        fclose( trm );
    }
    
    for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ )
    {
        fseek( ids[i][j], 0, SEEK_SET );
        fwrite( &idsCounts[i][j], 4, 1, ids[i][j] );
        fclose( ids[i][j] );
    }
    
    // Set inserts
    for ( int i = 0; i < 4; i++ )
    {
        FILE* ins = fns->getReadPointer( fns->tmpIns[0][i], true );
        
        CharId maxCount = max( idsCounts[i][0], max( idsCounts[i][1], max( idsCounts[i][2], idsCounts[i][3] ) ) );
        CharId thisMax = 255;
        uint8_t sapByte = 0;
        while ( maxCount > thisMax )
        {
            thisMax <<= 8;
            sapByte++;
        }
        
        CharId insCount = 2 + ( ( sapByte + 1 ) * 4 );
        fwrite( &insCount, 8, 1, ins );
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
        fwrite( &writeBuff, 1, insCount, ins );
        fclose( ins );
    }
    
    // Set base BWT
    FILE* bwt = fns->getWritePointer( fns->tmpBwt[0] );
    ReadId basePos[4]{0};
    CharId charCounts[5]{0};
    for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ ) basePos[i] += idsCounts[i][j];
    for ( int i = 0; i < 4; i++ ) charCounts[i] = basePos[i];
    bool writeEndBwt = false;
    CharId bwtCount = 0;
    fwrite( &id, 8, 1, bwt );
    fwrite( &writeEndBwt, 1, 1, bwt );
    fwrite( &bwtCount, 8, 1, bwt );
    fwrite( &charCounts, 8, 5, bwt );
    fwrite( &basePos, 4, 4, bwt );
    fclose( bwt );
    
    FILE* ends = fns->getWritePointer( fns->tmpEnd[0] );
    ReadId endcount = 0;
    fwrite( &endcount, 4, 1, ends );
    fclose( ends );
    
//    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << endl;
}

void BinaryReader::prep()
{
    trm = fns->getReadPointer( fns->tmpTrm, true );
    fread( &trmBegin, 2, 1, trm );
    fread( &minTrim, 1, 1, trm );
    ReadId inTrim;
    for ( uint8_t j = 0; j+minTrim < readLen; j++ )
    {
        fread( &inTrim, 4, 1, trm );
        trimCounts.push_back( inTrim );
    }
    fclose( trm );
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
    if ( cycle == 2 ) return;
    
    double readStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
    
    fread( chars, 1, charSize, chr );
    
    if ( anyEnds = ( cycle >= minTrim ) )
    {
        unordered_set<ReadId> ids;
        for( ReadId in : {110803679,112732249,90490184,97199286,25558601,152543493,155562156,136293585,152034118,55435202,106260666,50802570,95519741,4171125,56612370,121848256,9567765,139030520,112612339,72400370,122625440,35006278,139506091,10585004,143005037,145942559,151231624,159486086} ) ids.insert( in / 2 );
        if ( !ends ) ends = new uint8_t[ ( seqCount + 15 ) / 16 ]{0};
        int i = cycle - minTrim;
        ReadId id;
        for ( ReadId j = 0; j < trimCounts[i]; j++ )
        {
            fread( &id, 4, 1, trm );
            ends[id/8] |= endBitArray[id % 8];
            endCount++;
        }
    }
    
//    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << flush;
}

void BinaryReader::update()
{
    FILE* fp = fns->getReadPointer( fns->bin, true );
    fseek( fp, 10, SEEK_SET );
    fwrite( &cycle, 1, 1, fp );
    fclose( fp );
}

BinaryWriter::BinaryWriter( PreprocessFiles* filenames, uint8_t inLibCount, uint8_t inReadLen, bool revComp )
: fns( filenames ), libCount( inLibCount ), readLen( inReadLen ), readLens( inReadLen, 0 ), libCounts( NULL ), revComp( revComp )
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
//    fns->setBinaryWrite( bin, bwt, ends, ins, ids );
    bin = fns->getBinary( false, false );
    lineLen = 1 + ( readLen + 3 ) / 4;
    buffSize = 16777216 - ( 16777216 % lineLen );
    binBuff = new uint8_t[buffSize];
    
    if ( libCount ) libCounts = new ReadId[libCount]{0};
    seqsBegin = 21 + ( libCount * 12 );
    
    for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ ) charPlaceCounts[i][j].resize( readLen, 0 );
    
    uint8_t dummy8 = 0, revCal = revComp ? 2 : 0;
    uint16_t dummy16 = 0;
    uint32_t dummy32 = 0;
    
    fwrite( &seqsBegin, 1, 1, bin );             // Byte offset of first sequence
    fwrite( &id, 8, 1, bin );                    // ID number for this transform session
    fwrite( &readLen, 1, 1, bin );               // Read length
    fwrite( &cycle, 1, 1, bin );                 // Current cycles complete
    fwrite( &revCal, 1, 1, bin );                // Is calibrated
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
    fclose( bin );
    
    // Fill in missing binary variables
    bin = fns->getBinary( true, true );
    if ( revComp )
    {
        seqCount *= 2;
        vector<ReadId>revCounts[4][4];
        for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ ) revCounts[i][j].resize( readLen, 0 );
        for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ ) for ( int k = 1; k < readLen; k++ ) revCounts[i][j][k] = charPlaceCounts[3-j][3-i][readLen-k];
        for ( int i = 0; i < 4; i++ ) for ( int j = 0; j < 4; j++ ) for ( int k = 1; k < readLen; k++ ) charPlaceCounts[i][j][k] += revCounts[i][j][k];
    }
    
    uint8_t dummy;
    fseek( bin, seqsBegin, SEEK_SET );
    fread( &dummy, 1, 1, bin );
    assert( dummy );
    
    fseek( bin, 16, SEEK_SET );
    fwrite( &seqCount, 4, 1, bin );
    fseek( bin, 21, SEEK_SET );
    for ( int i ( 0 ); i < libCount; i++ )
    {
        fwrite( &libCounts[i], 4, 1, bin );
        fseek( bin, 8, SEEK_CUR );
    }
    fclose( bin );
    
    // Write counts to trim file
    FILE* trm = fns->getWritePointer( fns->tmpTrm );
    uint8_t minReadLen = readLen;
    for ( uint8_t i = 0; i < readLen; i++ ) if ( readLens[i] ) minReadLen = min( minReadLen, i );
    assert( minReadLen );
    uint16_t trimBegin = 3 + ( ( readLen-minReadLen ) * 4 );
    fwrite( &trimBegin, 2, 1, trm );
    fwrite( &minReadLen, 1, 1, trm );
    for ( uint8_t i = minReadLen; i < readLen; i++ ) fwrite( &readLens[i], 4, 1, trm );
    fclose( trm );
     
    // Set ids bucket limits
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            ReadId limit = 0;
            for ( int k = 1; k < readLen; k++ ) limit = max( limit, charPlaceCounts[i][j][k] );
            for ( int k = 1; k < readLen; k++ ) limit = max( limit, charPlaceCounts[3-j][3-i].end()[-k-1] );
            for ( int k : { 0, 1 } )
            {
                FILE* fp = fns->getWritePointer( fns->tmpIds[k][i][j] );
                fseek( fp, limit*4, SEEK_SET );
                fwrite( &limit, 4, 1, fp );
                fclose( fp );
            }
        }
        for ( int k : { 0, 1 } )
        {
            ReadId dummy = 0;
            FILE* fp = fns->getWritePointer( fns->tmpIds[k][i][4] );
            fwrite( &dummy, 4, 1, fp );
            fclose( fp );
        }
    }
    
    // Set ins bucket limits
    for ( int i = 0; i < 4; i++ )
    {
        ReadId limit = 0;
        for ( int k = 0; k < readLen; k++ )
        {
            ReadId counted = 0;
            for ( int j = 0; j < 4; j++ ) counted += charPlaceCounts[j][i][k];
            limit = max( limit, counted );
        }
        for ( int k : { 0, 1 } )
        {
            FILE* fp = fns->getWritePointer( fns->tmpIns[k][i] );
            fseek( fp, limit*4+8, SEEK_SET );
            fwrite( &limit, 4, 1, fp );
            fclose( fp );
        }
    }
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
    ReadId thisCount = revComp ? seqCount*2 : seqCount;
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
    // Check and write sequence length into one byte
    uint8_t line[lineLen];
    line[0] = read.length();
    if ( line[0] > readLen )
    {
        read = read.substr( 0, readLen );
        line[0] = readLen;
        cerr << "Error: Unexpectedly long read of length " << to_string( line[0] ) << " given set length of " << to_string( readLen ) << "." << endl;
        exit( EXIT_FAILURE );
    }
    
    // Encode characters into 2 bits per byte
    CharId p = 0;
    uint8_t l = 0;
    for ( uint8_t j ( 0 ); j < line[0]; j++ )
    {
        uint8_t c = charToInt[ read[j] ];
        assert( read[j] != 'N' );
        charPlaceCounts[l][c][j]++;
        uint8_t i = j & 0x3;
        l = c;
        if ( !i ) line[++p] = intToByte[i][c];
        else line[p] += intToByte[i][c];
    }
    fwrite( line, 1, lineLen, bin );
    seqCount++;
    readLens[ line[0] ]++;
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
        fclose( ids[i][4] );
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

