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

BinaryReader::BinaryReader( PreprocessFiles* filenames, bool revComp )
: fns( filenames ), revComp( revComp )
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
        cout << "Error: cannot resume as no cycles have been previously completed" << endl;
        exit( EXIT_FAILURE );
    }
    if ( !cycle )
    {
        init();
    }
    if ( cycle >= readLen + 1 )
    {
        cout << "Error: transformation has already been completed." << endl;
        exit( EXIT_FAILURE );
    }
    chr = fns->getReadPointer( fns->tmpChr, false );
    if ( cycle > 3 ) assert( false );
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
        uint8_t seq[readLen],* outs[readLen], base = readLen-1, nxt = readLen-3;
        for ( uint8_t j = 0; j < readLen; j++ ) outs[j] = new uint8_t[blockSize];
        for ( uint8_t j = 3; j < readLen; j++ ) seeks[j] = ( j-3 ) * charSize;
        
        chr = fns->getWritePointer( fns->tmpChr );
        fseek( chr, CharId( readLen-3 ) * charSize - 1, SEEK_SET );
        fwrite( line, 1, 1, chr );
        fclose( chr );
        chr = fns->getReadPointer( fns->tmpChr, true );

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
                        seeks[j] += p;
                    }
                    p = 0;
                }
                i = 0;
            }
            fread( line, 1, lineLen, bin );
            for ( uint8_t j = nxt; j < readLen; j++ ) seq[j] = byteToInt[ j&0x3 ][ line[1+j/4] ];
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
            for ( uint8_t j = 3; j < readLen; j++ ) outs[j][p] |= intToByte[i][ 3-seq[j] ];
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
        
        for ( uint8_t j = 0; j < readLen; j++ ) delete outs[j];
        
        fclose( chr );
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

//void BinaryReader::read()
//{
//    prevEndCount = endCount;
//    if ( ++cycle == readLen )
//    {
//        delete[] chars;
//        chars = NULL;
//        return;
//    }
//    anyEnds = false;
//    if ( cycle == 2 ) return;
//    
//    fseek( bin, seqsBegin, SEEK_SET );
//    ReadId readCount = seqCount / ( 1 + revComp ), p = -1, i = 0, pFwd = 1 + ( readLen - 1 - cycle ) / 4, pRev = 1 + cycle / 4;
//    uint8_t line[lineLen], c, iFwd = 1 + ( readLen - 1 - cycle ) / 4, iRev = cycle & 0x3;
//    
//    double readStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
//    
//    while ( readCount-- )
//    {
//        fread( line, 1, lineLen, bin );
//        c = byteToInt[ ( line[0] - 1 - cycle ) & 0x3 ][ line[ 1 + ( line[0] - 1 - cycle ) / 4 ] ];
//        
//        if ( !i ) chars[++p] = intToByte[0][ c ];
//        else chars[p] += intToByte[i][ c ];
//        if ( ++i == 4 ) i = 0;
//        
//        if ( revComp )
//        {
//            assert( false );
//            c = complement[ byteToInt[iRev][ line[pRev] ] ];
//            if ( !i ) chars[++p] = intToByte[0][ c ];
//            else chars[p] += intToByte[i][ c ];
//            if ( ++i == 4 ) i = 0;
//        }
//    }
//    
//    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << flush;
//}

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
    
    anyEnds = false;
    fread( chars, 1, charSize, chr );
    
//    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << flush;
}

//void BinaryReader::read()
//{
//    prevEndCount = endCount;
//    if ( ++cycle == readLen )
//    {
//        delete[] chars;
//        chars = NULL;
//        return;
//    }
//    if ( cycle == 2 ) return;
//    
//    anyEnds = false;
//    memset( chars, 0, ( seqCount + 3 ) / 4 );
//    if ( ends ) memset( ends, 0, ( seqCount + 15 ) / 16 );
//    
//    fseek( bin, seqsBegin, SEEK_SET );
//    CharId readsLeft = seqCount / ( 1 + revComp ), fileLeft = fileSize;
//    CharId pBuff = buffSize, pChar = 0, pFwd, pRev = 1 + cycle / 4;
//    ReadId readCount = 0;
//    uint8_t iChar = 0, iFwd, iRev = cycle & 0x3;
//    
//    double readStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
//    
//    while ( readsLeft )
//    {
////        if ( pBuff == buffSize ) rebuffBin( bin, buff, buffSize, fileLeft, pBuff );
//        if ( pBuff == buffSize )
//        {
//            CharId thisBuff = min( buffSize, fileLeft );
//            fread( buff, 1, thisBuff, bin );
//            fileLeft -= thisBuff;
//            pBuff = 0;
//        }
////        if ( !iChar ) chars[pChar] = 0;
//        if ( buff[pBuff] == cycle )
//        {
//            if ( !ends ) ends = new uint8_t[ ( seqCount + 15 ) / 16 ]{0};
//            ends[readCount/8] ^= endBitArray[readCount % 8];
//            anyEnds = true;
//            endCount++;
//        }
//        pFwd = 1 + ( buff[pBuff] - 1 - cycle ) / 4;
//        iFwd = ( buff[pBuff] - 1 - cycle ) & 0x3;
//        
//        chars[pChar] += intToByte[iChar][ byteToInt[iFwd][ buff[ pBuff + pFwd ] ] ];
//        if ( revComp ) chars[pChar] += intToByte[++iChar][ complement[ byteToInt[iRev][ buff[ pBuff + pRev ] ] ] ];
//        iter( pChar, iChar );
//        pBuff += lineLen;
//        readsLeft--;
//        readCount++;
//    }
//    
//    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << flush;
//}

void BinaryReader::update()
{
    FILE* fp = fns->getReadPointer( fns->bin, true );
    fseek( fp, 10, SEEK_SET );
    fwrite( &cycle, 1, 1, fp );
    fclose( fp );
}

//void BinaryReader::test()
//{
//    fseek( bin, seqsBegin, SEEK_SET );
//    
//    string prefix = "/home/glen/Genomes/SpBwt/sp_block";
//    vector<string> ffns;
//    for ( int i = 0; i < 10; ) ffns.push_back( prefix + "-"+ to_string( ++i ) );
//    
//    double tStart = clock();
//    auto t_start = chrono::high_resolution_clock::now();
//    ifstream ifs( fns->bin, ios::binary );
//    
//    buffSize = 10*1000*lineLen;
//    CharId fileLeft = fileSize, pBuff = buffSize;
//    CharId counted = 0;
//    CharId blockCount = 100000000;
//    CharId blockSize = blockCount*lineLen;
//    CharId blocks = 0;
//    char* buf = new char[buffSize];
//    
////    while ( fileLeft )
////    {
////        double blockStart = clock();
////        auto b_start = chrono::high_resolution_clock::now();
//////        ofstream ofs( ffns[blocks++], ios::binary );
////        FILE* ofp = fopen( ffns[blocks].c_str(), "wb" );
////        fseek( ofp, min( fileLeft, blockSize )-1, SEEK_SET );
////        fwrite( buf, 1, 1, ofp );
////        fclose( ofp );
////        ofp = fopen( ffns[blocks++].c_str(), "rb+" );
////        
////        for ( CharId j = blockSize/buffSize; fileLeft && j--; )
////        {
////            CharId thisBuff = min( buffSize, fileLeft );
////            ifs.read( buf, thisBuff );
////            fwrite( buf, 1, thisBuff, ofp );
//////            ofs.write( buf, thisBuff );
////            fileLeft -= thisBuff;
////        }
//////        ofs.close();
////        fclose( ofp );
////        cout << std::fixed << std::setprecision(2) << counted << " read: " << ( clock() - blockStart ) / CLOCKS_PER_SEC << " vs " << ( ( chrono::high_resolution_clock::now() - b_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
////    }
////    ifs.close();
////    cout << std::fixed << std::setprecision(2) << counted << " read: " << ( clock() - tStart ) / CLOCKS_PER_SEC << " vs " << ( ( chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
//    
//    double iStart = clock();
//    auto i_start = chrono::high_resolution_clock::now();
//    
//    fileLeft = fileSize;
//    
//    while ( fileLeft )
//    {
//        double blockStart = clock();
//        auto b_start = chrono::high_resolution_clock::now();
//        ifstream ifs( ffns[blocks++], ios::binary );
//        
//        for ( CharId j = blockSize/buffSize; fileLeft && j--; )
//        {
//            CharId thisBuff = min( buffSize, fileLeft );
//            ifs.read( buf, thisBuff );
//            fileLeft -= thisBuff;
//        }
//        ifs.close();
//        cout << std::fixed << std::setprecision(2) << counted << " read: " << ( clock() - blockStart ) / CLOCKS_PER_SEC << " vs " << ( ( chrono::high_resolution_clock::now() - b_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
//    }
//    
//    cout << std::fixed << std::setprecision(2) << counted << " read: " << ( clock() - iStart ) / CLOCKS_PER_SEC << " vs " << ( ( chrono::high_resolution_clock::now() - i_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
//    
////    for ( CharId i = seqCount; i--; )
////    {
////        if ( pBuff == buffSize )
////        {
////            CharId thisBuff = min( buffSize, fileLeft );
////            ifs.read( buf, thisBuff );
//////            fread( buff, 1, thisBuff, bin );
////            fileLeft -= thisBuff;
////            pBuff = 0;
////        }
////        
////        if ( !( ++counted % 10000000 ) )
////        {
////            
////            cout << std::fixed << std::setprecision(2) << counted << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
////        }
////        
////        pBuff += lineLen;
////    }
//////    
////    cout << std::fixed << std::setprecision(2) << " read: " << ( clock() - readStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << " ... " << endl;
//    exit( 1 );
//}

BinaryWriter::BinaryWriter( PreprocessFiles* filenames, uint8_t inLibCount, uint8_t inReadLen, bool revComp )
: fns( filenames ), libCount( inLibCount ), readLen( inReadLen ), libCounts( NULL ), revComp( revComp )
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
    // old
//    for ( int i ( 0 ); i < 4; i++ )
//    {
//        for ( int j( 0 ); j < 4; j++ )
//        {
//            idsCounts[i][j] = pIds[i][j] = 0;
//            idsBuff[i][j] = new ReadId[IDS_BUFFER];
//            fwrite( &idsCounts[i][j], 4, 1, ids[i][j] );
//        }
//    }
    
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
//    fseek( bin, 6811999995-1, SEEK_SET );
//    fwrite( &dummy8, 1, 1, bin );
//    fclose( bin );
//    bin = fopen( fns->bin.c_str(), "rb+" );
//    fseek( bin, seqsBegin, SEEK_SET );
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
//    cycle = 0;
//    fseek( bin, 10, SEEK_SET );
//    fwrite( &cycle, 1, 1, bin );
    fseek( bin, 16, SEEK_SET );
    fwrite( &seqCount, 4, 1, bin );
    fseek( bin, 21, SEEK_SET );
    for ( int i ( 0 ); i < libCount; i++ )
    {
        fwrite( &libCounts[i], 4, 1, bin );
        fseek( bin, 8, SEEK_CUR );
    }
    fclose( bin );
     
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

//void BinaryWriter::close()
//{
//    assert( libCount == currLib && cycle == 0 );
//    
//    // Clear remaining binary buffer
//    dumpBin();
//    fclose( bin );
//    
//    // Fill in missing binary variables
//    bin = fns->getBinary( true, true );
//    cycle = 1;
//    fseek( bin, 10, SEEK_SET );
//    fwrite( &cycle, 1, 1, bin );
//    fseek( bin, 16, SEEK_SET );
//    fwrite( &seqCount, 4, 1, bin );
//    fseek( bin, 21, SEEK_SET );
//    for ( int i ( 0 ); i < libCount; i++ )
//    {
//        fwrite( &libCounts[i], 4, 1, bin );
//        fseek( bin, 8, SEEK_CUR );
//    }
//    fclose( bin );
//    
//    // Write BWT, POS and SAP, fill in counts for BWT, POS, SAP and IDS
//    writeBwt();
//    writeEnd();
//    writeIns();
//    writeIds();
//}

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
//    if ( pBin == buffSize ) dumpBin();
    
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
}

//void BinaryWriter::write( string &read )
//{
//    if ( pBin == buffSize ) dumpBin();
//    
//    // Check and write sequence length into one byte
//    uint8_t seqLen = read.length();
//    binBuff[pBin] = seqLen;
//    if ( seqLen > readLen )
//    {
//        read = read.substr( 0, readLen );
//        seqLen = readLen;
////        cerr << "Error: Unexpectedly long read of length " << to_string( seqLen ) << " given set length of " << to_string( readLen ) << "." << endl;
////        exit( EXIT_FAILURE );
//    }
//    
//    // Encode characters into 2 bits per byte
//    CharId p = pBin + 1;
//    uint8_t i = 0;
//    for ( uint8_t j ( 0 ); j < seqLen; j++ )
//    {
//        assert( read[j] != 'N' );
//        write2Bit( binBuff, p, i, charToInt[ read[j] ] );
//    }
//    pBin += lineLen;
//    
//    // Get first and second characters from each end
//    uint8_t c[4];
//    c[0] = charToInt[ read[seqLen - 1] ];
//    c[1] = charToInt[ read[seqLen - 2] ];
//    c[2] = complement[ charToInt[ read[0] ] ];
//    c[3] = complement[ charToInt[ read[1] ] ];
//    
//    // Write first cycle of BWT
//    for ( int k = 0; k < 2 + revComp; k += 2 )
//    {
//        if ( pIds[ c[k] ][ c[k+1] ] == IDS_BUFFER ) dumpIds( c[k], c[k+1] );
//        idsBuff[ c[k] ][ c[k+1] ][ pIds[ c[k] ][ c[k+1] ]++ ] = seqCount++;
//        idsCounts[ c[k] ][ c[k+1] ]++;
//        charCounts[ c[k] ]++;
//    }
//}

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

