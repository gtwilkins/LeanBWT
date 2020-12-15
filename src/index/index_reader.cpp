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

#include "index_reader.h"
#include "index_structs.h"
#include <cassert>
#include <string.h>
#include <iostream>
#include "constants.h"

IndexReader::IndexReader( Filenames* fns )
{
    FILE* bin,* idx,* mer;
    assert( fns );
    fns->setIndex( bin, bwt, idx, mer );
    CharId binId, bwtId, idxId;
    fseek( bin, 1, SEEK_SET );
    fread( &binId, 8, 1, bin );
    fclose( bin );
    
    fread( &beginBwt, 1, 1, bwt );
    fread( &bwtId, 8, 1, bwt );
    fread( &bwtSize, 8, 1, bwt );
    fread( &charCounts[4], 8, 1, bwt );
    fread( &charCounts, 8, 4, bwt );
    
    fread( &beginIdx, 1, 1, idx );
    fread( &idxId, 8, 1, idx );
    fread( &bwtPerIndex, 4, 1, idx );
    fread( &indexPerMark, 4, 1, idx );
    
    if ( binId != bwtId || binId != idxId )
    {
        cerr << "Error: disagreement among data files. They may be corrupted, incomplete or from different sessions." << endl;
        exit( EXIT_FAILURE );
    }
    
    charRanks[0] = charCounts[4];
    charRanks[1] = charRanks[0] + charCounts[0];
    charRanks[2] = charRanks[1] + charCounts[1];
    charRanks[3] = charRanks[2] + charCounts[2];
    
    // Load index
    sizePerIndex = 37;
    fread( &indexSize, 8, 1, idx );
    fread( &markSize, 8, 1, idx );
    index_ = new uint8_t[indexSize * sizePerIndex];
    marks_ = new ReadId[markSize];
    buff = new uint8_t[bwtPerIndex*2];
    fread( index_, 1, indexSize * sizePerIndex, idx );
    fread( marks_, 4, markSize, idx );
    
    runFlag = 1 << 7;
    runMask = ~runFlag;
    memset( &isBaseRun, false, 256 );
    isBaseRun[255] = true;
    for ( int i ( 0 ); i < 4; i++ )
    {
        isBaseRun[i * 63 + 62] = true;
        memset( &decodeBaseChar[ i * 63 ], i, 63 );
        decodeBaseChar[ 252 + i ] = 4;
        decodeBaseRun[ 252 + i ] = i + 1;
        for ( int j ( 0 ); j < 63; j++ )
        {
            decodeBaseRun[ i * 63 + j ] = j + 1;
        }
    }
    
    CharCount ranks;
    setRank( 0, 0, ranks );
    memcpy( &baseCounts[0][0], &ranks.counts, 32 );
    for ( int i ( 0 ); i < 4; i++ )
    {
        setRank( i, charCounts[i], ranks );
        memcpy( &baseCounts[i + 1][0], &ranks.counts, 32 );
        setRank( i, baseCounts[0][i], ranks );
        memcpy( &midRanks[i][0], &ranks.counts, 32 );
    }
    
    if ( mer )
    {
        uint8_t kmer = 12;
        kmerLen = kmer;
        uint64_t merSize = pow( 4, kmer ) * 16;
        mers = new uint8_t[merSize];
        fread( mers, 1, merSize, mer );
    }
    else mers = NULL;
}

IndexReader::~IndexReader()
{
    if ( buff ) delete[] buff;
    if ( index_ ) delete[] index_;
    if ( marks_ ) delete[] marks_;
}

void IndexReader::countRange( uint8_t i, CharId rank, CharId count, CharCount &ranks, CharCount &counts )
{
    if ( i > 3 )
    {
        ranks.clear();
        counts.clear();
        return;
    }
    setRank( i, rank, ranks );
    setRank( i, rank + count, counts );
//    counts -= ranks;
    for ( int j ( 0 ); j < 4; j++ )
    {
        counts.counts[j] -= ranks.counts[j];
    }
    counts.endCounts -= ranks.endCounts;
}

void IndexReader::countRange( uint8_t i, CharId rank, CharId edge, CharId count, CharCount &ranks, CharCount &edges, CharCount &counts )
{
    if ( !count )
    {
        ranks.clear();
        edges.clear();
        counts.clear();
        return;
    }
    setRank( i, rank, ranks );
    setRank( i, rank + edge, edges );
    setRank( i, rank + edge + count, counts );
    counts -= edges;
    edges -= ranks;
//    for ( int j ( 0 ); j < 4; j++ )
//    {
//        counts.counts[j] -= edges.counts[j];
//        edges.counts[j] -= ranks.counts[j];
//    }
//    counts.endCounts -= edges.endCounts;
//    edges.endCounts -= ranks.endCounts;
}

void IndexReader::createSeeds( string &fn, int mer )
{
    FILE* fp = fopen( fn.c_str(), "wb" );
    assert( mer == 12 );
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            CharId rank, edge, count;
            setBaseAll( i, j, rank, edge, count );
            createSeeds( fp, j, 2, mer, rank, edge, count );
        }
    }
    fclose( fp );
}

void IndexReader::createSeeds( FILE* fp, int i, int it, int limit, CharId rank, CharId edge, CharId count )
{
    if ( it >= limit )
    {
        ReadId outEdges = edge, outCount = count;
        fwrite( &rank, 8, 1, fp );
        fwrite( &outEdges, 4, 1, fp );
        fwrite( &outCount, 4, 1, fp );
        return;
    }
    
    CharCount ranks, edges, counts;
    countRange( i, rank, edge, count, ranks, edges, counts );
    
    for ( int j = 0; j < 4; j++ ) createSeeds( fp, j, it+1, limit, ranks[j], edges[j], counts[j] );
}

int IndexReader::primeOverlap( uint8_t* q, CharId &rank, CharId &count )
{
    if ( q[0] > 3 || q[1] > 3 )
    {
        rank = count = 0;
        return 0;
    }
    assert( !mers );
    int ol = mers ? 12 : 2;
    rank = midRanks[ q[0] ][ q[1] ];
    count = baseCounts[ q[0] + 1 ][ q[1] ] - rank;
    return ol;
}

void IndexReader::primeOverlap( string &seq, vector<uint8_t> &q, CharId &rank, CharId &count, int &ol, bool drxn )
{
    ol = mers && seq.size() >= 12 ? 12 : 2;
    if ( drxn ) for ( int i = 0; i++ < ol; ) q.push_back( charToInt[ seq.end()[-i] ] );
    else for ( int i = 0; i < ol; i++ ) q.push_back( charToIntComp[ seq[i] ] );
    
    if ( mers )
    {
        assert( false );
    }
    else
    {
        rank = midRanks[ q[0] ][ q[1] ];
        count = baseCounts[ q[0] + 1 ][ q[1] ] - rank;
    }
}

int IndexReader::setBaseAll( vector<uint8_t> &q, CharId &rank, CharId &count )
{
    rank = count = 0;
    if ( q.size() < 12 || !mers ) assert( false ); 
    CharId p = 0;
    for ( int i = 0; i < kmerLen; i++ )
    {
        assert( q[i] < 4 );
        p += pow( 4, kmerLen - i - 1 ) * q[i] * 16;
    }
    ReadId inEdge, inCount;
    memcpy( &rank, &mers[p], 8 );
    memcpy( &inEdge, &mers[p+8], 4 );
    memcpy( &inCount, &mers[p+12], 4 );
    count = inEdge + inCount;
    
    return kmerLen;
}

void IndexReader::setBaseAll( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    rank = baseCounts[i][j];
    count = baseCounts[ i + 1 ][j] - rank;
}

void IndexReader::setBaseAll( uint8_t i, uint8_t j, CharId &rank, CharId &edge, CharId &count )
{
    rank = baseCounts[i][j];
    edge = midRanks[i][j] - rank;
    count = baseCounts[ i + 1 ][j] - edge - rank;
}

ReadId IndexReader::setBaseMap( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    ReadId base = 0;
    for ( int ii = 0; ii < i; ii++ ) base += baseCounts[0][ii];
    for ( int jj = 0; jj < j; jj++ ) base += midRanks[i][jj] - baseCounts[i][jj];
    rank = baseCounts[i][j];
    count = midRanks[i][j] - rank;
    return base;
}

void IndexReader::setBaseOverlap( uint8_t i, uint8_t j, CharId &rank, CharId &count )
{
    rank = midRanks[i][j];
    count = baseCounts[ i + 1 ][j] - rank;
}

void IndexReader::setRank( uint8_t i, CharId rank, CharCount &ranks )
{
    rank += charRanks[i];
    CharId rankMark = rank / indexPerMark;
    CharId rankIndex = marks_[rankMark];
    
    CharId totalCount = setRankIndex( rankIndex, ranks );
    assert( totalCount <= rank );
    if ( rank - totalCount >= bwtPerIndex && rankIndex + 1 < indexSize )
    {
        CharCount tmpRanks;
        CharId tmpTotal = setRankIndex( rankIndex + 1, tmpRanks );
        while ( tmpTotal <= rank )
        {
            memcpy( &ranks.counts, &tmpRanks.counts, 32 );
            memcpy( &ranks.endCounts, &tmpRanks.endCounts, 4 );
            totalCount = tmpTotal;
            ++rankIndex;
            if ( rank - totalCount < bwtPerIndex || rankIndex + 1 == indexSize ) break;
            tmpTotal = setRankIndex( rankIndex + 1, tmpRanks );
        }
    }
    
    CharId rankBwt = rankIndex * bwtPerIndex;
    CharId rankChunk = bwtPerIndex;
    uint8_t offset = index_[(rankIndex * sizePerIndex)+36];
    rankBwt -= offset;
    rankChunk += offset;
    CharId rankLeft = rank - totalCount;
    if ( rankLeft < rankChunk ) rankChunk = rankLeft;
    
    fseek( bwt, rankBwt + beginBwt, SEEK_SET );
    fread( buff, 1, rankChunk, bwt );
    
    uint8_t c;
    CharId p = 0, thisRun, addRun;
    
    while ( rankLeft )
    {
        c = decodeBaseChar[ buff[p] ];
        thisRun = decodeBaseRun[ buff[p] ];
        if ( isBaseRun[ buff[p++] ] )
        {
            addRun = buff[p] & runMask;
            uint8_t byteCount = 0;
            while ( buff[p++] & runFlag )
            {
                addRun ^= ( buff[p] & runMask ) << ( 7 * ++byteCount );
            }
            thisRun += addRun;
        }
        
        if ( thisRun > rankLeft ) thisRun = rankLeft;
        if ( c == 4 )
        {
            ranks.endCounts += thisRun;
        }
        else
        {
            ranks.counts[c] += thisRun;
        }
        rankLeft -= thisRun;
    }
}

CharId IndexReader::setRankIndex( CharId rankIndex, CharCount &ranks )
{
    CharId indexBegin = rankIndex * sizePerIndex;
    memcpy( &ranks.counts, &index_[indexBegin], 32 );
    memcpy( &ranks.endCounts, &index_[indexBegin+32], 4 );
    return ( ranks[0] + ranks[1] + ranks[2] + ranks[3] + ranks.endCounts );
}

