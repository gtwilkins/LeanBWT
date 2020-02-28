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

#include "index_writer.h"
#include <cassert>
#include <string.h>
#include <algorithm>
#include <iostream>
#include "timer.h"
#include "index_reader.h"

IndexWriter::IndexWriter( PreprocessFiles* fns, ReadId indexChunk, ReadId markChunk )
: bwtPerIndex( indexChunk ), countsPerMark( markChunk )
{
    fns->setIndexWrite( bwt, idx );
    fread( &bwtBegin, 1, 1, bwt );
    fread( &id, 8, 1, bwt );
    fread( &bwtSize, 8, 1, bwt );
    fread( &charCounts[4], 8, 1, bwt );
    fread( &charCounts, 8, 4, bwt );
    
    contFlag = 1 << 7;
    contMask = ~contFlag;
    indexSize = 1 + bwtSize / bwtPerIndex;
    markSize = 1 + ( charCounts[0] + charCounts[1] + charCounts[2] + charCounts[3] + charCounts[4] ) / countsPerMark;
    
    buff = new uint8_t[IDX_BUFFER];
    marks = new ReadId[ markSize ];
    
    currByte = 0;
    memset( &counts, 0, 40 );
    memset( &decodeBaseChar[0], 0, 63 );
    memset( &decodeBaseChar[63], 1, 63 );
    memset( &decodeBaseChar[126], 2, 63 );
    memset( &decodeBaseChar[189], 3, 63 );
    memset( &decodeBaseChar[252], 4, 4 );
    memset( &maxBaseRun, 63, 4 );
    maxBaseRun[4] = 4;
    for ( int i ( 0 ); i < 4; i++ )
    {
        decodeBaseRun[ 252 + i ] = i + 1;
        for ( int j ( 0 ); j < 63; j++ )
        {
            decodeBaseRun[ i * 63 + j ] = j + 1;
        }
    }
    
    writeIndex();
    fclose( bwt );
    fclose( idx );
//    writeMers( fns );
}

IndexWriter::~IndexWriter()
{
    if ( buff ) delete[] buff;
    if ( marks ) delete[] marks;
}

void IndexWriter::writeIndex()
{
    double indexStartTime = clock();
    
    uint8_t indexBegin = 69;
    fwrite( &indexBegin, 1, 1, idx );
    fwrite( &id, 8, 1, idx );
    fwrite( &bwtPerIndex, 4, 1, idx );
    fwrite( &countsPerMark, 4, 1, idx );
    fwrite( &indexSize, 8, 1, idx );
    fwrite( &markSize, 8, 1, idx );
    
    uint8_t currChar;
    uint8_t currRunBytes = 0;
    ReadId currChunkBytes = 0;
    CharId currRank = 0;
    CharId currRun, currAddRun;
    CharId indexCount = 1;
    CharId markCount = 0;
    bool startByte = true;
    ReadId p = IDX_BUFFER - 1;
    ReadId endCount = counts[4];
    
    fwrite( &counts, 8, 4, idx );
    fwrite( &endCount, 4, 1, idx );
    fwrite( &currRunBytes, 1, 1, idx );         // Dummy offset
    
    CharId bwtLeft = bwtSize + 1;
    while ( --bwtLeft )
    {
        if ( ++p == IDX_BUFFER )
        {
            fread( buff, 1, min( bwtLeft, IDX_BUFFER ), bwt );
            p = 0;
        }
        
        if ( startByte )
        {
            currChar = decodeBaseChar[ buff[p] ];
            currRun = decodeBaseRun[ buff[p] ];
            currRunBytes = 0;
            currAddRun = 0;
            startByte = currRun != maxBaseRun[currChar];
        }
        else
        {
            currAddRun ^= ( ( buff[p] & contMask ) << ( 7 * currRunBytes++ ) );
            startByte = !( buff[p] & contFlag );
        }
        
        if ( startByte )
        {
            if ( currChunkBytes >= bwtPerIndex )
            {
                currRunBytes -= ( currChunkBytes - bwtPerIndex );
                endCount = counts[4];
                fwrite( &counts, 8, 4, idx );
                fwrite( &endCount, 4, 1, idx );
                fwrite( &currRunBytes, 1, 1, idx );
                currChunkBytes -= bwtPerIndex;
                
                ReadId currMarks = ( currRank - 1 ) / countsPerMark;
                while ( markCount <= currMarks )
                {
                    marks[ markCount++ ] = indexCount - 1;
                }
                
                ++indexCount;
            }
            currRun += currAddRun;
            counts[currChar] += currRun;
            currRank += currRun;
        }
        
        ++currChunkBytes;
    }
    
    while ( markCount < markSize )
    {
        marks[ markCount++ ] = indexCount - 1;
    }
    
    if ( indexCount != indexSize )
    {
        cerr << endl << "Unexpected error constructing index." << endl;
        exit( EXIT_FAILURE );
    }
    
    fwrite( marks, 4, markCount, idx );
    
    cout << endl << "Indexing transformed data... completed!" << endl;
    cout << "Summary:" << endl;
    cout << "Indexed " << to_string( bwtSize ) << " BWT entries" << endl;
    cout << "Comprising " << to_string( counts[0] + counts[1] + counts[2] + counts[3] ) 
            << " sequence characters from " << to_string( counts[4] / 2 ) << " sequence reads " << endl;
    cout << "Created " << to_string( indexCount ) << " index points" << endl;
    cout << "Time taken: " << getDuration( indexStartTime ) << endl;
}

void IndexWriter::writeMers( PreprocessFiles* fns )
{
    IndexReader ir( fns );
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            
        }
    }
}
