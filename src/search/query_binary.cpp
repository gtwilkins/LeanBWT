/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
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

#include "query_binary.h"
#include "parameters.h"
#include <cassert>
#include <iostream>
#include <string.h>

extern Parameters params;

QueryBinaries::QueryBinaries( Filenames* fns )
{
    bin_ = fns->getBinary( true, false );
    ids_ = fns->getReadPointer( fns->ids, false );
    
    uint64_t binId, idsId;
    fread( &binBegin_, 1, 1, bin_ );
    fread( &binId, 8, 1, bin_ );
    fread( &idsBegin_, 1, 1, ids_ );
    fread( &idsId, 8, 1, ids_ );
    if ( binId != idsId )
    {
        cerr << endl << "Error: disagreement among input files." << endl;
        exit( EXIT_FAILURE );
    }
    
    uint8_t libCount, readLen, cycles;
    uint32_t coverage;
    fread( &readLen, 1, 1, bin_ );
    fread( &cycles, 1, 1, bin_ );
    params.readLen = readLen;
    lineLen_ = 1 + ( readLen + 3 ) / 4;
    if ( cycles != readLen + 1 )
    {
        cerr << endl << "Error: input data files appear either incomplete or corrupted." << endl;
        exit( EXIT_FAILURE );
    }
    fread( &params.isCalibrated, 1, 1, bin_ );
    fread( &coverage, 4, 1, bin_ );
    params.cover = (float)coverage / (float)100000;
    fread( &params.seqCount, 4, 1, bin_ );
    fread( &libCount, 1, 1, bin_ );
    ReadId countSoFar = 0;
    for ( int i ( 0 ); i < libCount; i++ )
    {
        Lib lib;
        uint16_t libMed, libMin, libMax;
        fread( &lib.count, 4, 1, bin_ );
        fread( &libMed, 2, 1, bin_ );
        fread( &libMin, 2, 1, bin_ );
        fread( &libMax, 2, 1, bin_ );
        fread( &lib.orientation, 1, 1, bin_ );
        fread( &lib.isPe, 1, 1, bin_ );
        lib.size = libMed;
        lib.minDist = libMin;
        lib.maxDist = libMax;
        lib.endCount = lib.count + countSoFar;
        countSoFar = lib.endCount;
        params.libs.push_back( lib );
    }
    
    params.set();
    set();
}


void QueryBinaries::decodeSequence( uint8_t* line, string &seq, uint8_t extLen, bool isRev, bool drxn )
{
    if ( isRev )
    {
        uint8_t i = drxn ? extLen : line[0];
        uint8_t j = drxn ? 0 : line[0] - extLen;
        while ( i -- > j )
        {
            seq += decodeRev[ i & 0x3 ][ line[ 1 + ( i / 4 ) ] ];
        }
    }
    else
    {
        uint8_t i = drxn ? line[0] - extLen : 0;
        uint8_t j = drxn ? line[0] : extLen;
        while ( i < j )
        {
            seq += decodeFwd[ i & 0x3 ][ line[ 1 + ( i / 4 ) ] ];
            ++i;
        }
    }
}


string QueryBinaries::getSequence( ReadId id )
{
    bool isRev = id & 0x1;
    string seq;
    uint8_t line[lineLen_];
    CharId seekId = CharId( id / 2 ) * lineLen_ + binBegin_;
    fseek( bin_, seekId, SEEK_SET );
    fread( &line, 1, lineLen_, bin_ );
    decodeSequence( line, seq, line[0], isRev, 1 );
    return seq;
}

vector<ReadId> QueryBinaries::getIds( CharId rank, CharId count )
{
    vector<ReadId> readIds( count );
    ReadId idBuff[count];
    CharId seekId = rank * 4 + idsBegin_;
    fseek( ids_, seekId, SEEK_SET );
    fread( &idBuff, 4, count, ids_ );
    for ( int i ( 0 ); i < count; i++ )
    {
        readIds[i] = idBuff[i];
    }
    return readIds;
}


void QueryBinaries::set()
{
    for ( int i ( 0 ); i < 256; i++ )
    {
        if ( i < 64 ){ decodeFwd[0][i] = 'A'; decodeRev[0][i] = 'T'; }
        else if ( i < 128 ){ decodeFwd[0][i] = 'C'; decodeRev[0][i] = 'G'; }
        else if ( i < 192 ){ decodeFwd[0][i] = 'G'; decodeRev[0][i] = 'C'; }
        else{ decodeFwd[0][i] = 'T'; decodeRev[0][i] = 'A'; }
        
        if ( ( i & 63 ) < 16 ){ decodeFwd[1][i] = 'A'; decodeRev[1][i] = 'T'; }
        else if ( ( i & 63 ) < 32 ){ decodeFwd[1][i] = 'C'; decodeRev[1][i] = 'G'; }
        else if ( ( i & 63 ) < 48 ){ decodeFwd[1][i] = 'G'; decodeRev[1][i] = 'C'; }
        else{ decodeFwd[1][i] = 'T'; decodeRev[1][i] = 'A'; }
        
        if ( ( i & 15 ) < 4 ){ decodeFwd[2][i] = 'A'; decodeRev[2][i] = 'T'; }
        else if ( ( i & 15 ) < 8 ){ decodeFwd[2][i] = 'C'; decodeRev[2][i] = 'G'; }
        else if ( ( i & 15 ) < 12 ){ decodeFwd[2][i] = 'G'; decodeRev[2][i] = 'C'; }
        else{ decodeFwd[2][i] = 'T'; decodeRev[2][i] = 'A'; }
        
        if ( ( i & 3 ) == 0 ){ decodeFwd[3][i] = 'A'; decodeRev[3][i] = 'T'; }
        else if ( ( i & 3 ) == 1 ){ decodeFwd[3][i] = 'C'; decodeRev[3][i] = 'G'; }
        else if ( ( i & 3 ) == 2 ){ decodeFwd[3][i] = 'G'; decodeRev[3][i] = 'C'; }
        else{ decodeFwd[3][i] = 'T'; decodeRev[3][i] = 'A'; }
    }
}
