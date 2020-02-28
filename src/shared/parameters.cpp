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

#include "parameters.h"
#include "filenames.h"
#include <cassert>
#include <iostream>

bool Lib::doAddMarker( bool thisRev, int &pairRev, bool drxn )
{
    if ( orientation == 0 ) // UK
    {
        return true;
    }
    else if ( orientation == 1 ) // FF
    {
        pairRev = thisRev;
    }
    else if ( orientation == 2 ) // FR
    {
        if ( thisRev == drxn )
        {
            pairRev = !thisRev;
            return true;
        }
    }
    else if ( orientation == 3 ) // RF
    {
        if ( thisRev != drxn )
        {
            pairRev = !thisRev;
            return true;
        }
    }
    return false;
}

int Lib::getPair( SeqNum &id )
{
    uint8_t readVer = id & 0x3;
    id = id - readVer;
    // Unknown
    if ( orientation == 0 )
    {
        id += !readVer;
        return 2;
    }
    
    // FF
    if ( orientation == 1 ) id += ( readVer > 1 ? 0 : 2 ) + ( readVer & 0x1 ? 1 : 0 ); // 0 <-> 2; 1 <-> 3
    // FR or RF
    else id += ( readVer > 1 ? 0 : 2 ) + ( readVer & 0x1 ? 0 : 1 ); // 0 <-> 3; 1 <-> 2
    
    return orientation < 3 ? !( readVer & 0x1 ) : readVer & 0x1;
}

bool Lib::getPair( SeqNum &id, int &pairDrxn )
{
    if ( id >= endCount ) return false;
    pairDrxn = getPair( id );
    return true;
}

void Lib::setMinMax()
{
    minDist = max( ( ( size / 5 ) - 100 ), 1 );
    maxDist = ( size * 1.2 ) + 200;
};

Parameters::Parameters()
: rna( false )
{
    isSet = haploid = isCalibrated = false;
    drxns[0] = drxns[1] = true;
}

void Parameters::checkReady()
{
    if ( !isCalibrated )
    {
        cerr << "Error: sequence data has not been calibrated. This can be done with the <calibrate> command." << endl;
        exit( EXIT_FAILURE );
    }
    bool anyPe = false;
    for ( Lib &lib : libs )
    {
        anyPe = anyPe || lib.isPe;
    }
    if ( !anyPe )
    {
        cerr << "Error: no paired end libraries detected in data files." << endl;
        exit( EXIT_FAILURE );
    }
    if ( !isSet )
    {
        cerr << "Unexpected error loading sequence data files." << endl;
        exit( EXIT_FAILURE );
    }
}

bool Parameters::isReadPe( ReadId id )
{
    Lib* lib = getLib( id );
    return lib && lib->isPe;
}

 bool Parameters::isReadMp( ReadId id )
 {
    Lib* lib = getLib( id );
    return lib && !lib->isPe;
 }

int32_t Parameters::getFurthestMpDist( int32_t coord, bool drxn )
{
    return coord + ( drxn ? maxMpMax : -maxMpMax );
}

int32_t Parameters::getFurthestPeDist( int32_t coord, bool drxn )
{
    return coord + ( drxn ? maxPeMax : -maxPeMax );
}

int32_t Parameters::getFurthestMpMean( int32_t coord, bool drxn )
{
    return coord + ( drxn ? maxMpMean : -maxMpMean );
}

int32_t Parameters::getFurthestPeMean( int32_t coord, bool drxn )
{
    return coord + ( drxn ? maxPeMean : -maxPeMean );
}

Lib* Parameters::getLib( const SeqNum &readId )
{
    for ( Lib &lib : libs )
    {
        if ( readId < lib.endCount )
        {
            if ( lib.orientation == 0 )
            {
//                SeqNum pairId = readId % 4;
//                readId = readId - pairId + !pairId;
                return NULL;
            }
            return &lib;
        }
    }
    return NULL;
}

int32_t Parameters::getLibSize( ReadId readId )
{
    for ( Lib &lib : libs ) if ( readId < lib.endCount ) return lib.size;
    return 0;
}

ReadId Parameters::getPairId( const SeqNum readId )
{
    Lib* lib = getLib( readId );
    if ( lib )
    {
        SeqNum readVer = readId & 0x3;
        SeqNum pairId = readId - readVer;
        if ( lib->orientation == 0 )
        {
            pairId += !readVer;
        }
        else if ( lib->orientation == 1 )
        {
            pairId += ( readVer > 1 ? 2 : 0 ) + ( readVer %2 ? 0 : 1 );
        }
        else
        {
            pairId += ( readVer > 1 ? 0 : 2 ) + ( readVer %2 ? 0 : 1 );
        }
        return pairId;
    }
    return readId;
}

ReadId Parameters::getRevId( const ReadId readId )
{
    return readId - ( readId & 0x1 ) + !( readId & 0x1 );
}

void Parameters::set()
{
    CharId totalPe = 0, totalPeMean = 0;
    maxPeMax = maxPeMean = maxMpMax = maxMpMean = 0;
    branchMinHits = 1;
    for ( Lib &lib : libs )
    {
        lib.setMinMax();
        maxMpMean = max( maxMpMean, lib.size );
        maxMpMax = max( maxMpMax, lib.maxDist );
        if ( lib.isPe )
        {
            maxPeMean = max( maxPeMean, lib.size );
            maxPeMax = max( maxPeMax, lib.maxDist );
            totalPe += lib.count;
            totalPeMean += (CharId)lib.size * (CharId)lib.count;
            branchMinHits += ( (float)lib.count / (float)seqCount ) * cover
                    * float( max( 0, lib.size - readLen - 10 ) ) / float( readLen * 6 );
        }
    }
    
    if ( totalPe ) avgPeMean = totalPeMean / totalPe;
    if ( isCalibrated )
    {
        peCover = ( (float)totalPe / (float)seqCount ) * cover;
        peRatio = peCover / cover;
        queryMpCutoff = 5 + ( peRatio > 0.9 ? 1 : ( log( 5000.0 ) / log( 1.0 / ( 1.0 - peRatio ) ) ) );
        readSpacing = (float)readLen / cover;
        isSet = true;
    }
}

void Parameters::setLimits( int32_t &limit )
{
    limit = max( 1000, min( 200000, limit ) );
    locusLimits[0] = -limit;
    locusLimits[1] = limit;
}

int Parameters::setPair( ReadId &id )
{
    Lib* lib = getLib( id );
    if ( !lib ) return 2;
    return lib->getPair( id );
}

bool Parameters::setPairId( ReadId &id, bool pairDrxn )
{
    if ( !isReadPe( id ) || pairDrxn == bool(id & 0x1) ) return false;
    id = getPairId( id );
    return true;
}

int32_t Parameters::shortLen()
{
    return max( readLen, maxPeMean - readLen );
}

void Parameters::write( FILE* fp )
{
    
}
