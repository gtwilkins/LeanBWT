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

#ifndef TRANSFORM_FUNCTIONS_H
#define TRANSFORM_FUNCTIONS_H

#include <algorithm>
#include <cassert>
#include "transform_structs.h"
#include "transform_constants.h"

inline void iter( CharId &p, uint8_t &i )
{
    if ( i == 3 )
    {
        p++;
        i = 0;
    }
    else
    {
        i++;
    }
}

//inline void rebuffBin( FILE* &fBin, uint8_t* binBuff, CharId bufferSize, CharId &fileLeft, CharId &p )
//{
//    bufferSize = min( bufferSize, fileLeft );
//    fread( binBuff, 1, bufferSize, fBin );
//    fileLeft -= bufferSize;
//    p = 0;
//}

//inline void rebuffBwt( FILE* &fBwt, uint8_t* bwtBuff, CharId fileLeft, CharId &p )
//{
//    fread( bwtBuff, 1, min( fileLeft, BWT_BUFFER ), fBwt );
//    p = 0;
//}

//inline void rebuffIds( FILE* &fIds, ReadId* idsBuff, CharNum &fileLeft, CharNum &p )
//{
//    fread( idsBuff, 4, min( fileLeft, IDS_BUFFER ), fIds );
//    p = 0;
//}
//
//inline void rebuffIds( FILE* &fIds, FILE* &fPos, ReadId* idsBuff, CharNum* posBuff, CharNum fileLeft, CharNum &p )
//{
//    fread( idsBuff, 4, min( fileLeft, IDS_BUFFER ), fIds );
//    fread( posBuff, 8, min( fileLeft, IDS_BUFFER ), fPos );
//    p = 0;
//}

inline void writeBwtBuff( FILE* &fBwt, uint8_t* bwtBuff, CharId &p )
{
    fwrite( bwtBuff, 1, p, fBwt );
    p = 0;
}

inline void writeIdsBuff( FILE* &fIds, ReadId* idsBuff, CharId &p )
{
    fwrite( idsBuff, 1, p, fIds );
    p = 0;
}

inline void writePosBuff( FILE* &fIds, FILE* &fPos, ReadId* idsBuff, CharId* posBuff, CharId &p )
{
    fwrite( idsBuff, 4, p, fIds );
    fwrite( posBuff, 8, p, fPos );
    p = 0;
}

inline uint8_t read2Bit( uint8_t* buff, CharId &p , uint8_t &i )
{
    if ( !i )
    {
        i++;
        return buff[p] >> 6;
    }
    else if ( i == 1 )
    {
        i++;
        return ( buff[p] >> 4 ) & 3;
    }
    else if ( i == 2 )
    {
        i++;
        return ( buff[p] >> 2 ) & 3;
    }
    else
    {
        i = 0;
        return buff[p++] & 3;
    }
}

inline void write2Bit( uint8_t* buff, CharId &p, uint8_t &i, uint8_t c )
{
    if ( !i )
    {
        buff[p] = c << 6;
        i++;
    }
    else if ( i == 1 )
    {
        buff[p] += c << 4;
        i++;
    }
    else if ( i == 2 )
    {
        buff[p] += c << 2;
        i++;
    }
    else
    {
        buff[p] += c;
        p++;
        i = 0;
    }
}

//inline void writeEncoded( uint8_t* buff, CharId &p, uint8_t c )
//{
//    if ( incrementEncoded[ buff[p] ] == c )
//    {
//        buff[p]++;
//    }
//    else
//    {
//        buff[++p] = baseEncode[c];
//    }
//}

//inline void readPos( CharId in, CharId &pos, uint8_t &c )
//{
//    c = in & 0x3;
//    pos = in >> 2;
//}

//inline uint8_t readChar( uint8_t* chars, ReadId id )
//{
//    return byteToInt[ id & 3 ][ chars[ id / 4 ] ];
//}

//inline ReadId readId( uint8_t* chars, CharId* posBuff, ReadId* idsBuff, CharId &p, CharId &pos, uint8_t* c )
//{
//    ReadId id = idsBuff[p];
//    readPos( posBuff[p], pos, c[0] );
//    c[1] = readChar( chars, id );
//    p++;
//    return id;
//}

//inline ReadId readId( CharId* posBuff, ReadId* idsBuff, CharId &p, CharId &pos, uint8_t &c )
//{
//    ReadId id = idsBuff[p];
//    readPos( posBuff[p], pos, c );
//    p++;
//    return id;
//}

#endif /* TRANSFORM_FUNCTIONS_H */

