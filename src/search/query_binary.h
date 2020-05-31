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

#ifndef QUERY_BINARY_H
#define QUERY_BINARY_H

#include "types.h"
#include "filenames.h"

class QueryBinaries
{
public:
    QueryBinaries( Filenames* fns );
    ~QueryBinaries(){};
    vector<ReadId> getIds( CharId ranks, CharId counts );
    string getSequence( ReadId id );
    
private:
    void decodeSequence( uint8_t* line, string &seq, uint8_t extLen, bool isRev, bool drxn );
    void set();
    
    FILE* bin_,* ids_;
    uint8_t binBegin_, idsBegin_, lineLen_;
    
    char decodeFwd[4][256];
    char decodeRev[4][256];
};

#endif /* QUERY_BINARY_H */

