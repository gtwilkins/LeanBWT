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

#ifndef QUERY_FLAY_H
#define QUERY_FLAY_H

#include "types.h"
#include "constants.h"
#include "query_binary.h"
#include "index_reader.h"
#include "query_extension.h"

class QueryFlay
{
    struct Match
    {
        Match( ReadId rank, ReadId count, int coord, int len, uint8_t c ): rank_( rank ), count_( count ), coord_( coord ), len_( len ), c_( c ){};
        ReadId rank_, count_;
        int coord_, len_;
        uint8_t c_;
    };
    
    void add( ReadId base, ReadId count, int i, int len, uint8_t c );
    void query( CharId rank, CharId count, ReadId base, int i, int len );
    
    IndexReader* ir_;
    vector<uint8_t> q_;
    vector<bool> ols_[4];
    vector<Match> matches_, flays_;
    int minLen_[2];
    
public:
    QueryFlay( string seq, IndexReader* ir, int minFlay, bool drxn );
    vector< pair<ReadId, pair<int, int> > > fill( string seq, QueryBinaries* qb );
    vector<Exts*> flay( string seq, QueryBinaries* qb, bool drxn );
    void test( string seq, QueryBinaries* qb );
};


#endif /* QUERY_FLAY_H */

