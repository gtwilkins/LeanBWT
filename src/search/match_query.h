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

#ifndef MATCH_QUERY_H
#define MATCH_QUERY_H

#include "types.h"
#include "index_reader.h"
#include "query_binary.h"
#include "query_structs.h"
#include "shared_structs.h"

struct MatchRead
{
    MatchRead( ReadId id ):id( id ){};
    string seq;
    ReadId id;
    int coord[2];
};

class MatchQuery
{
    bool query( CharId rank, CharId count, uint8_t c, int i, int j, int len, int errLeft, int d );
    void match( int errors );
    
    IndexReader* ir_;
    vector<uint8_t> q_[2];
    vector<int> blocks_[2];
    vector<QueryHit> hits_[2];
    int len_;
    
public:
    MatchQuery( string seq, IndexReader* ir, int errors );
//    vector<MatchRead> yield( QueryBinaries* qb );
    vector<Read> yield( QueryBinaries* qb );
};

#endif /* MATCH_QUERY_H */

