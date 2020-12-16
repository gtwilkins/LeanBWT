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


#ifndef QUERY_OVERLAP_H
#define QUERY_OVERLAP_H

#include "index_reader.h"
#include "query_binary.h"
#include "query_structs.h"
#include "query_extension.h"

class QueryOverlap
{
    QueryOverlap( string seq, IndexReader* ir, int minOl, bool drxn );
    void query( CharId rank, CharId count, int i );
    IndexReader* ir_;
    vector<uint8_t> q_;
    vector<QueryHit> hits_;
    CharCount ranks_, counts_;
    int minOl_;
public:
    static ReadId countOverlaps( string seq, IndexReader* ir, int minOl, bool drxn );
    static Exts* getExtensions( string seq, IndexReader* ir, QueryBinaries* qb, int minOl, bool drxn );
};

#endif /* QUERY_OVERLAP_H */

