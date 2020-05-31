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

#include "query_structs.h"
#include <cassert>

QueryHit::QueryHit( ReadId rank, ReadId count, int coord )
: rank_( rank ), count_( count ), coord_( coord )
{}

QueryHit::QueryHit( ReadId rank, ReadId count, int coord, vector<QueryHit>& hits )
: QueryHit( rank, count, coord )
{
    int i = 0, hit = -1;
    if ( coord == 56 )
    {
        int x = 0;
    }
    for ( ; i < hits.size() && hits[i].coord_ <= coord; i++ ) if ( coord == hits[i].coord_ )
    {
        if ( hits[i].rank_+hits[i].count_ < rank || rank+count < hits[i].rank_ ) continue;
        if ( hits[i].rank_ <= rank && rank+count <= hits[i].rank_+hits[i].count_ ) return;
        if ( hit < 0 ) hit = i;
        ReadId limits[2]{ min( hits[i].rank_, rank ), max( hits[i].rank_+hits[i].count_, rank+count ) };
        hits[hit].rank_ = rank = limits[0];
        hits[hit].count_ = count = limits[1] - limits[0];
        if ( hit != i ) hits.erase( hits.begin() + i-- );
        
    }
    if ( hit < 0 ) hits.insert( hits.begin()+i, *this );
}

