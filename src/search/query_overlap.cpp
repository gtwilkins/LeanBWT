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

#include "query_overlap.h"
#include "constants.h"

QueryOverlap::QueryOverlap( string seq, IndexReader* ir, int minOl, bool drxn )
: ir_( ir ), minOl_( minOl )
{
    for ( int i = 0; i < seq.size(); i++ ) q_.push_back( drxn ? charToInt[ seq.end()[-i-1] ] : charToIntComp[ seq[i] ] );
    
    CharId rank, count;
    ir->setBaseOverlap( q_[0], q_[1], rank, count );
    query( rank, count, 1 );
}

ReadId QueryOverlap::countOverlaps( string seq, IndexReader* ir, int minOl, bool drxn )
{
    ReadId overlaps = 0;
    QueryOverlap qo( seq, ir, minOl, drxn );
    for ( QueryHit& qh : qo.hits_ ) overlaps += qh.count_;
    return overlaps;
}

void QueryOverlap::query( CharId rank, CharId count, int i )
{
    ir_->countRange( q_[i++], rank, count, ranks_, counts_ );
    for ( int j = 0; j < 4; j++ ) if ( counts_[j] && ( i >= q_.size() || q_[i] > 3 || q_[i] == j ) ) query( ranks_[j], counts_[j], i );
    if ( counts_.endCounts && i >= minOl_ ) hits_.push_back( QueryHit( ranks_.endCounts, counts_.endCounts, i ) );
}