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
#include "overlap_query.h"
#include "constants.h"

//OverlapQuery::OverlapQuery( string seq, IndexReader* ir, bool drxn, int len )
//: ir_( ir )
//{
//    vector<uint8_t> q;
//    for ( int i = 0; i < seq.size(); i++ ) q.push_back( drxn ? charToInt[ seq.end()[-i-1] ] : charToIntComp[ seq[i] ] );
//    
//    CharId rank, count;
//    ir->setBaseOverlap( q[0], q[1], rank, count );
//    query( rank, count, 1 );
//}
//
//void OverlapQuery::query( CharId rank, CharId count, int i )
//{
//    CharCount ranks, counts;
//    ir_->countRange( q_[i++], rank, count, ranks, counts );
//    for ( int j = 0; j < 4; j++ ) if ( i >= q_.size() || q_[i] > 3 || q_[i] == j ) query( ranks[j], counts[j], i );
//    if ( counts.endCounts && i >= cutoff_ ) hits_.push_back( QueryHit( ranks.endCounts, counts.endCounts, i ) );
//}