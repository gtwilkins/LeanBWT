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

#include "match_query.h"
#include "constants.h"
#include "query_binary.h"
#include "parameters.h"
#include <cassert>

extern Parameters params;

MatchQuery::MatchQuery( string seq, IndexReader* ir, int errors )
: ir_( ir ), len_( seq.size() )
{
    q_[0].resize( seq.size(), 0 );
    q_[1].resize( seq.size(), 0 );
    for ( int i = 0; i < seq.size(); i++ )
    {
        q_[0][i] = 3 - charToInt[ seq[i] ];
        q_[1][i] = charToInt[ seq.end()[-i-1] ];
    }
    match( errors );
}

bool MatchQuery::query( CharId rank, CharId count, uint8_t c, int i, int j, int len, int errLeft, int d )
{
    CharCount ranks, counts;
    ir_->countRange( c, rank, count, ranks, counts );
    i++;
    
    if ( ( ++len > min( len_, 50 ) ) && counts.endCounts ) QueryHit( ranks.endCounts, counts.endCounts, d ? len_-i : i, hits_[d] );
    
    if ( j < blocks_[d].size() && i >= blocks_[d][j+1] && ++j ) errLeft++;
    
    if ( i < len_ ) for ( int k = 0; k < 4; k++ ) if ( counts[k] )
    {
        if ( k == q_[d][i] ) query( ranks[k], counts[k], k, i, j, len, errLeft, d );
        else if ( errLeft ) query( ranks[k], counts[k], k, i, j, len, errLeft-1, d );
    }
}

void MatchQuery::match( int errors )
{
    int blockSize = min( min( len_, 100 ), max ( 6, 100 / ( errors+1 ) ) );
    int blockCount = len_ / blockSize;
    int blockExtra = len_ - ( blockSize * blockCount );
    int extraBegin = max( 0, ( ( blockCount+1 ) / 2 ) - ( ( blockExtra+1 ) / 2 ) );
    int dBlocks[2]{ ( 1 + blockCount ) / 2, max( 1, blockCount / 2 ) };
    assert( blockExtra >= 0 );
    
    blocks_[0].push_back( 0 );
    for ( int i = 0; i < blockCount; i++ ) blocks_[0].push_back( blocks_[0].back() + blockSize + ( i > extraBegin && blockExtra && blockExtra-- ) );
    for ( int i = 0; i < blocks_[0].size(); i++ ) blocks_[1].push_back( len_ - blocks_[0].end()[-i-1] );
    
    for ( int d : { 0, 1 } ) for ( int i = 0; i < dBlocks[d]; i++ )
    {
        CharId rank, count;
        ir_->setBaseAll( q_[d][ blocks_[d][i] ], q_[d][ blocks_[d][i+1] ], rank, count );
        
        bool seed = !i && min( errors, blockCount ) > 2;
        query( rank, count, q_[d][ blocks_[d][i+1] ], blocks_[d][i]+1, i, 1, seed, d );
        if ( seed ) for ( int j : { 0, 1 } ) for ( int k = 0; k < 4; k++ ) if ( k != q_[d][j] )
        {
            ir_->setBaseAll( j ? q_[d][0] : k, j ? k : q_[d][1], rank, count );
            query( rank, count, j ? k : q_[d][1], blocks_[d][1], 1, 1, 0, d );
        }
        if ( seed ) i++;
    }
}

vector<MatchRead> MatchQuery::yield( QueryBinaries* qb )
{
    vector<MatchRead> reads;
    unordered_set<ReadId> used;
    for ( int d : { 0, 1 } ) for ( QueryHit& qh : hits_[d] ) for ( ReadId id : qb->getIds( qh.rank_, qh.count_ ) )
    {
        if ( !used.insert( id = ( d ? id : params.getRevId( id ) ) ).second ) continue;
        reads.push_back( MatchRead( id ) );
        reads.back().seq = qb->getSequence( id );
        reads.back().coord[0] = reads.back().coord[1] = qh.coord_;
        reads.back().coord[d] += ( d ? reads.back().seq.size() : -reads.back().seq.size() );
    }
    return reads;
}
