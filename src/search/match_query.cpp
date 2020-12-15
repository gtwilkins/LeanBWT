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
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <iomanip>
#include <algorithm>

extern Parameters params;

MatchQuery::MatchQuery( string seq, IndexReader* ir, int errors )
: ir_( ir ), len_( seq.size() ), failure_( false )
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
    
    
    auto t_start = std::chrono::high_resolution_clock::now();
    for ( int d : { 0, 1 } ) for ( int i = 0; !failure_ && i < dBlocks[d]; i++ )
    {
        CharId rank, count;
        ir_->setBaseAll( q_[d][ blocks_[d][i] ], q_[d][ blocks_[d][i]+1 ], rank, count );
        
        bool seed = !i && min( errors, blockCount ) > 2;
        query( rank, count, q_[d][ blocks_[d][i]+1 ], blocks_[d][i]+1, i, 1, seed, d );
        if ( seed ) for ( int j : { 0, 1 } ) for ( int k = 0; k < 4; k++ ) if ( k != q_[d][j] )
        {
            ir_->setBaseAll( j ? q_[d][0] : k, j ? k : q_[d][1], rank, count );
            query( rank, count, j ? k : q_[d][1], blocks_[d][1], 1, 1, 0, d );
        }
        if ( seed ) i++;
        if ( double( ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC ) > 3 ) failure_ = true;
    }
}

bool MatchQuery::query( CharId rank, CharId count, uint8_t c, int i, int j, int len, int errLeft, int d )
{
    CharCount ranks, counts;
    ir_->countRange( c, rank, count, ranks, counts );
    i++;
    
    if ( ( ++len > min( len_, 50 ) ) && counts.endCounts ) QueryHit( ranks.endCounts, counts.endCounts, d ? len_-i : i, hits_[d] );
    
    if ( j < blocks_[d].size() && i >= blocks_[d][j+1] && ++j ) errLeft++;
    
    for ( int k = 0; k < 4; k++ ) if ( counts[k] )
    {
        if ( i >= len_ || k == q_[d][i] ) query( ranks[k], counts[k], k, i, j, len, errLeft, d );
        else if ( errLeft ) query( ranks[k], counts[k], k, i, j, len, errLeft-1, d );
    }
}

//vector<MatchRead> MatchQuery::yield( QueryBinaries* qb )
//{
//    vector<MatchRead> reads;
//    unordered_set<ReadId> used;
//    for ( int d : { 0, 1 } ) for ( QueryHit& qh : hits_[d] ) for ( ReadId id : qb->getIds( qh.rank_, qh.count_ ) )
//    {
//        if ( !used.insert( id = ( d ? id : params.getRevId( id ) ) ).second ) continue;
//        reads.push_back( MatchRead( id ) );
//        reads.back().seq = qb->getSequence( id );
//        reads.back().coord[0] = reads.back().coord[1] = qh.coord_;
//        reads.back().coord[d] += ( d ? reads.back().seq.size() : -reads.back().seq.size() );
//    }
//    return reads;
//}

vector<Read> MatchQuery::yield( QueryBinaries* qb )
{
    vector<Read> reads;
    unordered_set<ReadId> used;
    for ( int d : { 0, 1 } ) for ( QueryHit& qh : hits_[d] ) for ( ReadId id : qb->getIds( qh.rank_, qh.count_ ) )
    {
        if ( !used.insert( id = ( d ? id : params.getRevId( id ) ) ).second ) continue;
        reads.push_back( Read( qb ? qb->getSequence( id ) : "", id, qh.coord_, qh.coord_ ) );
        reads.back().coords_[d] += ( d ? reads.back().seq_.size() : -reads.back().seq_.size() );
    }
    for ( int i = 0; i > reads.size(); i++ ) if ( params.isReadMp( reads[i].id_ ) ) reads.erase( reads.begin() + i-- );
    return reads;
}

MatchedQuery::MatchedQuery( string header, string seq, IndexReader* ir, QueryBinaries* qb, int errors )
: header_( header ), seq_( seq )
{
    Coords coords[2];
    for ( Read r : MatchQuery( seq_, ir, errors ).yield( qb ) )
    {
        int it = seq.find( r.seq_ );
        if ( it != string:: npos ) exact_.push_back( Read( r.seq_, r.id_, it, it + r.coords_.len() ) );
        else if ( ( it = r.seq_.find( seq ) ) != string::npos ) exact_.push_back( Read( r.seq_, r.id_, -it, r.coords_.len()-it ) );
        else
        {
            int best = 0;
            for ( pair<Coords, Coords> align : Coords::align( seq_, r.seq_, 10 ) ) if ( align.first[1]-align.first[0] > best )
            {
                coords[0] = align.first;
                coords[1] = align.second;
                best = align.first[1]-align.first[0];
            }
            bool cut = true;
            if ( best && ( !coords[1][0] || !coords[0][0] ) && coords[0][1] < seq.size() )
            {
                string q = "CTG";
                for ( int i = 0; i < min( 3, (int)r.seq_.size() - coords[1][1] ); i++ ) if ( r.seq_[i+coords[1][1]] != q[i] ) cut = false;
                if ( cut )
                {
                    r.seq_ = r.seq_.substr( 0, coords[1][1] );
                    exact_.push_back( Read( r.seq_, r.id_, coords[0][0], coords[0][1] ) );
                    continue;
                }
            }
            if ( best && ( coords[1][1] == r.seq_.size() || coords[0][1] == seq.size() ) && coords[0][0] )
            {
                string q = "GAC";
                for ( int i = 0; i < min( 3, coords[1][0] ); i++ ) if ( r.seq_[coords[1][0]-1-i] != q[i] ) cut = false;
                if ( cut )
                {
                    r.seq_ = r.seq_.substr( coords[1][0] );
                    exact_.push_back( Read( r.seq_, r.id_, coords[0][0], coords[0][1] ) );
                    continue;
                }
            }
            if ( best )
            {
                if ( !coords[0][0] && coords[1][1] == r.seq_.size() ) exact_.push_back( Read( r.seq_, r.id_, coords[0][1]-r.coords_.len(), coords[0][1] ) );
                else if ( !coords[1][0] && coords[0][1] == seq.size() ) exact_.push_back( Read( r.seq_, r.id_, coords[0][0], coords[0][0]+r.coords_.len() ) );
                else inexact_.push_back( MatchRead( r.id_, r.seq_, coords[0], coords[1] ) );
            }
            else unmatched_.push_back( r );
        }
    }
    Read::sort( exact_, true, 0 );
    sort( inexact_.begin(), inexact_.end(), []( MatchRead& a, MatchRead& b ){ return a.query_[0]-a.read_[0] < b.query_[0]-b.read_[0]; } );
    Read::sort( unmatched_, true, 0 );
}

void MatchedQuery::compete( vector<MatchedQuery>& queries )
{
    unordered_set<ReadId> exact;
    unordered_map<ReadId, int> best;
    for ( MatchedQuery& mq : queries ) for ( Read& r : mq.exact_ ) exact.insert( r.id_ );
    for ( MatchedQuery& mq : queries ) for ( int i = 0; i < mq.inexact_.size(); i++ )
    {
        if ( exact.find( mq.inexact_[i].id_ ) == exact.end() )
        {
            auto ins = best.insert( make_pair( mq.inexact_[i].id_, mq.inexact_[i].read_.len() ) );
            if ( !ins.second && ins.first->second < mq.inexact_[i].read_.len() ) ins.first->second = mq.inexact_[i].read_.len();
        }
        else mq.inexact_.erase( mq.inexact_.begin() + i-- );
    }
    for ( MatchedQuery& mq : queries ) for ( int i = 0; i < mq.inexact_.size(); i++ )
    {
        auto it = best.find( mq.inexact_[i].id_ );
        if ( it != best.end() && it->second > mq.inexact_[i].read_.len() ) mq.inexact_.erase( mq.inexact_.begin() + i-- );
    }
    for ( MatchedQuery& mq : queries ) for ( int i = 0; i < mq.unmatched_.size(); i++ )
    {
        if ( best.find( mq.inexact_[i].id_ ) != best.end() || exact.find( mq.inexact_[i].id_ ) != exact.end() ) mq.unmatched_.erase( mq.unmatched_.begin() + i-- );
    }
}
