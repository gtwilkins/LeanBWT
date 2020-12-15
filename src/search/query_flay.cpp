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

#include "query_flay.h"
#include "index_reader.h"
#include "parameters.h"
#include "shared_functions.h"
#include <cassert>

extern Parameters params;

QueryFlay::QueryFlay( string seq, IndexReader* ir, int minFlay, bool drxn )
: ir_( ir )
{
    for ( int i = 0; i < 4; i++ ) ols_[i].resize( seq.size(), false );
    minFlay = max( min( (int)params.readLen, minFlay ), 30 );
    minLen_[0] = max( minFlay, 30 );
    minLen_[1] = max( minLen_[0] / 2, 30 );
    
    assert( minFlay < seq.size() );
    for ( int i = 0; i < seq.size(); i++ ) q_.push_back( drxn ? charToIntComp[ seq[i] ] : charToInt[ seq.end()[-i-1] ] );
    for ( int i = 0; i < seq.size()-minLen_[1]; i++ )
    {
        CharId rank, count;
        ReadId base = ir->setBaseMap( q_[i], q_[i+1], rank, count );
        query( rank, count, base, i+1, 2 );
    }
}

void QueryFlay::add( ReadId base, ReadId count, int i, int len, uint8_t c )
{
    if ( !ols_[c][i] && !( ols_[c][i] = ( len >= minLen_[0] ) ) ) return;
    int j = 0;
    for ( ; j < flays_.size() && ( flays_[j].coord_ < i || ( i == flays_[j].coord_ && flays_[j].c_<= c ) ); j++ );
    flays_.insert( flays_.begin()+j, Match( base, count, i, len, c ) );
}

void QueryFlay::query( CharId rank, CharId count, ReadId base, int i, int len )
{
    CharCount ranks, counts;
    ir_->countRange( q_[i++], rank, count, ranks, counts );
    
    
    for ( int j = 0; j < 4; j++ ) if ( counts[j] )
    {
        if ( q_[i] == j ) query( ranks[j], counts[j], base, i, len+1 );
        else if ( len >= minLen_[1] ) add( params.seqCount - base - counts[j], counts[j], i, len, j );
        base += counts[j];
    }
    
    if ( counts.endCounts ) matches_.push_back( Match( params.seqCount - base - counts.endCounts, counts.endCounts, i, len, 4 ) );
}

vector<Exts*> QueryFlay::flay( string seq, QueryBinaries* qb, bool drxn )
{
    vector<Exts*> exts;
    for ( int i = 0; i < flays_.size(); i++ )
    {
        if ( !i || flays_[i].coord_ != flays_[i-1].coord_ || flays_[i].c_ != flays_[i-1].c_ ) exts.push_back( new Exts( seq, flays_[i].coord_, drxn ) );
        
        vector<ReadId> ids;
        for ( ReadId id : qb->getIds( flays_[i].rank_, flays_[i].count_ ) ) ids.push_back( drxn ? id : params.getRevId( id ) );
        
        for ( int j = ids.size(); j-- > 0; ) exts.back()->add( exts.back()->exts_, qb->getSequence( ids[j] ), ids[j], flays_[i].len_, drxn );
    }
    for ( int i = 0; i < exts.size(); i++ ) if ( exts[i]->cull( params.readLen*2/3, 3, drxn ) )
    {
        delete exts[i];
        exts.erase( exts.begin() + i-- );
    }
    for ( Exts* es : exts ) for ( Ext* e : es->exts_ ) e->set( es->seq_, drxn );
    
    return exts;
}

void QueryFlay::test( string seq, QueryBinaries* qb )
{
    for ( Match& m : matches_ )
    {
        for ( ReadId id : qb->getIds( m.rank_, m.count_ ) )
        {
            string s = qb->getSequence( id );
            assert( seq.find( s ) != string::npos );
        }
    }
    
}

