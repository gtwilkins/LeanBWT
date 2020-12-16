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

#include "query_extension.h"
#include <algorithm>
#include <cassert>

extern Parameters params;

Ext::Ext( string& seq, ReadId id, int ol, bool drxn )
: ext_( drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol ) ), reads_{ ExtRead( id, seq.size()-ol, ol ) }
{}

Ext::Ext( Ext* base, int skim, bool drxn )
: ext_( base->ext_ ), exts_( base->exts_ ), reads_( base->reads_.begin()+skim, base->reads_.end() )
{
    int ext = base->reads_[skim-1].ext_;
    assert( skim < base->reads_.size() );
    base->reads_.erase( base->reads_.begin()+skim, base->reads_.end() );
    base->exts_.clear();
    base->ext_ = drxn ? base->ext_.substr( 0, ext ) : base->ext_.substr( base->ext_.size()-ext );
    for ( int i = 0; i < base->redundant_.size(); i++)
    {
        if ( reads_[0].ol_ <= base->redundant_[i].ol_ ) redundant_.push_back( base->redundant_[i] );
        if ( ext < base->redundant_[i].ext_ ) base->redundant_.erase( base->redundant_.begin()+i-- );
    }
    shift( ext, drxn );
}

Ext::~Ext()
{
    for ( Ext* e : exts_ ) delete e;
}


bool Ext::add( ReadId id, int ext, int ol, int ins )
{
    if ( ( ext <= reads_.back().ext_ ) && ( ( reads_.back().ext_+reads_.back().ol_ ) > ( ext+ol ) ) ) redundant_.push_back( ExtRead( id, ext, ol ) );
    else reads_.insert( reads_.begin()+ins, ExtRead( id, ext, ol ) );
    return true;
}

bool Ext::absorb( bool drxn )
{
    if ( exts_.size() != 1 ) return false;
    for ( ExtRead& er : exts_[0]->reads_ ) reads_.push_back( ExtRead( er.id_, er.ext_ + ext_.size(), er.ol_ - ext_.size() ) );
    for ( ExtRead& er : exts_[0]->redundant_ ) redundant_.push_back( ExtRead( er.id_, er.ext_ + ext_.size(), er.ol_ - ext_.size() ) );
    ext_ = drxn ? ext_ + exts_[0]->ext_ : exts_[0]->ext_ + ext_;
    exts_ = exts_[0]->exts_;
    return true;
}

int Ext::count( bool inclMp, int pairDrxn )
{
    count_ = 0;
    for ( Ext* e : exts_ ) count_ = max( count_, e->count( inclMp, pairDrxn ) );
    
    Lib* lib;
    for ( ExtRead& er : reads_ ) 
    {
        if ( ( lib = params.getLib( er.id_ ) ) ? ( inclMp || lib->isPe ) && ( pairDrxn == 2 || lib->getDrxn( er.id_ ) == pairDrxn )
                                               : ( pairDrxn == 2 && inclMp ) ) count_++;
    }
    
    return count_;
}

void Ext::cull()
{
    int maxReads = 0;
    for ( Ext* e : exts_ ) maxReads = max( maxReads, e->count_ );
    int cutoff = min( 4, maxReads / 100 + ( maxReads > 9 ) );
    for ( int i = 0; i < exts_.size(); i++ ) if ( exts_[i]->count_ <= cutoff )
    {
        delete exts_[i];
        exts_.erase( exts_.begin() + i-- );
    }
    for ( Ext* e : exts_ ) e->cull();
}

void Ext::cull( vector<Ext*>& exts, int minReads, bool setCounts, bool drxn )
{
    if ( setCounts ) for ( Ext* e : exts ) e->count( true, 2 );
    for ( int i = 0; i < exts.size(); i++ ) if ( exts[i]->count_ < minReads )
    {
        delete exts[i];
        exts.erase( exts.begin() + i-- );
    }
    for ( Ext* e : exts ) if ( !e->exts_.empty() ) cull( e->exts_, minReads, false, drxn );
    for ( Ext* e : exts ) e->absorb( drxn );
}


void Ext::sanitise( vector<Ext*>& exts, int minOl )
{
    for ( int i = 0; i < exts.size(); i++ ) if ( exts[i]->reads_[0].ol_ < minOl )
    {
        delete exts[i];
        exts.erase( exts.begin() + i-- );
    }
    for ( Ext* ext : exts ) Ext::sanitise( ext->exts_, minOl );
}

void Ext::set( string& seq, bool drxn )
{
    while ( exts_.size() == 1 )
    {
        Ext* ext = exts_[0];
        for ( ExtRead& er : ext->reads_ ) reads_.push_back( ExtRead( er.id_, er.ext_+ext_.size(), er.ol_-ext_.size() ) );
        for ( ExtRead& er : ext->redundant_ ) redundant_.push_back( ExtRead( er.id_, er.ext_+ext_.size(), er.ol_-ext_.size() ) );
        ext_ = drxn ? ext_ + ext->ext_ : ext->ext_ + ext_;
        exts_ = ext->exts_;
        ext->exts_.clear();
        delete ext;
    }
    seq_ = drxn ? seq.substr( seq.size()-reads_[0].ol_ ) + ext_ : ext_ + seq.substr( 0, reads_[0].ol_ );
    for ( Ext* e : exts_ ) e->set( seq_, drxn );
}

void Ext::shift( int ext, bool drxn )
{
    assert( ext < ext_.size() );
    ext_ = drxn ? ext_.substr( ext ) : ext_.substr( 0, ext_.size()-ext );
    for ( ExtRead& er : reads_ )
    {
        er.ext_ -= ext;
        er.ol_ += ext;
    }
    for ( ExtRead& er : redundant_ )
    {
        er.ext_ -= ext;
        er.ol_ += ext;
    }
}

Exts::Exts( string& base, int coord, bool drxn )
: seq_( drxn ? base.substr( 0, coord ) : base.substr( base.size()-coord ) ), coord_( drxn ? coord : base.size() - coord )
{
    
}

Exts::~Exts()
{
    for ( Ext* e : exts_ ) delete e;
}

bool Exts::add( vector<Ext*>& exts, string seq, ReadId id, int ol, bool drxn )
{
    int ext = seq.size()-ol;
    bool added = false;
    for ( Ext* e : exts )
    {
        int i = 0, j = 0, limit = min( (int)e->ext_.size(), ext );
        while ( i < limit && ( drxn ? seq[i+ol] == e->ext_[i] : seq.end()[-i-ol-1] == e->ext_.end()[-i-1] ) ) i++;
        
        // Perfect match
        if ( i+ol == seq.size() && ( added = true ) && e->add( id, ext, ol, e->reads_.size() ) ) continue;
        
        // Matched and surpassed the existing ext seq
        if ( i == e->ext_.size() )
        {
            // Continue to pre-existing branches
            if ( !e->exts_.empty() ) return add( e->exts_, seq, id, i+ol, drxn );
            
            // append unbranched extend
            e->ext_ = drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol );
            return e->add( id, ext, ol, e->reads_.size() );
        }
        
        while ( j < e->reads_.size() && e->reads_[j].ext_ <= i ) j++;
        
        // Fork after read
        if ( j && i+ol < seq.size() )
        {
            Ext* alt = new Ext( e, j, drxn );
            e->exts_ = { alt, new Ext( seq, id, ol, drxn ) };
            e->exts_.back()->shift( e->ext_.size(), drxn );
            assert( !e->exts_[0]->ext_.empty() && !e->exts_[1]->ext_.empty() );
            return true;
        }
    }
    
    // Create new branch
    if ( !added ) exts.push_back( new Ext( seq, id, ol, drxn ) );
    
    return true;
}

bool Exts::cull( int minOl, int minReads, int pairDrxn )
{
    Ext::sanitise( exts_, minOl );
    for ( Ext* e : exts_ ) e->count( true, pairDrxn );
    for ( int i = 0 ; i < exts_.size(); i++ ) if ( exts_[i]->count_ < minReads )
    {
        delete exts_[i];
        exts_.erase( exts_.begin() + i-- );
    }
    for ( Ext* e : exts_ ) e->cull();
    return exts_.empty();
}

