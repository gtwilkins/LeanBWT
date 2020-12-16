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

#ifndef QUERY_EXTENSION_H
#define QUERY_EXTENSION_H

#include "types.h"
#include "parameters.h"

struct ExtRead
{
    ExtRead( ReadId id, int ext, int ol ): id_( id ), ext_( ext ), ol_( ol ){};
    ReadId id_;
    int ext_, ol_;
};

struct Ext
{
    Ext( string& seq, ReadId id, int ol, bool drxn );
    Ext( Ext* base, int skim, bool drxn );
    ~Ext();
    static void add( vector<Ext>& exts, string seq, ReadId id, int ol, bool drxn );
    bool add( ReadId id, int ext, int ol, int ins );
    bool absorb( bool drxn );
    int count( bool inclMp, int pairDrxn );
    void cull();
    static void cull( vector<Ext*>& exts, int minReads, bool setCounts, bool drxn );
    static void sanitise( vector<Ext*>& exts, int minOl );
    void set( string& seq, bool drxn );
    void shift( int ext, bool drxn );
    
    vector<Ext*> exts_;
    vector<ExtRead> reads_, redundant_;
    string seq_, ext_, ol_;
    int count_;
};

struct Exts
{
    Exts( string& base, int coord, bool drxn );
    ~Exts();
    bool add( vector<Ext*>& exts, string seq, ReadId id, int ol, bool drxn );
    bool cull( int minOl, int minReads, int pairDrxn );
    vector<Ext*> exts_;
    string seq_;
    int coord_;
};


#endif /* QUERY_EXTENSION_H */

