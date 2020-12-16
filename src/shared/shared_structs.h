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

#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H

#include <vector>
#include "types.h"
#include "query_extension.h"

struct Coords
{
    Coords():Coords( 0, 0 ){};
    Coords( int32_t i, int32_t j ){ coords_[0] = i; coords_[1] = j; };
    int32_t& operator[]( int i ){ return coords_[i]; };
    void operator+=( int i ){ coords_[0] += i; coords_[1] += i; };
    static vector< pair<Coords, Coords> > align( string a, string b, int minAlign );
    static bool conflict( pair<Coords, Coords>& a, pair<Coords, Coords>& b );
    bool conflict( Coords& alt, bool drxn );
    int len(){ return coords_[1] - coords_[0]; };
    int32_t coords_[2];
};

struct Mapped
{
    Mapped( ReadId id, Coords coords ): id_( id ), coords_( coords ){};
    static unordered_set<ReadId> getIds( vector<Mapped>& reads, int coord, bool coordDrxn, bool readDrxn );
    static void sort( vector<Mapped>& reads, bool ascending, bool coordDrxn );
    ReadId id_;
    Coords coords_;
};

struct Read
{
    Read( string seq, ReadId id, int i, int j ): seq_( seq ), id_( id ), coords_( Coords( i, j ) ){};
    static unordered_set<ReadId> getIds( vector<Read>& reads, int coord, bool coordDrxn, bool readDrxn );
    static void sort( vector<Read>& reads, bool ascending, bool coordDrxn );
    string seq_;
    ReadId id_;
    Coords coords_;
};

struct Scythe
{
    Scythe( Coords base, Coords alt ): keep_( true ){ coords_[0] = base; coords_[1] = alt; };
    int bridge( Scythe& r );
    Coords coords_[2];
    unordered_set<ReadId> ids_[2];
    vector<int> cover_[2];
    int limit_[2];
    bool keep_;
};

struct Node
{
    Node( string seq ): seq_( seq ), coords_( 0, seq.size() ){ };
    Node( Read read );
    ~Node();
    void add( Read read, int off );
    void add( string seq, bool rezero, bool drxn );
    int countReads( int counted, int limit, bool drxn );
    int countReads( unordered_set<Node*>& pathed, bool doSort, bool drxn );
    static void cull( vector<Node*>& nodes );
    void get( unordered_set<Node*>& nodes, bool inclNode, bool drxn );
    void getConnected( unordered_set<Node*>& nodes );
    bool getDiploid( vector<Node*> paths[2] );
    vector<Node*> getHaploid( bool drxn );
    static vector<Node*> getNodes( vector<Read>& reads, bool coordinated, bool& failure );
    void getPathCounts( unordered_map<Node*, int>& counts, unordered_set<Node*>* block, bool drxn, bool force=false );
    static vector< pair< vector<Node*>, int > > getPaths( vector<Node*>& nodes );
    static string getSeq( vector<Node*>& path );
    bool isSubstantial( int minReads );
    bool isSubstantial( vector<int>& cover, unordered_set<Node*>& pathed, int ol, int minReads, bool drxn );
    bool isSubstantial( int reads, bool drxn );
    bool isWeak( int cutoff, bool drxn );
    void rebase();
    bool removeEdge( Node* node, bool drxn );
    int scythe( Scythe& scythe, Node* base, unordered_set<ReadId>& blocked, int i, bool drxn);
    static void simplify( vector<Node*>& nodes );
    static void sortEdges( vector<Node*>& nodes );
    static bool trim( unordered_set<Node*>& pathed, unordered_set<Node*>& safe, vector<Node*>& nodes, int cutoff );
    void trim( int len, bool drxn );
    void test( vector<Node*>& nodes, bool drxn );
    static void test( vector<Node*>& nodes );
    bool test( bool doTrim=false );
    bool unfork( bool drxn );
    string seq_;
    vector< pair<Node*, int> > edges_[2];
    vector<Mapped> reads_;
    Coords coords_;
};

struct Locus
{
    Locus( string a, string b );
    ~Locus();
    bool add( Node* node, bool extended[2] );
    bool confirm( float similarity, int readLen, int minLen, float maxRatio );
    void cull( vector<Node*>& nodes );
    bool extend( unordered_set<Node*>& pathed, vector<Node*>& nodes, int readLen, bool drxn );
    void fill( vector<Read> reads, vector<Read>& unmapped, int i, bool extended[2] );
    float getSimilarity();
    bool setCoverage( int readLen, float maxRatio );
    bool setLen( int minLen );
    bool setSimilarity( float minSimilarity );
    Node* node_[2];
    float similarity_, coverage_;
    int len_;
};


#endif /* SHARED_STRUCTS_H */

