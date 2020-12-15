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

#include "shared_structs.h"
#include "shared_functions.h"
#include "local_alignment.h"
#include "parameters.h"
#include "query_extension.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <iomanip>

vector< pair<Coords, Coords> > Coords::align( string a, string b, int minAlign )
{
    int len = minAlign / 2;
    vector< pair<Coords, Coords> > aligns;
    if ( len < a.size() ) for ( int i = 0; i < a.size(); )
    {
        int j = min( i, (int)a.size()-len );
        string q = a.substr( j, len );
        size_t it = b.find( q );
        while ( it != string::npos )
        {
            Coords x[2]{ Coords( j, j+len ), Coords( it, it+len ) };
            while ( x[0][0] && x[1][0] && a[ x[0][0]-1 ] == b[ x[1][0]-1 ] ){ x[0][0]--; x[1][0]--; };
            while ( x[0][1] < a.size() && x[1][1] < b.size() && a[ x[0][1] ] == b[ x[1][1] ] ){ x[0][1]++; x[1][1]++; };
            bool good = x[0][1]-x[0][0] >= minAlign;
            for ( pair<Coords, Coords>& align : aligns ) if ( align.first[0] == x[0][0] && align.first[1] == x[0][1] && align.second[0] == x[1][0] ) good = false;
            if ( good ) aligns.push_back( make_pair( x[0], x[1] ) );
            it = b.find( q, it+1 );
        }
        i += len;
    }
    
    return aligns;
}

bool Coords::conflict( pair<Coords, Coords>& a, pair<Coords, Coords>& b )
{
    return b.first[0] < a.first[1] || b.second[0] < a.second[1];
}

unordered_set<ReadId> Mapped::getIds( vector<Mapped>& reads, int coord, bool coordDrxn, bool readDrxn )
{
    unordered_set<ReadId> ids;
    for ( Mapped& read : reads ) if ( readDrxn ? coord < read.coords_[coordDrxn] : read.coords_[coordDrxn] < coord ) ids.insert( read.id_ );
    return ids;
}

void Mapped::sort( vector<Mapped>& reads, bool ascending, bool coordDrxn )
{
    if ( ascending ) std::sort( reads.begin(), reads.end(), [&]( Mapped& a, Mapped& b ){ return a.coords_[coordDrxn] < b.coords_[coordDrxn]; } );
    else std::sort( reads.begin(), reads.end(), [&]( Mapped& a, Mapped& b ){ return a.coords_[coordDrxn] > b.coords_[coordDrxn]; } );
}

unordered_set<ReadId> Read::getIds( vector<Read>& reads, int coord, bool coordDrxn, bool readDrxn )
{
    unordered_set<ReadId> ids;
    for ( Read& read : reads ) if ( readDrxn ? coord < read.coords_[coordDrxn] : read.coords_[coordDrxn] < coord ) ids.insert( read.id_ );
    return ids;
}

void Read::sort( vector<Read>& reads, bool ascending, bool coordDrxn )
{
    if ( ascending ) std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] < b.coords_[coordDrxn]; } );
    else std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] > b.coords_[coordDrxn]; } );
}

int Scythe::bridge( Scythe& r )
{
    int ids = 0, cover = 1 + max( cover_[0][ limit_[0]-1 ], r.cover_[0][ r.limit_[0] ] );
    for ( ReadId id : ids_[0] ) if ( r.ids_[0].find( id ) != r.ids_[0].end() ) ids++;
    for ( int i = limit_[0]; i < r.limit_[0]; i++ ) cover = max( cover, cover_[0][i] );
    return min( ids, cover );
}

Node::Node( Read read )
: seq_( read.seq_ ), reads_{ Mapped( read.id_, read.coords_ ) }, coords_( read.coords_ )
{
    
}

Node::~Node()
{
    for ( int d : { 0, 1 } ) for ( pair<Node*, int>& edge : edges_[d] ) for ( int i = 0; i < edge.first->edges_[!d].size(); i++ )
    {
        if ( edge.first->edges_[!d][i].first == this ) edge.first->edges_[!d].erase( edge.first->edges_[!d].begin() + i-- );
    }
}

void Node::add( Read read, int off )
{
    reads_.push_back( Mapped( read.id_, Coords( coords_[0] + off, coords_[0] + read.seq_.size() + off ) ) );
}

void Node::add( string seq, bool rezero, bool drxn )
{
    seq_ = drxn ? seq_ + seq : seq + seq_;
    coords_[drxn] += drxn ? seq.size() : -seq.size();
    if ( rezero && coords_[0] ) for ( Mapped& read : reads_ ) read.coords_ += -coords_[0];
    if ( rezero && coords_[0] ) coords_ += -coords_[0];
}

int Node::countReads( int counted, int limit, bool drxn )
{
    counted += reads_.size();
    if ( counted >= limit ) return limit;
    int maxCounted = counted;
    for ( pair<Node*, int>& edge : edges_[drxn] ) maxCounted = max( maxCounted, edge.first->countReads( counted, limit, drxn ) );
    return maxCounted;
}

int Node::countReads( unordered_set<Node*>& pathed, bool doSort, bool drxn )
{
    if ( !pathed.insert( this ).second ) return 0;
    int maxCount = 0;
    vector< pair< pair<Node*, int>, int> > counts;
    bool unsorted = false;
    for ( pair<Node*, int>& edge : edges_[drxn] )
    {
        int thisCount = edge.first->countReads( pathed, doSort, drxn );
        if ( !counts.empty() && counts.back().second < thisCount ) unsorted = true;
        maxCount = max( maxCount, thisCount );
        if ( doSort ) counts.push_back( make_pair( edge, thisCount ) );
    }
    if ( doSort && unsorted )
    {
        sort( counts.begin(), counts.end(), []( pair< pair<Node*, int>, int>& a, pair< pair<Node*, int>, int>& b){ return a.second > b.second;} );
        edges_[drxn].clear();
        for ( pair< pair<Node*, int>, int>& edge : counts ) edges_[drxn].push_back( edge.first );
    }
    pathed.erase( this );
    return maxCount + reads_.size();
}

void Node::cull( vector<Node*>& nodes )
{
    Node::simplify( nodes );
    Node::sortEdges( nodes );
    unordered_set<Node*> safe;
    for ( pair<vector<Node*>, int > path : Node::getPaths( nodes ) ) if ( find( nodes.begin(), nodes.end(), path.first[0] ) != nodes.end() )
    {
        if ( !path.first[0]->edges_[0].empty() ) continue;
        unordered_set<Node*> pathed;
        int reads = 0;
        for ( Node* node = path.first[0]; node && pathed.find( node ) == pathed.end(); )
        {
            pathed.insert( node );
            reads += node->reads_.size();
            node = node->edges_[1].empty() ? NULL : node->edges_[1][0].first;
        }
        safe.insert( path.first.begin(), path.first.end() );
        if ( path.second > 1 ) Node::trim( pathed, safe, nodes, 1 + ( reads > 9 ) + ( reads > 29 ) + ( reads / 100 ) );
    }
    Node::simplify( nodes );
}

void Node::get( unordered_set<Node*>& nodes, bool inclNode, bool drxn )
{
    if ( inclNode && !nodes.insert( this ).second ) return;
    for ( pair<Node*, int>& edge : edges_[drxn] ) edge.first->get( nodes, true, drxn );
}

void Node::getConnected( unordered_set<Node*>& nodes )
{
    nodes.insert( this );
    for ( int d : { 0, 1 } ) for ( pair<Node*, int >& edge : edges_[d] ) if ( nodes.insert( edge.first ).second ) edge.first->getConnected( nodes );
}

bool Node::getDiploid( vector<Node*> paths[2] )
{
    for ( int d : { 0, 1 } ) paths[d].clear();
    for ( int d : { 0, 1 } ) if ( edges_[d].size() > 2 ) return false;
    for ( int d : { 0, 1 } ) if ( paths[0].empty() && !edges_[d].empty() )
    {
        paths[0].push_back( edges_[d][0].first );
        paths[1].push_back( edges_[d][edges_[d].size() > 1].first );
    }
    if ( paths[0].empty() ) paths[0] = paths[1] = { this };
    
    unordered_set<Node*> used( paths[0].begin(), paths[0].end() );
    used.insert( paths[1].begin(), paths[1].end() );
    
    if ( !paths[0].empty() ) for ( int d : { 0, 1 } ) for ( bool good = true; good; )
    {
        Node* node[2]{ ( d ? paths[0].back() : paths[0][0] ) , ( d ? paths[1].back() : paths[1][0] ) };
        Node* edge[2]{ NULL, NULL };
        if ( node[0]->edges_[d].empty() || node[1]->edges_[d].empty() || node[0]->edges_[d].size() > 2 || node[1]->edges_[d].size() > 2 ) break;
        if ( node[0] == node[1] && node[0]->edges_[d].size() != 2 ) break;
        
        vector<Node*> edges[2];
        for ( int i : { 0, 1 } ) for ( pair<Node*, int>& e : node[i]->edges_[d] ) edges[i].push_back( e.first );
        
        for ( int i : { 0, 1 } ) if ( edges[i].size() <= edges[!i].size() ) for ( Node* e : edges[i] )
        {
            if ( find( edges[!i].begin(), edges[!i].end(), e ) == edges[!i].end() ) good = false;
        }
        for ( int i : { 0, 1 } ) if ( edges[i].size() == 1 ) edge[i] = edges[i][0];
        for ( int i : { 0, 1 } ) if ( edges[i].size() == 2 ) edge[i] = edges[i][ edges[i][0] == edge[!i] ];
        
        if ( !used.insert( edge[0] ).second || ( edge[0] != edge[1] && !used.insert( edge[1] ).second ) )
        {
            //loop?!
            return false;
        }
        
        if ( good ) for ( int i : { 0, 1 } ) paths[i].insert( d ? paths[i].end() : paths[i].begin(), edge[i] );
    }
    
    return used.find( this ) != used.end();
}

vector<Node*> Node::getHaploid( bool drxn )
{
    vector<Node*> path = { this };
    for ( Node* node = this; node = ( node->edges_[drxn].size() == 1 ? node->edges_[drxn][0].first: NULL ); )
    {
        if ( find( path.begin(), path.end(), node ) != path.end() ) break;
        path.insert( drxn ? path.end() : path.begin(), node );
    }
    return path;
}

vector<Node*> Node::getNodes( vector<Read>& reads, bool coordinated, bool& failure )
{
    auto t_start = std::chrono::high_resolution_clock::now();
    vector<Node*> nodes;
    sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.seq_.size() > b.seq_.size(); } );
    unordered_map<ReadId, int> matches;
    for ( Read& read : reads ) matches.insert( make_pair( read.id_, read.coords_[0] ) );
    for ( int i = 0; i < reads.size(); i++ )
    {
        Node* node = new Node( reads[i] );
        size_t it;
        for ( int j = i+1; j < reads.size(); j++ ) if ( ( it = node->seq_.find( reads[j].seq_ ) ) != string::npos )
        {
            node->add( reads[j], it );
            reads.erase( reads.begin() + j-- );
        }
        nodes.push_back( node );
        if ( double( ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC ) > 3 )
        {
            failure = true;
            for ( Node* node : nodes ) delete node;
            nodes.clear();
            return nodes;
        }
    }
    sort( nodes.begin(), nodes.end(), [&]( Node* a, Node* b ){ return a->coords_[0] < b->coords_[0]; } );

    int jj = 0;
    for ( int i = 0; i < nodes.size(); i++ )
    {
        for ( ; jj < nodes.size() && nodes[jj]->coords_[1] < nodes[i]->coords_[0]; jj++ );
        vector< pair<Node*, int> > ols;
        for ( int j = jj; j < nodes.size(); j++ ) if ( i != j )
        {
            int minOl = min( nodes[i]->seq_.size(), nodes[j]->seq_.size() ) * 7 / 10;
            if ( coordinated ) minOl = min( minOl, max( nodes[i]->coords_[1] - nodes[j]->coords_[0] - 50, 40 ) );
            int ol = 0;
            if ( ol = mapSeqOverlap( nodes[i]->seq_, nodes[j]->seq_, minOl ) ) ols.push_back( make_pair( nodes[j], ol ) );
        }
        sort( ols.begin(), ols.end(), []( pair<Node*, int>& a , pair<Node*, int>& b ){ return a.second > b.second; } );
        for ( int j = 0; j < ols.size(); j++ )
        {
            nodes[i]->edges_[1].push_back( ols[j] );
            ols[j].first->edges_[0].push_back( make_pair( nodes[i], ols[j].second ) );

            for ( int k = j+1; k < ols.size(); k++ )
            {
                int len = 0, limit = min( ols[j].first->seq_.size() - ols[j].second, ols[k].first->seq_.size() - ols[k].second );
                while ( len < limit && ols[j].first->seq_[ ols[j].second+len ] == ols[k].first->seq_[ ols[k].second+len ] ) len++;
                if ( len == limit ) ols.erase( ols.begin() + k-- );
            }
        }
    }
    return nodes;
}

void Node::getPathCounts( unordered_map<Node*, int>& counts, unordered_set<Node*>* block, bool drxn, bool force )
{
    int maxCount = 0;
    for ( pair<Node*, int>& edge : edges_[drxn] ) if ( !block || block->find( edge.first ) == block->end() )
    {
        auto it = counts.find( edge.first );
        if ( it == counts.end() && force ) continue;
        if ( it == counts.end() ) return;
        maxCount = max( maxCount, it->second );
    }
    if ( !counts.insert( make_pair( this, maxCount+reads_.size() ) ).second ) return;
    for ( pair<Node*, int>& edge : edges_[!drxn] ) if ( !block || block->find( edge.first ) == block->end() ) edge.first->getPathCounts( counts, block, drxn );
}
vector< pair< vector<Node*>, int > > Node::getPaths( vector<Node*>& nodes )
{
    vector< pair< vector<Node*>, int > > paths;
    for ( Node* node : nodes ) if ( node->edges_[0].empty() || node->edges_[0][0].first->edges_[1][0].first != node )
    {
        vector<Node*> path;
        int reads = 0;
        while ( node && find( path.begin(), path.end(), node ) == path.end() )
        {
            reads += node->reads_.size();
            path.push_back( node );
            node = node->edges_[1].empty() || node->edges_[1][0].first->edges_[0][0].first != node ? NULL : node->edges_[1][0].first;
        }
        paths.push_back( make_pair( path, reads ) );
    }
    sort( paths.begin(), paths.end(), []( pair< vector<Node*>, int >& a, pair< vector<Node*>, int >& b){ return a.second > b.second; } );
    return paths;
}

string Node::getSeq( vector<Node*>& path )
{
    if ( path.empty() ) return "";
    string seq = path[0]->seq_;
    for ( int i = 1; i < path.size(); i++ )
    {
        int ol = 0;
        for ( pair<Node*, int> e : path[i]->edges_[0] ) if ( e.first == path[i-1] ) ol = e.second;
        seq += path[i]->seq_.substr( ol );
    }
    return seq;
}

bool Node::isSubstantial( int minReads )
{
    vector<int> cover( seq_.size(), 0 );
    int diff = -coords_[0];
    for ( Mapped& read : reads_ ) for ( int i = read.coords_[0]; i < read.coords_[1]; i++ ) cover[i+diff]++;
    unordered_set<Node*> pathed;
    for ( int d : { 0, 1 } ) for ( pair<Node*, int> edge : edges_[d] ) if ( edge.first->isSubstantial( cover, pathed, edge.second, minReads, d ) );
    for ( int i : cover ) if ( i >= minReads ) return true;
    return false;
}

bool Node::isSubstantial( vector<int>& cover, unordered_set<Node*>& pathed, int ol, int minReads, bool drxn )
{
    if ( ol <= 0 || !pathed.insert( this ).second ) return false;
    for ( Mapped& read : reads_ ) for ( int i = abs( coords_[!drxn]-read.coords_[!drxn] ); i < ol; i++ ) cover[ drxn ? i : cover.size()-i-1 ]++;
    ol -= seq_.size();
    for ( int i : cover ) if ( i >= minReads ) return true;
    for ( pair<Node*, int> edge : edges_[drxn] ) if ( edge.first->isSubstantial( cover, pathed, ol+edge.second, minReads, drxn ) ) return true; 
    return false;
}

bool Node::isSubstantial( int reads, bool drxn )
{
    reads -= reads_.size();
    if ( reads <= 0 ) return true;
    for ( pair<Node*, int> edge : edges_[drxn] ) if ( edge.second > 40 && edge.first->isSubstantial( reads, drxn ) ) return true;
    return false;
}

bool Node::isWeak( int cutoff, bool drxn )
{
    int strong = 0, weak = seq_.size();
    Mapped::sort( reads_, !drxn, !drxn );
    for ( int i = 1; i < min( 3, (int)reads_.size() ); i++ ) weak = min( weak, abs( reads_[i].coords_[drxn] - reads_[i-1].coords_[!drxn] ) );
    for ( pair<Node*, int > edge : edges_[drxn] ) if ( edge.first->reads_.size() > 1 || !edge.first->edges_[drxn].empty() )
    {
        int ol = edge.second;
        if ( edge.first->reads_.size() == 2 )
        {
            ol = min( abs( edge.first->reads_[0].coords_[0] - edge.first->reads_[1].coords_[1] )
                    , abs( edge.first->reads_[0].coords_[1] - edge.first->reads_[1].coords_[0] ) );
        }
        else if ( edge.first->reads_.size() == 1 )
        {
            ol = 0;
            for ( pair<Node*, int > e : edge.first->edges_[drxn] ) ol = max( ol, e.second );
        }
        strong = max( strong, min( ol, edge.second ) );
    }
    return min( weak, strong ) < cutoff;
}

void Node::rebase()
{
    if ( !coords_[0] ) return;
    for ( Mapped& read : reads_ ) read.coords_ += -coords_[0];
    coords_ += -coords_[0];
    if ( !reads_.empty() ) test();
}

bool Node::removeEdge( Node* node, bool drxn )
{
    for ( int i = 0; i < edges_[drxn].size(); i++ ) if ( edges_[drxn][i].first == node )
    {
        edges_[drxn].erase( edges_[drxn].begin() + i-- );
        return true;
    }
    return false;
}

int Node::scythe( Scythe& scythe, Node* base, unordered_set<ReadId>& blocked, int i, bool drxn )
{
    int diff = -coords_[0], coord = scythe.coords_[i][drxn] + coords_[0], weak = 0, limit = drxn ? seq_.size() : 0, unscythed = 0;
    scythe.limit_[i] = drxn ? 0 : seq_.size();
    scythe.cover_[i] = vector<int>( seq_.size(), 0 );
    
    for ( Mapped& read : reads_ )
    {
        if ( drxn ? read.coords_[1] <= coord : coord <= read.coords_[0] ) unscythed++;
        else if ( i || blocked.find( read.id_ ) == blocked.end() )
        {
            limit = drxn ? min( limit, read.coords_[1] ) : max( limit, read.coords_[0] );
            for ( int j = read.coords_[0]; j < read.coords_[1]; j++ ) scythe.cover_[i][j+diff]++;
            scythe.ids_[i].insert( read.id_ );
        }
    }
    for ( Mapped& read : reads_ ) if ( drxn ? read.coords_[1] < limit : limit < read.coords_[0] )
    {
        scythe.limit_[i] = drxn ? max( scythe.limit_[i], read.coords_[1]+diff ) : min( scythe.limit_[i], read.coords_[0]+diff );
    }
    weak = 1 + scythe.cover_[i][ scythe.limit_[i] + ( drxn ? -1 : 0 ) ];
    if ( drxn ) for ( int j = scythe.limit_[i]; j < scythe.cover_[i].size(); j++ ) weak = max( weak, scythe.cover_[i][j] );
    else for ( int j = scythe.limit_[i]; j-- > 0; ) weak = max( weak, scythe.cover_[i][j] );
    weak = min( weak, (int)scythe.ids_[i].size() );
    
    if ( weak && base && i && unscythed > 2 )
    {
        int cut = 0;
        coord = scythe.coords_[0][drxn] + base->coords_[0];
        for ( Mapped& read : base->reads_ ) if ( drxn ? coord < read.coords_[1] : read.coords_[0] < coord ) cut++;
        if ( cut < min( 3, weak ) )
        {
            weak = cut;
            scythe.keep_ = false;
        }
    }
    
    return weak;
}

void Node::simplify( vector<Node*>& nodes )
{
    for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->edges_[1].size() == 1 && nodes[i]->edges_[1][0].first->edges_[0].size() == 1 )
    {
        Node* node[2]{ nodes[i], nodes[i]->edges_[1][0].first };
        int ol = nodes[i]->edges_[1][0].second;
        int off = node[1]->coords_[0] - node[0]->coords_[1] + ol;
        if ( off ) for ( Mapped& read : node[0]->reads_ ) read.coords_ += off;
        node[1]->coords_[0] = node[0]->coords_[0] + off;
        node[1]->seq_ = node[0]->seq_ + node[1]->seq_.substr( ol );
        node[1]->reads_.insert( node[1]->reads_.begin(), node[0]->reads_.begin(), node[0]->reads_.end() );
        node[1]->edges_[0] = node[0]->edges_[0];
        for ( pair<Node*, int>& edge : node[1]->edges_[0] ) for ( pair<Node*, int>& re : edge.first->edges_[1] ) if ( re.first == node[0] ) re.first = node[1];
//        for ( pair<Node*, int>& edge : node[1]->edges_[0] ) edge.first->edges_[1].push_back( make_pair( node[1], edge.second ) );
        delete nodes[i];
        nodes.erase( nodes.begin() + i-- );
    }
    for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->edges_[0].empty() && nodes[i]->edges_[1].empty() && nodes[i]->reads_.size() <= 2 )
    {
        delete nodes[i];
        nodes.erase( nodes.begin() + i-- );
    }
}

void Node::sortEdges( vector<Node*>& nodes )
{
    unordered_map<Node*, int> counts[2];
    for ( Node* node : nodes ) for ( int d : { 0, 1 } ) if ( node->edges_[d].empty() ) node->getPathCounts( counts[d], NULL, d );
    for ( Node* node : nodes ) for ( int d : { 0, 1 } ) if ( node->edges_[d].size() > 1 ) 
    {
        bool bad = false;
        vector< pair< pair<Node*, int>, int> > edges;
        for ( pair<Node*, int>& edge : node->edges_[d] ) if ( counts[d].find( edge.first ) == counts[d].end() ) bad = true;
        if ( bad )
        {
            unordered_set<Node*> pathed;
            node->countReads( pathed, true, d );
        }
        else
        {
            for ( pair<Node*, int>& edge : node->edges_[d] ) edges.push_back( make_pair( edge, counts[d][edge.first] ) );
            sort( edges.begin(), edges.end(), []( pair< pair<Node*, int>, int>& a, pair< pair<Node*, int>, int>& b ){ return a.second > b.second;} );
            node->edges_[d].clear();
            for ( pair< pair<Node*, int>, int>& edge : edges ) node->edges_[d].push_back( edge.first );
        }
    }
}

bool Node::trim( unordered_set<Node*>& pathed, unordered_set<Node*>& safe, vector<Node*>& nodes, int cutoff )
{
    unordered_set<Node*> connected;
    for ( Node* node : pathed ) node->getConnected( connected );
    
    unordered_map<Node*, int> counts[2];
    for ( Node* node : connected ) if ( pathed.find( node ) == pathed.end() ) for ( int d : { 0, 1 } )
    {
        bool unedged = true;
        for ( pair<Node*, int> edge : node->edges_[d] ) if ( pathed.find( edge.first ) == pathed.end() ) unedged = false;
        if ( unedged ) node->getPathCounts( counts[d], &pathed, d );
    }
    for ( int again = min( counts[0].size(), counts[1].size() ) < connected.size() - pathed.size(); again-- > 0; )
    {
        vector<Node*> forks[2], orphans[2];
        for ( int d : { 0, 1 } ) for ( Node* node : connected ) if ( pathed.find( node ) == pathed.end() && counts[d].find( node ) == counts[d].end() )
        {
            int ol = 0;
            vector< pair<Node*, int> > edges;
            for ( pair<Node*, int> edge : node->edges_[d] )
            {
                if ( pathed.find( edge.first ) != pathed.end() || counts[d].find( edge.first ) != counts[d].end() ) ol = max( ol, edge.second );
                else edges.push_back( edge );
            }
            if ( !ol && !edges.empty() ) orphans[d].push_back( node );
            sort( edges.begin(), edges.end(), []( pair<Node*, int>& a, pair<Node*, int>& b ){ return a.second > b.second; } );
            if ( ol ) for ( int i = 0; i < edges.size(); i++ )
            {
                if ( edges[i].second <= ol-20 )
                {
                    node->removeEdge( edges[i].first, d );
                    edges[i].first->removeEdge( node, !d );
                    edges.erase( edges.begin() + i-- );
                }
                else ol = max( ol, edges[i].second );
            }
            if ( edges.empty() && ( again = 1 ) ) node->getPathCounts( counts[d], &pathed, d );
            else if ( ol ) forks[d].push_back( node );
        }
        if ( again )
        {
            connected.clear();
            for ( Node* node : pathed ) node->getConnected( connected );
            for ( int d : { 0, 1 } ) for ( auto it = counts[d].begin(); it != counts[d].end(); )
            {
                if ( connected.find( it->first ) == connected.end() ) it = counts[d].erase( it );
                else it++;
            }
        }
        else for ( int d : { 0, 1 } ) for ( Node* node : forks[d] ) if ( counts[d].find( node ) == counts[d].end() && ( again = 1 ) )
        {
            node->getPathCounts( counts[d], &pathed, d, true );
        }
        
        if ( !again ) for ( int d : { 0, 1 } ) for ( Node* node : orphans[d] ) node->getPathCounts( counts[d], &pathed, d, true );
        else again = min( counts[0].size(), counts[1].size() ) < connected.size() - pathed.size();
        
    }
    assert( counts[0].size() == counts[1].size() );
    for ( pair<Node*, int > q : counts[0] ) if ( safe.find( q.first ) == safe.end() )
    {
        assert( counts[1].find( q.first ) != counts[1].end() );
        int reads = -q.first->reads_.size() + q.second + counts[1][q.first];
        if ( reads <= cutoff )
        {
            delete q.first;
            nodes.erase( remove( nodes.begin(), nodes.end(), q.first ), nodes.end() );
        }
    }
}

void Node::trim( int len, bool drxn )
{
    assert( len >= 0 );
    if ( len < 1 ) return;
    seq_.erase( drxn ? seq_.size()-len : 0, len );
    coords_[drxn] += drxn ? -len : len;
    for ( int i = 0; i < reads_.size(); i++ ) if ( reads_[i].coords_[0] < coords_[0] || coords_[1] < reads_[i].coords_[1] ) reads_.erase( reads_.begin() + i--);
}

void Node::test( vector<Node*>& nodes, bool drxn )
{
    if ( find( nodes.begin(), nodes.end(), this ) != nodes.end() ) assert( false );
    nodes.push_back( this );
    for ( pair<Node*, int>& edge : edges_[drxn] ) edge.first->test( nodes, drxn );
    nodes.pop_back();
}

void Node::test( vector<Node*>& nodes )
{
    for ( int d : { 0, 1 } ) for ( Node* node : nodes ) for ( pair<Node*, int> edge : node->edges_[d] )
    {
        bool good = false;
        for ( pair<Node*, int> re : edge.first->edges_[!d] ) if ( re.first == node ) good = true;
        assert( good );
    }
}

bool Node::test( bool doTrim )
{
    int coords[2]{ coords_[1], coords_[0] };
    for ( Mapped& read : reads_ )
    {
        coords[0] = min( coords[0], read.coords_[0] );
        coords[1] = max( coords[1], read.coords_[1] );
    }
    assert( coords_[0] <= coords[0] && coords[1] <= coords_[1] );
    if ( doTrim )
    {
        trim( coords_[1] - coords[1], 1 );
        trim( coords[0] - coords_[0], 0 );
    }
    else if ( !( coords[0] == coords_[0] && coords[1] == coords_[1] ) )
    {
        return false;
    }
    return true;
}

bool Node::unfork( bool drxn )
{
    if ( edges_[drxn].size() < 2 ) return false;
    vector< pair<Node*, int> > branches = edges_[drxn];
    sort ( branches.begin(), branches.end(), []( pair<Node*, int>& a, pair<Node*, int>& b ){ return a.second > b.second; } );
    
    bool unforked = false;
    for ( int i = 0; i < branches.size(); i++ ) if ( branches[i].first->isSubstantial( 3, drxn ) )
    {
        int diff[2];
        for ( int j = i+1; j < branches.size(); j++ ) if ( ( diff[0] = branches[i].second - branches[j].second ) >= 10 )
        {
            bool bad = false;
            for ( pair<Node*, int> edge : branches[j].first->edges_[!drxn] ) if ( ( edge.first != this ) && edge.first->isSubstantial( 3, !drxn ) )
            {
                if ( ( diff[1] = edge.second - branches[j].second ) < 10 ) continue;
            }
            if ( bad ) removeEdge( branches[j].first, drxn );
            if ( bad ) branches[j].first->removeEdge( this, !drxn );
        }
        break;
    }
    
    return unforked;
}

Locus::Locus( string a, string b )
: len_( min( a.size(), b.size()) )
{
    node_[0] = new Node( a );
    node_[1] = new Node( b );
}

Locus::~Locus()
{
    for ( int i : { 0, } ) delete node_[i];
}

bool Locus::add( Node* node, bool extended[2] )
{
    int ii = -1, jj = -1, dd, limit = 4, cut;
    vector<Scythe> scythes[2][2];
    unordered_set<ReadId> baseIds[2];
    for ( int i : { 0, 1 } ) for ( Mapped& read : node_[i]->reads_ ) baseIds[i].insert( read.id_ );
    for ( int i : { 0, 1 } ) for ( pair<Coords, Coords>& aligned : Coords::align( node_[i]->seq_, node->seq_, 60 ) ) for ( int d : { 0, 1 } )
    {
        if ( d ? !aligned.second[0] : aligned.second[1] == node->seq_.size() ) continue;
        Scythe scythe( aligned.first, aligned.second );
        if ( node->scythe( scythe, node_[i], baseIds[!i], 1, d ) >= min( 2, (int)node->reads_.size() ) ) continue;
        scythes[i][d].push_back( scythe );
        if ( ( cut = node_[i]->scythe( scythes[i][d].back(), NULL, baseIds[!i], 0, !d ) ) > limit ) continue;
        ii = i;
        jj = scythes[i][d].size()-1;
        dd = d;
        limit = cut;
    }
    
    Scythe* paired[2]{ NULL, NULL };
    for ( int i : { 0, 1 } ) for ( Scythe& l : scythes[i][0] ) for ( Scythe& r : scythes[i][1] ) if ( l.coords_[0][0] < r.coords_[0][0] )
    {
        if ( ( cut = l.bridge( r ) ) > limit ) continue;
        ii = i;
        paired[0] = &l;
        paired[1] = &r;
        limit = cut;
    }
    
    if ( !paired[0] && !paired[1] && ii >= 0 ) paired[dd] = &scythes[ii][dd][jj];
    if ( !paired[0] && !paired[1] ) return false;
    
    for ( int i : { 0, 1 } ) node_[i]->test();
    
    bool keep[2]{ false, false };
//    int base[2]{ node_[ii]->coords_[0] + ( paired[0] && paired[0]->keep_ ? paired[0]->coords_[0][1] : 0 )
//               , node_[ii]->coords_[0] + ( paired[1] && paired[1]->keep_ ? paired[1]->coords_[0][0] : (int)node_[ii]->seq_.size() ) };
    Coords base[2]{ Coords( 0, 0 ), Coords( node_[ii]->seq_.size(), node_[ii]->seq_.size() ) }, uncut( node_[ii]->seq_.size(), node_[ii]->seq_.size() );
    for ( int i : { 0, 1 } ) if ( paired[i] && paired[i]->keep_ ) base[i] = paired[i]->coords_[0];
    for ( int i : { 0, 1 } ) base[i] += node_[ii]->coords_[0];
    for ( int i = 0; i < node_[ii]->reads_.size(); i++ )
    {
        bool trash = true;
        if ( node_[ii]->reads_[i].coords_[1] <= base[0][1] && !( trash = false ) && node_[ii]->reads_[i].coords_[0] <= base[0][0] ) keep[0] = true;
        if ( base[1][0] <= node_[ii]->reads_[i].coords_[0] && !( trash = false ) && base[1][1] <= node_[ii]->reads_[i].coords_[1] ) keep[1] = true;
        if ( trash ) node_[ii]->reads_.erase( node_[ii]->reads_.begin() + i-- );
        else for ( int d : { 0, 1 } ) uncut[d] = min( uncut[d], abs( node_[ii]->reads_[i].coords_[d] - node_[ii]->coords_[d] ) );
    }
    
    int ll[2]{ 0, keep[0] ? paired[0]->coords_[0][0] : 0 };
    int mm[2]{ keep[0] ? paired[0]->coords_[1][0] : 0, keep[1] ? paired[1]->coords_[1][1] : (int)node->seq_.size() };
    int rr[2]{ keep[1] ? paired[1]->coords_[0][1] : (int)node_[ii]->seq_.size(), (int)node_[ii]->seq_.size() };
    int alt[2]{ node->coords_[0] + ( keep[0] ? paired[0]->coords_[1][0] : 0 ), ( keep[1] ? node->coords_[0] + paired[1]->coords_[1][1] : node->coords_[1] ) };
    int diffMin = keep[0] && keep[1] ? paired[1]->coords_[0][0] : 0;
    int diffBase = keep[0] && keep[1] ? ( mm[1]-mm[0] ) - ( rr[0]-ll[1] ) : 0;
//    int diffAlt = node_[ii]->coords_[0] - node->coords_[0] + ( paired[0] ? paired[0]->coords_[0][1] - paired[0]->coords_[1][1]
//                                                                         : paired[1]->coords_[0][0] - paired[1]->coords_[1][0] );
    int diffAlt = node_[ii]->coords_[0] - node->coords_[0] + ( keep[0] || !paired[1] ? paired[0]->coords_[0][1] - paired[0]->coords_[1][1]
                                                                                     : paired[1]->coords_[0][0] - paired[1]->coords_[1][0] );
    assert( ll[0]<=ll[1] && mm[0]<=mm[1] && rr[0]<=rr[1] );
    node_[ii]->seq_ = node_[ii]->seq_.substr( 0, ll[1] ) + node->seq_.substr( mm[0], mm[1]-mm[0] ) + node_[ii]->seq_.substr( rr[0], rr[1]-rr[0] );
    int cutted[2]{0}, added = 0;
    if ( diffBase ) for ( Mapped& read : node_[ii]->reads_ ) if ( read.coords_[0] >= diffMin ) read.coords_ += diffBase;
    for ( Mapped& read : node->reads_ ) if ( alt[0] <= read.coords_[0] && read.coords_[1] <= alt[1] )
    {
        node_[ii]->reads_.push_back( Mapped( read.id_, Coords( read.coords_[0]+diffAlt, read.coords_[1]+diffAlt ) ) );
        added++;
    }
    if ( !keep[0] ) node_[ii]->coords_[0] = node->coords_[0] + diffAlt;
    node_[ii]->coords_[1] = keep[1] ?  node_[ii]->coords_[1] + diffBase : node->coords_[1] + diffAlt;
    for ( int d : { 0, 1 } ) if ( keep[d] && uncut[d] ) node_[ii]->trim( uncut[d], d );
    for ( int i : { 0, 1 } ) node_[i]->test();
    for ( int d : { 0, 1 } ) extended[d] = !keep[d];
    assert( added );
    
    return true;
    
//    struct Tmp
//    {
//        Tmp( Coords bc, Coords ac, unordered_set<ReadId> bi, unordered_set<ReadId> ai ):baseCoords( bc ), altCoords( ac ), baseIds( bi ), altIds( ai ){};
//        Coords baseCoords, altCoords;
//        unordered_set<ReadId> baseIds, altIds;
//    };
//    
//    
//    vector<Tmp> aligns[2][2];
//    Tmp* best[2]{ NULL, NULL };
//    for ( int i : { 0, 1 } ) for ( pair<Coords, Coords>& aligned : Coords::align( node_[i]->seq_, node->seq_, 60 ) ) for ( int d : { 0, 1 } )
//    {
//        unordered_set<ReadId> ids = Mapped::getIds( node->reads_, node->coords_[0]+aligned.second[d], d, d );
//        if ( ids.size() < 3 ) aligns[i][d].push_back( Tmp( aligned.first, aligned.second, Mapped::getIds( node_[i]->reads_, aligned.first[!d], !d, !d ), ids ) );
//    }
//
//    int j, cuts, merges = 0;
//    
//    
//    // Try single merge
//    for ( int i : { 0, 1 } ) for ( int d : { 0, 1 } ) for ( Tmp& align : aligns[i][d] )
//    {
//        int counted = 0;
//        for ( ReadId id : align.baseIds ) if ( baseIds[!i].find( id ) == baseIds[!i].end() ) counted++;
//        if ( counted > 2 || ( best[0] && cuts < ( counted+align.altIds.size() ) ) ) continue;
//        j = i;
//        dd = d;
//        cuts = align.baseIds.size()+align.altIds.size();
//        best[0] = &align;
//        merges = 1;
//    }
//    
//    // Try double merge
//    for ( int i : { 0, 1 } ) for ( Tmp& l : aligns[i][0] ) for ( Tmp& r : aligns[i][1] ) if ( l.baseCoords[0] < r.baseCoords[0] )
//    {
//        int missed = 0;
//        for ( ReadId id : l.baseIds ) if ( r.baseIds.find( id ) != r.baseIds.end() && baseIds[!i].find( id ) == baseIds[!i].end() ) missed++;
//        if ( missed > 2  ) continue;
//        missed += l.altIds.size()+r.altIds.size();
//        if ( best[1] && cuts < missed ) assert( false );
//        if ( best[1] && cuts < missed ) continue;
//        j = i;
//        cuts = missed;
//        best[0] = &l;
//        best[1] = &r;
//        merges = 2;
//    }
//    
//    if ( merges == 2 )
//    {
//        assert( paired[0] && paired[1] );
//        node_[j]->seq_ = node_[j]->seq_.substr( 0, best[0]->baseCoords[1] )
//                       + node->seq_.substr( best[0]->altCoords[1], best[1]->altCoords[1] - best[0]->altCoords[1] )
//                       + node_[j]->seq_.substr( best[1]->baseCoords[1] );
//        int diff = ( best[1]->baseCoords[0]-best[0]->baseCoords[1] ) - ( best[1]->altCoords[0]-best[0]->altCoords[1] );
//        for ( int i = 0; i < node_[j]->reads_.size(); i++ )
//        {
//            if ( best[1]->baseIds.find( node_[j]->reads_[i].id_ ) != best[1]->baseIds.end() ) node_[j]->reads_[i].coords_ += diff;
//            else if ( best[1]->baseIds.find( node_[j]->reads_[i].id_ ) != best[1]->baseIds.end() ) node_[j]->reads_.erase( node_[j]->reads_.begin() + i-- );
//        }
//        diff = best[0]->baseCoords[0]-node->coords_[0];
//        for ( Mapped& read : node->reads_ ) node_[j]->reads_.push_back( Mapped( read.id_, Coords( read.coords_[0]+diff, read.coords_[1]+diff ) ) );
//        node_[j]->coords_[1] = node_[j]->seq_.size();
//        return true;
//    }
//    else assert( !paired[0] && !paired[1] );
//    
//    if ( merges == 1 )
//    {
//        assert( ii >= 0 && jj >= 0 );
//        node_[j]->seq_ = dd ? node->seq_.substr( 0, best[0]->altCoords[0] ) + node_[j]->seq_.substr( best[0]->baseCoords[0] )
//                             : node_[j]->seq_.substr( 0, best[0]->baseCoords[1] ) + node->seq_.substr( best[0]->altCoords[1] );
//        int diff = dd ? -node->coords_[0] : best[0]->baseCoords[0] - best[0]->altCoords[0];
//        if ( dd ) for ( Mapped& read : node_[j]->reads_ ) read.coords_ += best[0]->altCoords[0];
//        for ( Mapped& read : node->reads_ ) node_[j]->reads_.push_back( Mapped( read.id_, Coords( read.coords_[0]+diff, read.coords_[1]+diff ) ) );
//        node_[j]->coords_[1] = node_[j]->seq_.size();
//        extended[!dd] = true;
//        return true;
//    }
//    else assert( ii < 0 && jj < 0 );
//    
//    return false;
}

bool Locus::confirm( float similarity, int readLen, int minLen, float maxRatio )
{
    len_ = min( node_[0]->seq_.size(), node_[1]->seq_.size() );
    similarity_ = getSimilarity();
    if ( similarity_ < similarity ) return false;
    return setCoverage( readLen, maxRatio );
}

void Locus::cull( vector<Node*>& nodes )
{
    for ( Node* node : nodes ) for ( int d : { 0, 1 } ) if ( !node->edges_[d].empty() )
    {
        int minOl = node->seq_.size(), ol = 0, coords[2];
        for ( pair<Node*, int> edge : node->edges_[d] ) minOl = min( minOl, edge.second );
        for ( int i : { 0, 1 } ) if ( mapSeqEnd( node->seq_, node_[i]->seq_, minOl, coords, d ) ) ol = max( ol, coords[1]-coords[0] );
        if  ( ol ) for ( int i = 0; i < node->edges_[d].size(); i++ ) if ( node->edges_[d][i].second <= ol )
        {
            node->edges_[d][i].first->removeEdge( node, !d );
            node->edges_[d].erase( node->edges_[d].begin() + i-- );
        }
    }
}

bool Locus::extend( unordered_set<Node*>& pathed, vector<Node*>& nodes, int readLen, bool drxn )
{
    bool extended = false;
    Coords base[2], ext[2];
    Node* branch[2]{ NULL, NULL };
    int coord[2]{ node_[0]->coords_[drxn] - node_[0]->coords_[0], node_[1]->coords_[drxn] - node_[1]->coords_[0] };
    int limit[2]{ coord[0], coord[1] };
    for ( int i : { 0, 1 } ) node_[i]->rebase();
    for ( int i : { 0, 1 } ) Mapped::sort( node_[i]->reads_, true, drxn );
    for ( int i : { 0, 1 } ) if ( node_[i]->reads_.size() > 1 ) coord[i] = node_[i]->reads_[ drxn ? node_[i]->reads_.size()-2 : 1 ].coords_[drxn];
    for ( int i : { 0, 1 } ) for ( Node* node : nodes ) for ( pair<Coords, Coords> hit : Coords::align( node_[i]->seq_, node->seq_, readLen*.6 ) )
    {
        if ( hit.first.len() >= readLen ) pathed.insert( node );
        if ( drxn ? hit.first[1] < coord[i] : coord[i] < hit.first[0] ) continue;
        if ( drxn ? (int)node->seq_.size()-readLen < hit.second[0] : hit.second[1] < readLen ) continue;
        if ( branch[i] && base[i][drxn] == limit[i] && hit.first[drxn] != limit[i] ) continue;
        if ( branch[i] && ( ( base[i][drxn] == limit[i] ) == ( hit.first[drxn] == limit[i] ) ) && base[i].len() > hit.first.len() ) continue;
        branch[i] = node;
        base[i] = hit.first;
        ext[i] = Coords( hit.second[0] + node->coords_[0], hit.second[1] + node->coords_[0] );
    }
    for ( int i : { 0, 1 } ) if ( branch[i] )
    {
        int reads = 0;
        for ( Mapped& read : branch[i]->reads_ )
        {
            if ( drxn ? ext[i][0] <= read.coords_[0] && ext[i][1] < read.coords_[1] : read.coords_[0] < ext[i][0] && read.coords_[1] <= ext[i][1] ) reads++;
        }
        if ( reads < 2 ) for ( pair<Node*, int> edge : branch[i]->edges_[drxn] ) if ( edge.first->countReads( reads, 2, drxn ) >= 2 ) reads = 2;
        if ( reads < 2 ) branch[i] = NULL;
    }
    
    vector<Node*> path[2], tmp[2];
    if ( branch[0] && branch[1] ) for ( int i : { 0, 1 } ) if ( path[0].empty() && path[1].empty() && branch[i]->getDiploid( tmp ) )
    {
        for ( int j : { 0, 1 } ) if ( branch[j] == ( drxn ? tmp[j][0] : tmp[j].back() ) )
        {
            path[j] = tmp[j];
        }
    }
    for ( int i : { 0, 1 } ) if ( branch[i] && path[i].empty() )
    {
        path[i] = branch[i]->getHaploid( drxn );
    }
    for ( int i : { 0, 1 } ) if ( branch[i] && !path[i].empty() )
    {
        path[i].erase( drxn ? path[i].begin() : path[i].end()-1 );
        pathed.insert( path[i].begin(), path[i].end() );
        if ( branch[i]->coords_[drxn] != ext[i][drxn] || !path[i].empty() ) extended = true;
        node_[i]->trim( abs( node_[i]->coords_[drxn] - base[i][drxn] ), drxn );
        node_[i]->add( drxn ? branch[i]->seq_.substr( ext[i][1]-branch[i]->coords_[0] )
                            : branch[i]->seq_.substr( 0, ext[i][0]-branch[i]->coords_[0] ), true, drxn );
        if ( path[i].empty() ) continue;
        int ol = 0;
        for ( pair<Node*, int> edge : branch[i]->edges_[drxn] ) if ( edge.first == ( drxn ? path[i][0] : path[i].back() ) ) ol = edge.second;
        string seq = Node::getSeq( path[i] );
        node_[i]->add( drxn ? seq.substr( ol ) : seq.substr( 0, seq.size()-ol ), true, drxn );
    }
    return extended;
    
//    Node* branch[2]{ NULL, NULL },* bridge[2]{ NULL, NULL };
//    int coord[2], tmp[2];
//    for ( int i : { 0, 1 } ) for ( Node* node : pathed ) if ( mapSeqEnd( node->seq_, node_[i]->seq_, readLen, tmp, drxn ) )
//    {
//        if ( branch[i] && ( drxn ? tmp[1] < coord[i] : coord[i] < tmp[0] ) ) continue;
//        coord[i] = tmp[drxn];
//        branch[i] = node;
//    }
//    
//    for ( int i : { 0, 1 } ) for ( Node* node : pathed ) for ( pair<Coords, Coords> hit : Coords::align( node_[i]->seq_, node->seq_, readLen ) )
//    {
//        if ( branch[i] && ( drxn ? hit.first[1] < coord[i] : coord[i] < hit.first[0] ) ) continue;
//        coord[i] = hit.first[drxn];
//        branch[i] = node;
//    }
//    for ( int i : { 0, 1 } ) if ( branch[i] && coord[i] == node_[i]->coords_[drxn] ) bridge[i] = branch[i];
//    for ( int i : { 0, 1 } ) if ( branch[i] && !bridge[i] )
//    {
//        Coords base, ext;
//        for ( Node* node : nodes ) if ( pathed.find( node ) == pathed.end() ) for ( pair<Coords, Coords> hit : Coords::align( node_[i]->seq_, node->seq_, readLen*.6 ) )
//        {
//            if ( hit.second[1]-hit.second[0] >= min( readLen, (int)node->seq_.size() ) ) pathed.insert( node );
//            if ( drxn ? hit.first[1] < coord[i] : coord[i] < hit.first[0] ) continue;
//            if ( !node->edges_[!drxn].empty() ) continue;
//            if ( bridge[i] && ( drxn ? base[0] < hit.first[0] : hit.first[1] < base[1] ) ) continue;
//            if ( Mapped::getIds( node->reads_, hit.second[!drxn] + node->coords_[0], !drxn, !drxn ).size() > 1 ) continue;
//            if ( drxn ? hit.second[1] < node->seq_.size() && hit.first[1] < node_[i]->seq_.size() : hit.second[0] && hit.first[0] )
//            {
//                // too many reads would be lost from the branch
//                if ( Mapped::getIds( node_[i]->reads_, hit.first[drxn] + node_[i]->coords_[0], drxn, drxn ).size() > 1 ) continue;
//            }
//            base = hit.first;
//            ext = hit.second;
//            bridge[i] = node;
//        }
//        if ( bridge[i] && ( drxn ? ext[1] < bridge[i]->seq_.size() : ext[0] ) )
//        {
//            node_[i]->trim( drxn ? node_[i]->seq_.size() - base[1] : base[0], drxn );
//            node_[i]->add( drxn ? bridge[i]->seq_.substr( ext[1] ) : bridge[i]->seq_.substr( 0, ext[0] ), true, drxn );
//            extended = true;
//        }
//        if ( !bridge[0] || !bridge[1] ) continue;
//        vector<Node*> path[2];
//        if ( !bridge[i]->getDiploid( path ) ) break;
//        if ( bridge[i] != ( drxn ? path[i][0] : path[i].back() ) ) path[0].swap( path[1] );
//        if ( bridge[0] != ( drxn ? path[0][0] : path[0].back() ) || bridge[1] != ( drxn ? path[1][0] : path[1].back() ) ) break;
//        
//        int ol[2]{0}, coords[2][2];
//        for ( int j : { 0, 1 } ) pathed.insert( path[j].begin(), path[j].end() );
//        for ( int j : { 0, 1 } ) for ( pair<Node*, int> edge : ( drxn ? path[j][1] : path[j].end()[-2] )->edges_[!drxn] ) if ( edge.first == bridge[j] ) ol[j] = edge.second;
//        for ( int j : { 0, 1 } ) path[j].erase( drxn ? path[j].begin() : path[j].end()-1 );
//        string seq[2]{ Node::getSeq( path[0] ), Node::getSeq( path[1] ) };
//        for ( int j : { 0, 1 } ) assert( mapSeqEnd( node_[j]->seq_, seq[j], ol[j], coords[j], drxn ) );
//        for ( int j : { 0, 1 } ) node_[j]->add( drxn ? seq[j].substr( coords[j][1] ) : seq[j].substr( 0, coords[j][0] ), true, drxn );
//        node_[0]->test();
//        node_[1]->test();
//        return true;
//    }
//    return extended;
}

void Locus::fill( vector<Read> reads, vector<Read>& unmapped, int i, bool extended[2] )
{
    unordered_set<ReadId> used;
    Exts exts[2]{ Exts( node_[i]->seq_, node_[i]->seq_.size(), 0 ), Exts( node_[i]->seq_, node_[i]->seq_.size(), 1 ) };
    vector< pair<Read, int> > ols[2];
    unordered_map<ReadId, int> mapped;
    for ( Read& read : unmapped ) used.insert( read.id_ );
    for ( Read& read : reads )
    {
        bool hit = false;
        for ( size_t it = node_[i]->seq_.find( read.seq_ ); it != string::npos && ( hit = true ); )
        {
            auto ins = mapped.insert( make_pair( read.id_, it ) );
            if ( ins.second || ins.first->second < it )
            {
                node_[i]->reads_.push_back( Mapped( read.id_, Coords( it, it + read.coords_[1] - read.coords_[0] ) ) );
                ins.first->second = it;
                break;
            }
            it = node_[i]->seq_.find( read.seq_, ins.first->second+1 );
        }
        if ( !hit && ( node_[!i]->seq_.find( read.seq_ ) != string::npos ) ) hit = true;
        else if ( !hit )
        {
            int ol = 0;
            for ( int d : { 0, 1 } ) if ( ( ol = mapSeqOverlap( node_[i]->seq_, read.seq_, 50, d ) ) && ( hit = true ) )
            {
                ols[d].push_back( make_pair( read, ol ) );
//                exts[d].add( exts[d].exts_, read.seq_, read.id_, ol, d );
            }
            if ( !hit ) for ( int d : { 0, 1 } ) if ( ol = mapSeqOverlap( read.seq_, node_[i]->seq_, 50, d ) ) hit = true;
            if ( read.coords_[0] < 0 || node_[i]->seq_.size() < read.coords_[1] ) hit = true;

        }
        if ( !hit && used.insert( read.id_ ).second ) unmapped.push_back( read );
    }
    if ( node_[i]->seq_.size() < 300 ) for ( int d : { 0, 1 } )
    {
        sort( ols[d].begin(), ols[d].end(), []( pair<Read, int>& a, pair<Read, int>& b ){ return a.second > b.second; });
        for ( pair<Read, int>& read : ols[d] ) exts[d].add( exts[d].exts_, read.first.seq_, read.first.id_, read.second, d );
        Ext::cull( exts[d].exts_, 2, true, d );
        if ( exts[d].exts_.empty() || exts[d].exts_.size() > 2 ) continue;
        Ext* alt = NULL;
        int minOl = 30, ol;
        for ( Ext* e : exts[d].exts_ )
        {
            e->set( node_[i]->seq_, d );
            if ( ( ol = mapSeqOverlap( node_[!i]->seq_, e->seq_, minOl, d ) ) && ( alt = e ) ) minOl = ol;
        }
        if ( !alt && exts[d].exts_.size() > 1 ) continue;
        for ( Ext* e : exts[d].exts_ ) if ( e != alt || exts[d].exts_.size() == 1 )
        {
            for ( ExtRead& er : e->reads_ ) node_[i]->reads_.push_back( 
                    Mapped( er.id_, Coords( node_[i]->coords_[d] - ( d ? er.ol_ : er.ext_ ), node_[i]->coords_[d] + ( d ? er.ext_ : er.ol_ ) ) ) );
            node_[i]->add( e->ext_, true, d );
            break;
        }
        if ( alt )
        {
            int diff = ol - ( alt->reads_[0].ol_ );
            for ( ExtRead& er : alt->reads_ ) node_[!i]->reads_.push_back( 
                    Mapped( er.id_, Coords( node_[!i]->coords_[d] - ( d ? er.ol_+diff : er.ext_-diff ), node_[!i]->coords_[d] + ( d ? er.ext_-diff : er.ol_+diff ) ) ) );
            node_[!i]->add( d ? alt->seq_.substr( ol ) : alt->seq_.substr( 0, alt->seq_.size()-ol ), true, d );
        }
        
        extended[d] = true;
    }
    node_[i]->test( true );
}

float Locus::getSimilarity()
{
    vector< pair<Coords, Coords> > coords = Coords::align( node_[0]->seq_, node_[1]->seq_, 15 );
    if ( coords.empty() ) return 0;
    sort( coords.begin(), coords.end(), []( pair<Coords, Coords>& a, pair<Coords, Coords>& b ){ return a.first[0] < b.first[0]; } ); 
    for ( int i = 0; i < coords.size(); i++ ) for ( int j = i+1; j < coords.size(); j++ ) if ( Coords::conflict( coords[i], coords[j] ) )
    {
        int lens[2]{ coords[i].first.len(), coords[j].first.len() };
        for ( int k = j+1; k < coords.size(); k++ )
        {
            bool conflicted[2]{ Coords::conflict( coords[i], coords[k] ), Coords::conflict( coords[j], coords[k] ) };
            if ( conflicted[0] && !conflicted[1] ) lens[1] += coords[k].first.len();
            if ( !conflicted[0] && conflicted[1] ) lens[0] += coords[k].first.len();
        }
        if ( lens[1] > lens[0] )
        {
            coords.erase( coords.begin() + i-- );
            break;
        }
        else coords.erase( coords.begin() + j-- );
    }
    coords.push_back( make_pair( Coords( node_[0]->seq_.size(), node_[0]->seq_.size() ), Coords( node_[1]->seq_.size(), node_[1]->seq_.size() ) ) );
    string base[2], seq[2];
    Coords c[2][2];
    for ( int i = 0; i < coords.size(); i++ )
    {
        c[1][0] = coords[i].first;
        c[1][1] = coords[i].second;
        for ( int j : { 0, 1 } ) seq[j] = node_[j]->seq_.substr( c[0][j][1], c[1][j][0] - c[0][j][1] );
        if ( !seq[0].empty() && !seq[1].empty() && ( seq[0].size() + seq[1].size() ) > 2 ) LocalAlignment( seq[0], seq[1], false, false ).realign( seq[0], seq[1], false );
        for ( int j : { 0, 1 } ) if ( seq[j].empty() && !seq[!j].empty() ) seq[j] = string( seq[!j].size(), '-' );
        for ( int j : { 0, 1 } ) base[j] += seq[j];
        for ( int j : { 0, 1 } ) base[j] += node_[j]->seq_.substr( c[1][j][0], c[1][j][1] - c[1][j][0] );
        for ( int j : { 0, 1 } ) c[0][j] = c[1][j];
    }
    assert( base[0].size() == base[1].size() );
    int ends[2]{ 0, (int)base[0].size() }, hit = 0, miss = 0;
    for ( int i : { 0, 1 } ) if ( !ends[0] && base[i][0] == '-' ) while( ends[0] < base[i].size() && base[i][ ends[0] ] == '-' ) ends[0]++;
    for ( int i : { 0, 1 } ) if ( ends[1] == base[i].size() && base[i][ ends[1]-1 ] == '-' ) while( ends[1] > 0 && base[i][ ends[1]-1 ] == '-' ) ends[1]--;
    for ( int i = ends[0]; i < ends[1]; i++ ) ( base[0][i] == base[1][i] ? hit : miss )++;
    return float( hit * 100 ) / (float)max( 1, hit+miss );
}

bool Locus::setCoverage( int readLen, float maxRatio )
{
    for ( int i : { 0, 1 } ) if ( node_[i]->seq_.size()+1 < readLen*2 ) return false;
    float cover[2]{0};
    unordered_map<ReadId, int> ids;
    for ( int i : { 0, 1 } ) for ( Mapped& read : node_[i]->reads_ )
    {
        auto ins = ids.insert( make_pair( read.id_, 1 ) );
        if ( !ins.second ) ins.first->second++;
    }
    for ( int i : { 0, 1 } )
    {
        int len = node_[i]->seq_.size(), diff = -node_[i]->coords_[0];
        float base[len]{0};
        for ( Mapped& read : node_[i]->reads_ ) assert( read.coords_[0]+diff >= 0 && read.coords_[1]+diff <= len );
        for ( Mapped& read : node_[i]->reads_ )
        {
            float multi = 1 / (float)ids[read.id_];
            for ( int j = read.coords_[0]+diff; j < read.coords_[1]+diff; j++ ) base[j] += multi;
        }
        vector<float> based;
        for ( int j = readLen-1; j < len+1-readLen; j++ ) based.push_back( base[j] );
        sort( based.begin(), based.end() );
        cover[i] = ( based[ ( based.size()-1 ) / 2 ] + based[ based.size() / 2 ] ) / 2;
    }
    coverage_ = cover[0] + cover[1];
    if ( max( ( cover[1] / cover[0] ), ( cover[0] / cover[1] ) ) > maxRatio ) return false;
    return true;
}

bool Locus::setLen( int minLen )
{
    len_ = min( node_[0]->seq_.size(), node_[1]->seq_.size() );
    return len_ >= minLen;
}

bool Locus::setSimilarity( float minSimilarity )
{
    similarity_ = getSimilarity();
    return similarity_ >= minSimilarity;
}

