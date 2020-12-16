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

#include "coverage.h"
#include "filenames.h"
#include "match_query.h"
#include <ctime> 
#include <iostream>
#include <string.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include "query_overlap.h"
#include <sys/stat.h>
#include <chrono>
#include <iomanip>

extern Parameters params;

Coverage::Coverage( int argc, char** argv )
:ir_( NULL ), qb_( NULL )
{
    queried_ = unmatched_ = nondiploid_ = overtimed_ = miscovered_ = dissimilar_ = short_ = testFail_ = 0;
    string ifn, ofn;
    Filenames* fns = NULL;
    
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
        }
        else if ( !strcmp( argv[i], "-i" ) )
        {
            if ( !ifn.empty() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            ifn = argv[++i];
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( !ofn.empty() )
            {
                cerr << "Error: multiple outputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            ofn = argv[++i];
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one output prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            string prefix = argv[++i];
            if ( prefix[0] != '/' )
            {
                string curr = getcwd( NULL, 0 );
                if ( !prefix.empty() && prefix[0] == '.' ) prefix = prefix.substr( 1 );
                if ( !prefix.empty() && prefix[0] != '/' ) prefix = "/" + prefix;
                if ( prefix.empty() || curr.empty() )
                {
                    cerr << "Invalid file prefix. Please use the absolute path." << endl;
                }
                prefix = curr + prefix;
            }
            
            fns = new Filenames( prefix );
        }
    }
    
    ir_ = new IndexReader( fns );
    qb_ = new QueryBinaries( fns );
    
    if ( ofn.empty() ) ofn = "./match_result.fa";
    if ( ifn.empty() )
    {
        cerr << "Error: no query sequence(s) provided." << endl;
        exit( EXIT_FAILURE );
    }
    
    ifstream ifs( ifn );
    string line, seq;
    int i = 0;
    while ( getline( ifs, line ) && covers_.size() < 1000 )
    {
        i++;
        if ( !line.empty() && line.size() >= 300 ) seed( line );
    }
    ifs.close();
    
    ofstream ofs( ofn );
    for ( float cover : coverage_ ) ofs << std::fixed << std::setprecision(2) << cover << "\n";
    sort( coverage_.begin(), coverage_.end() );
    cout << "Final median coverage:   " << ( coverage_[ ( coverage_.size()-1 ) / 2 ] + coverage_[ coverage_.size() / 2 ] ) / 2 << endl;
    cout << "Total queried coding sequences:   " << queried_ << endl;
    cout << "Unmatched queries:   " << unmatched_ << endl;
    cout << "Non-diploid loci:   " << nondiploid_ - testFail_ << endl;
    cout << "Loci with allelic dissimilarity in coverage (+/- 50%):   " << miscovered_ << endl;
    cout << "Loci with allelic dissimilarity in sequence (<97% similarity):   " << dissimilar_ << endl;
    cout << "Short loci:   " << short_ << endl;
    ofs.close();
}

bool Coverage::confirm( Locus* locus, unordered_set<Node*>& pathed, vector<Node*>& nodes )
{
    bool failure = false, extended[2]{ false, false }, merged[2]{ false, false };
    vector<Read> unmapped;
    for ( int i : { 0, 1 } ) locus->node_[i]->reads_.clear();
    for ( int i : { 0, 1 } ) locus->node_[i]->rebase();
    for ( int i : { 0, 1 } )
    {
        MatchQuery mq( locus->node_[i]->seq_, ir_, 5 );
        if ( mq.failure_ )
        {
            cout << "Mq fail" << endl;
            return false;
        }
        locus->fill( mq.yield( qb_ ), unmapped, i, extended );
    }
    vector<Node*> alts = Node::getNodes( unmapped, false, failure );
    if ( failure )
    {
        cout << "Not got" << endl;
        return false;
    }
    locus->cull( alts );
    cull( alts );
    for ( int i : { 0, 1 } ) locus->node_[i]->test();
    if ( !extended[0] && !extended[1] ) for ( int i = 0; i < alts.size(); i++ ) if ( locus->add( alts[i], merged ) || !alts[i]->isSubstantial( 4 ) )
    {
        delete alts[i];
        alts.erase( alts.begin() + i-- );
    }
//    if ( !alts.empty() && !extended[0] && !extended[1] )
//    {
//        cout << ">A" << endl << locus->node_[0]->seq_ << endl;
//        cout << ">B" << endl << locus->node_[1]->seq_ << endl;
//        for ( Node* node : alts ) cout << ">" << node->reads_.size() << endl << node->seq_ << endl;
//        int x = 0;
//    }
    for ( Node* node : alts ) delete node;
    for ( int i : { 0, 1 } ) if ( !locus->node_[i]->test() )
    {
        testFail_++;
        return false;
    }
    if ( !alts.empty() ) return false;
    for ( int d : { 0, 1 } ) if ( merged[d] ) merged[d] = locus->extend( pathed, nodes, params.readLen, d );
    if ( extended[0] || extended[1] || merged[0] || merged[1] ) return confirm( locus, pathed, nodes );
    for ( int i : { 0, 1 } ) locus->node_[i]->test();
    return true;
}

void Coverage::cull( vector<Node*>& nodes )
{
    Node::cull( nodes );
    
    for ( int again = 1; again--; ) for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->reads_.size() < ( 1 + ( nodes[i]->seq_.size() * 2 / params.readLen ) ) )
    {
        for ( int d : { 0, 1 } ) if ( nodes[i]->isWeak( params.readLen * .6, d ) )
        {
            string seq = nodes[i]->seq_;
            if ( seq.size() > params.readLen ) seq = d ? seq.substr( seq.size()-params.readLen ) : seq.substr( 0, params.readLen );
            bool bad = QueryOverlap::countOverlaps( seq, ir_, params.readLen*2/3, d ) < 3;
            if ( !bad && ( seq.size() == nodes[i]->seq_.size() || QueryOverlap::countOverlaps( seq, ir_, params.readLen*2/3, !d ) > 2 ) ) continue;
            delete nodes[i];
            nodes.erase( nodes.begin() + i-- );
            break;
        }
    }
    
    Node::simplify( nodes );
    
    for ( int again = 1; again--; ) for ( int i = 0; i < nodes.size(); i++ ) if ( nodes[i]->reads_.size() < 3 )
    {
        bool bad[2]{ nodes[i]->edges_[0].empty(), nodes[i]->edges_[1].empty() };
        unordered_set<Node*> pathed[2];
        int counts[2]{ nodes[i]->countReads( 0, 10, 0 ), nodes[i]->countReads( 0, 10, 1 ) };
        for ( int d : { 0, 1 } ) for ( pair<Node*, int> edge : nodes[i]->edges_[d] ) for ( pair<Node*, int> re : edge.first->edges_[!d] )
        {
            if ( re.first != nodes[i] && re.first->isSubstantial( counts[!d]+2, !d ) ) bad[d] = true;
        }
        if ( !bad[0] || !bad[1] ) continue;
        delete nodes[i];
        nodes.erase( nodes.begin() + i-- );
    }
    
    Node::simplify( nodes );
}

vector<Node*> Coverage::getNodes( vector<Read>& reads )
{
    bool failure = false;
    if ( reads.size() >= 2000 ) return vector<Node*>();
    vector<Node*> nodes = Node::getNodes( reads, true, failure );
    if ( failure ) assert( nodes.empty() );
    else cull( nodes );
    return nodes;
}

vector<Locus*> Coverage::getSeeds( vector<Node*>& nodes )
{
    auto t_start = std::chrono::high_resolution_clock::now();
    unordered_set<Node*> added, tried;
    vector< pair< vector<Node*>, int> > paths = Node::getPaths( nodes );
    vector<Locus*> loci;
    bool tooShort = false, miscovered = false, nondiploid = false, dissimilar, overtime = false;
    for ( int i = 0; !nondiploid && !overtime && i < paths.size(); i++ )
    {
        vector<Node*> path[2];
        for ( Node* node : paths[i].first ) if ( tried.find( node ) == tried.end() && node->getDiploid( path ) )
        {
            unordered_set<Node*> pathed, unpathed;
            for ( int j : { 0, 1 } ) pathed.insert( path[j].begin(), path[j].end() );
            Coords coords = path[0][0]->coords_;
            for ( Node* n : pathed ) coords[0] = min( coords[0], n->coords_[0] );
            for ( Node* n : pathed ) coords[1] = max( coords[1], n->coords_[1] );
            for ( int j : { 0, 1 } ) path[j][0]->get( unpathed, false, 0 );
            for ( int j : { 0, 1 } ) path[j].back()->get( unpathed, false, 1 );
            
            Locus* locus = new Locus( Node::getSeq( path[0] ), Node::getSeq( path[1] ) );
            
            bool good = true;
            for ( Node* n : pathed ) if ( added.find( n ) != added.end() ) good = false;
            if ( good && !( good = confirm( locus, pathed, nodes ) ) )  nondiploid = true;
            if ( good && !( good = locus->setLen( 300 ) ) ) tooShort = true;
            if ( good && !( good = locus->setCoverage( params.readLen, 1.5 ) ) ) miscovered = true;
            if ( good && !( good = locus->setSimilarity( 97 ) ) ) dissimilar = true;
            if ( good ) for ( Node* n : pathed ) if ( !added.insert( n ).second ) good = false;
            tried.insert( pathed.begin(), pathed.end() );
            
            if ( good ) loci.push_back( locus );
            else delete locus;
            if ( nondiploid ) break;
            if ( double( ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC ) > 5 ) nondiploid = overtime = true;
        }
    }
    sort( loci.begin(), loci.end(), []( Locus* a, Locus* b ){ return a->len_ > b->len_; } );
    for ( Locus* locus : loci ) for ( int i : { 0, 1 } ) locus->node_[i]->rebase();
    for ( Locus* locus : loci ) for ( int i : { 0, 1 } ) for ( Mapped& read : locus->node_[i]->reads_ )
    {
        string seq = qb_->getSequence( read.id_ );
        size_t it = locus->node_[i]->seq_.find( seq );
        if ( it == string::npos || read.coords_[0] != it ) assert( false );
    }
    if ( overtime ) for ( Locus* locus : loci ) delete locus;
    if ( overtime ) loci.clear();
    
    if ( loci.empty() )
    {
        if ( overtime ) overtimed_++;
        else if ( nondiploid ) nondiploid_++;
        else if ( miscovered ) miscovered_++;
        else if ( dissimilar ) dissimilar_++;
        else if ( tooShort ) short_++;
        else
        {
            assert( false );
        }
    }
    
    return loci;
}

void Coverage::seed( string seq )
{
    auto t_start = std::chrono::high_resolution_clock::now();
    MatchQuery mq( seq, ir_, 10 );
    vector<Read> reads = mq.yield( qb_ );
    if ( mq.failure_ || reads.size() > ( seq.size() + 1 - params.readLen ) * 20 )
    {
        nondiploid_++;
        reads.clear();
    }
    queried_++;
    vector<Node*> nodes = getNodes( reads );
    if ( nodes.empty() )
    {
        unmatched_++;
        return;
    }
    
    vector<Locus*> loci = getSeeds( nodes );
    cout << "Mapped in:   " << std::fixed << std::setprecision(2) << double( ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC ) << endl;
    sort( loci.begin(), loci.end(), []( Locus* a, Locus* b ){ return a->len_ > b->len_;} );
    if ( !loci.empty() )
    {
        covers_.push_back( loci[0]->coverage_ );
        coverage_.push_back( loci[0]->coverage_ );
        sort( covers_.begin(), covers_.end() );
        cout << covers_.size() << ")   Coverage: " << loci[0]->coverage_ << "   Length: " << loci[0]->len_ << "   Median: " 
                << ( covers_[ ( covers_.size()-1 ) / 2 ] + covers_[ covers_.size() / 2 ] ) / 2 << endl;
    }
    for ( Locus* locus : loci )
    {
        double reads = 0, shorts = 0;
        for ( int i = 0; i < 2; i++ )
        {
            for ( Mapped& read : locus->node_[i]->reads_)
            {
                string seq = qb_->getSequence( read.id_ );
                size_t it = locus->node_[i]->seq_.find( seq );
                assert( it != string::npos );
                if ( read.coords_[1] - read.coords_[0] < 100 ) shorts++;
                reads++;
            }
        }
        cout << ( reads - shorts ) * 100 / reads << endl;
    }
    for ( Node* node : nodes ) delete node;
    for ( Locus* locus : loci ) delete locus;
}

void Coverage::printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "    coverage" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "    Estimates the error-free read coverage of a shotgun genome sequencing dataset." << endl;
    cout << endl << "Usage:" << endl;
    cout << "    locass index [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "    -p    Prefix for BWT data files." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "    -o    Output file name (default: ./match_result.fa)." << endl;
    cout << "    -i    Input fasta file name." << endl;
}

