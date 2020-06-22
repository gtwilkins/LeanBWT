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

#include "match.h"
#include "filenames.h"
#include "match_query.h"
#include "parameters.h"
#include "timer.h"
#include "shared_functions.h"
#include <iostream>
#include <string.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>

extern Parameters params;

Match::Match( int argc, char** argv )
{
    string ifn, ofn, seq;
    int errors = 0;
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
            if ( !ifn.empty() || !seq.empty() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            ifn = argv[++i];
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( !ifn.empty() )
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
        else if ( !strcmp( argv[i], "-e" ) )
        {
            errors = stoi( argv[++i] );
            if ( errors > 15 || errors < 0 )
            {
                cerr << "Error: invalid mismatch rate set at: " << errors << "%, must be between 0-15%" << endl;
                exit( EXIT_FAILURE );
            }
        }
        else if ( !strcmp( argv[i], "-s" ) ) seq = argv[++i];
    }
    
    IndexReader* ir = new IndexReader( fns );
    QueryBinaries* qb = new QueryBinaries( fns );
    
    if ( ofn.empty() ) ofn = "./match_result.fa";
    
    ofstream ofs( ofn );
    
    if ( !ifn.empty() )
    {
        assert( false );
    }
    else if ( !seq.empty() )
    {
        string header = "query";
        match( seq, header, ir, qb, ofs.good() ? &ofs : NULL, min( 15, errors ) );
        if ( ofs.good() ) ofs.close();
    }
    else
    {
        cerr << "Error: no query sequence(s) provided." << endl;
        exit( EXIT_FAILURE );
    }
}

void Match::match( string& q, string& header, IndexReader* ir, QueryBinaries* qb, ofstream* ofs, int errors )
{
    vector<Read> reads = MatchQuery( q, ir, errors ).yield( qb );
    Read::sort( reads, true, 0 );
    int base = !reads.empty() ? max( -reads[0].coords_[0], 0 ) : 0;
    if ( ofs ) ( *ofs ) << ">" << header << "|matched:" << reads.size() << endl << string( base, '-' ) << "reads" << q << endl;
    else cout << ">" << header << "|matched:" << reads.size() << endl << string( base, '-' ) << "reads" << q << endl;
    Lib* lib;
    unordered_set<ReadId> used;
    bool addPairs = false;
    for ( Read& r : reads ) if ( used.find( r.id_ ) == used.end() )
    {
        string header = ">read_" + to_string( r.id_ );
        string seq = r.seq_;
        int d, coord = r.coords_[0]+base;
        if ( addPairs )
        {
            lib = params.getLib( r.id_ );
            if ( lib && lib->isPe )
            {
                ReadId id = r.id_;
                lib->getPair( id, d );
                string alt = qb->getSequence( id );
                int ol = mapSeqOverlap( seq, alt, 15, d );
                if ( ol )
                {
                    if ( !d ) coord -= alt.size()-ol;
                    seq = (d?seq:alt) + (d?alt:seq).substr( ol );
                }
                else
                {
                    seq = (d?seq:alt) + string( 100, '-' ) + (d?alt:seq);
                    if ( !d ) coord -= 250;
                }
                used.insert( id );
            }
        }
        assert( coord >= 0 );
        seq = string( coord, '-' ) + seq;
        if ( ofs ) ( *ofs ) << header << endl << seq << endl;
        else cout << header << endl << seq << endl;
    }
}

void Match::test( IndexReader* ir, QueryBinaries* qb, int tests, int errors )
{
    srand( time(NULL) );
    vector<MatchQuery> matches;
    vector< pair<ReadId, string> > queries;
    while ( queries.size() < tests )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = qb->getSequence( id );
        if ( seq.size() == params.readLen ) queries.push_back( make_pair( id, seq ) );
    }
    
    double startTime = clock();
    for ( int i = 0; i < queries.size(); i++ ) matches.push_back( MatchQuery( queries[i].second, ir, errors ) );
    cout << "Errors: " << errors << ", queryies: " << tests << ", time taken: " << getDuration( startTime ) << endl;
}

void Match::printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "    match" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "    Matches reads to one or more query sequences. By default, matches are exact. To find inexact matches, set an allowed mismatch rate (up to 15%)." << endl;
    cout << endl << "Usage:" << endl;
    cout << "    locass index [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "    -p    Prefix for BWT data files." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "    -o    Output file name (default: ./match_result.fa)." << endl;
    cout << "    -i    Input sequence query file (mutually exclusive with -s)." << endl;
    cout << "    -s    Input sequence query (mutually exclusive with -i)." << endl;
    cout << "    -m    Allowed mismatches per 100 bases for inexact matching (default: 0, maximum: 15)." << endl;
}
