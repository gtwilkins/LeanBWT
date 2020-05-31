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
        else if ( !strcmp( argv[i], "-m" ) )
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
    
    for ( int i = 0; i < 20; i++ ) test( ir, qb, 2, i );
    assert( false );
    
    ofstream ofs;
    if ( !ofn.empty() ) ofs.open( ofn );
    
    if ( !ifn.empty() )
    {
        assert( false );
    }
    else if ( !seq.empty() )
    {
        string header = "query";
        match( seq, header, ir, qb, ofs.good() ? &ofs : NULL, min( 15, errors ) );
        assert( false );
    }
    else assert( false );
}

void Match::match( string& seq, string& header, IndexReader* ir, QueryBinaries* qb, ofstream* ofs, int errors )
{
    vector<MatchRead> reads = MatchQuery( seq, ir, errors ).yield( qb );
    sort( reads.begin(), reads.end(), []( MatchRead& a, MatchRead& b ){ return a.coord[0] < b.coord[0]; } );
    int base = !reads.empty() ? max( -reads[0].coord[0], 0 ) : 0;
    cout << ">" << header << "|matched:" << reads.size() << endl;
    cout << string( base, '-' ) << seq << endl;
    for ( MatchRead& mr : reads ) cout << ">read_" << mr.id << endl << string( mr.coord[0]+base, '-' ) << mr.seq << endl;
    assert( false );
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
    
}
