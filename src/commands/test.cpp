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

#include "test.h"
#include "parameters.h"
#include "constants.h"
#include "timer.h"
#include "index_writer.h"
#include <iostream>
#include <string.h>
#include <cassert>
#include <unistd.h>
#include <algorithm>

extern Parameters params;

Test::Test( int argc, char** argv )
{
    Filenames* fns = NULL;
    int testCount = 100000;
    
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
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
        else if ( !strcmp( argv[i], "-c" ) )
        {
            testCount = stoi( argv[++i] );
        }
    }
    
    ir_ = new IndexReader( fns );
    qb_ = new QueryBinaries( fns );
    
    srand( time(NULL) );
    int success = 0, failed = 0;
    double startTime = clock();
    for ( int i = 0; i < testCount; i++ )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = qb_->getSequence( id );
        bool rc = rand() % 2;
        vector<int> q( seq.size(), 0 );
        for ( int i = 0; i < seq.size(); i++ ) q[i] = rc ? 3 - charToInt[ seq[i] ] : charToInt[ seq.end()[-i-1] ];
        
        CharId rank, count;
        ir_->setBaseAll( q[0], q[1], rank, count );
        ( query( q, rank, count, 1 ) ? success : failed )++;
    }
    
    cout << "Tested a total of " << testCount << " reads as queries." << endl;
    if ( failed ) cout << success << " were found successfully, but " << failed << " were not." << endl;
    else cout << "All " << success << " were successfully found in the BWT." << endl;
    cout << "Total time taken: " << getDuration( startTime ) << endl;
}

bool Test::query( vector<int>& q, CharId rank, CharId count, int i )
{
    CharCount ranks, counts;
    ir_->countRange( q[i++], rank, count, ranks, counts );
    if ( i == q.size() ) return counts.endCounts;
    if ( counts[ q[i] ] ) return query( q, ranks[ q[i] ], counts[ q[i] ], i );
    assert( false );
    return false;
}


void Test::printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "    test" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "    Tests that the BWT has been correctly constructed by querying it with randomly-sampled reads." << endl;
    cout << endl << "Usage:" << endl;
    cout << "    locass index [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "    -p    Prefix for BWT data files." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "    -c    Number of reads to query (default: 10000)." << endl;
}
