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

#include "overlap.h"
#include "filenames.h"
#include "index_writer.h"
#include <iostream>
#include <string.h>
#include <unistd.h>

Overlap::Overlap( int argc, char** argv )
{
    string ifn, ofn, seq;
    int errors = 0;
    Filenames* fns = NULL;
    PreprocessFiles* ppf = NULL;
    
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
            
            ppf = new PreprocessFiles( prefix, true );
            fns = new Filenames( prefix );
        }
        else if ( !strcmp( argv[i], "-s" ) ) seq = argv[++i];
    }
    IndexWriter( ppf, 128, 20000 );
    IndexReader* ir = new IndexReader( fns );
    
    test( ir );
}

void Overlap::test( IndexReader* ir )
{
    CharId rank, edge, count, reads = 0;
    
    time_t my_time = time(NULL);
    cout << reads << "   " << ctime( &my_time );
    
    ir->setBaseAll( 1, 2, rank, edge, count );
    test( ir, rank, edge, count, reads, 2 );
}

void Overlap::test( IndexReader* ir, CharId rank, CharId edge, CharId count, CharId& reads, uint8_t c )
{
    CharCount ranks, edges, counts;
    ir->countRange( c, rank, edge, count, ranks, edges, counts );
    for ( int i = 0; i < 4; i++ ) if ( edges[i] )
    {
        if ( counts[i] ) test( ir, ranks[i], edges[i], counts[i], reads, i );
        else
        {
            if ( reads / 10000 < ( reads + edges[i] ) / 10000 )
            {
                time_t my_time = time(NULL);
                cout << reads << "   " << ctime( &my_time );
                if ( ( reads / 10000 ) >= 10 ) exit( 0 );
            }
            reads += edges[i];
        }
    }
}

void Overlap::printUsage()
{
    
}