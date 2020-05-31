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

#include "index.h"
#include "transform_structs.h"
#include <iostream>
#include <sys/stat.h>
#include <string.h>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <unistd.h>

Index::Index( int argc, char** argv )
{
    ifstream infile;
    string prefix;
    PreprocessFiles* fns = NULL;
    bool isResume = false;
    bool didInput = false;
    bool doRevComp = true;
    int minScore = 0;
    
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
        }
        else if ( !strcmp( argv[i], "-i" ) )
        {
            if ( didInput )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            didInput = true;
            infile.open( argv[++i] );
            if ( !infile.good() )
            {
                cerr << "Failed to open file: \"" << argv[i] << "\"" << endl;
                exit( EXIT_FAILURE );
            }
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one output prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            prefix = argv[++i];
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
        }
        else if ( !strcmp( argv[i], "-s" ) ) minScore = stoi( argv[++i] );
        else if ( !strcmp( argv[i], "--resume" ) ) isResume = true;
        else if ( !strcmp( argv[i], "--no-rev-comp" ) ) doRevComp = false;
        else
        {
            cerr << "Unrecognised argument: \"" << argv[i] << "\"" << endl << endl;
            printUsage();
            exit( EXIT_FAILURE );
        }
    }
    
    if ( argc <= 2 )
    {
        cerr << "Error: no arguements supplied, see usage:" << endl << endl;
        printUsage();
        exit( EXIT_FAILURE );
    }
    
    double preprocessStartTime = clock();
    
    if ( prefix.empty() )
    {
        cerr << "Error: no output prefix supplied." << endl;
        exit( EXIT_FAILURE );
    }
    
    fns = new PreprocessFiles( prefix, true );
    
    if ( isResume && didInput )
    {
        cerr << "Error: resume (--resume) and input (-i) are mutually exclusive arguments." << endl;
        cerr << "If you wish to resume (--resume), simply specify the output prefix (-o) previously used." << endl;
        exit( EXIT_FAILURE );
    }
    else if ( didInput )
    {
        newTransform( fns, minScore, infile, doRevComp );
    }
    else if ( isResume )
    {
        resumeTransform( fns );
    }
    else
    {
        cerr << "Error: did not specify either an input (-i), or the resume flag (--resume)." << endl << endl;
        printUsage();
        exit( EXIT_FAILURE );
    }
    
    cout << "Preprocessing step 3 of 3: indexing transformed data..." << endl;
    IndexWriter idx( fns, 1024, 20000 );
    
    cout << endl << "Preprocessing completed!" << endl;
    cout << "Total time taken: " << getDuration( preprocessStartTime ) << endl;
}

void Index::newTransform( PreprocessFiles* fns, int minScore, ifstream &infile, bool revComp )
{
    uint8_t fileCount = 0, pairedLibCount = 0;
    
    if ( infile.is_open() && infile.good() )
    {
        vector< vector<ReadFile*> > libs;
        string line;
        int lineNum = 1;
        ReadFile* readFile = NULL;
        while( getline( infile, line ) )
        {
	    if ( line.empty() || line[0] == '#' ) continue;
            vector<string> args;
            
            size_t it = line.find( ' ' );
            while ( it != line.npos && args.size() < 4 )
            {
                args.push_back( line.substr( 0, it ) );
                line = line.substr( it + 1 );
                it = line.find( ' ' );
            }
            
            args.push_back( line );
            
            if ( ( args.size() == 2 || args.size() == 3 ) && args[0] == "paired" )
            {
                vector<ReadFile*> lib;
                readFile = new ReadFile( args[1], 0, minScore );
                fileCount++;
                lib.push_back( readFile );
                if ( args.size() == 3 )
                {
                    readFile = new ReadFile( args[2], 0, minScore );
                    fileCount++;
                }
                lib.push_back( readFile );
                libs.push_back( lib );
                pairedLibCount++;
            }
            else if ( args.size() == 2 && args[0] == "single" )
            {
                readFile = new ReadFile( args[1], 0, minScore );
                fileCount++;
                vector<ReadFile*> lib = { readFile };
                libs.push_back( lib );
            }
            else
            {
                cerr << "Error: invalid input form in line " << lineNum << " of input file." << endl << endl;
                cerr << "Expecting inputs to be in one of the following forms:" << endl;
                cerr << "\tpaired [separated_pair_file_1] [separated_pair_file_2]" << endl;
                cerr << "\tpaired [interleaved_pair_file]" << endl;
                cerr << "\tsingle [unpaired_file]" << endl;
                exit( EXIT_FAILURE );
            }
            
            lineNum++;
        }
        
        if ( pairedLibCount > 5 )
        {
            cerr << "Error: Excessive library count of " << pairedLibCount << ". Maximum supported is 5." << endl;
            exit( EXIT_FAILURE );
        }
        
        cout << "Preprocessing step 1 of 3: reading input files..." << endl << endl;
        Transform::load( fns, libs, pairedLibCount, revComp );
        Transform::run( fns, revComp );
    }
    else
    {
        cerr << "Error: No input file provided." << endl;
        exit( EXIT_FAILURE );
    }
}

void Index::resumeTransform( PreprocessFiles* fns )
{
    cout << "Resuming preprocessing..." << endl << endl;
    Transform::run( fns, false );
}

void Index::printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "\tindex" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "\tTransforms and indexes sequence read data in preparation for assembly." << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tlocass index [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "\t-i\tInput text file containing a list of sequence read files. See notes for details." << endl;
    cout << "\t-p\tOutput prefix for transformed sequence files." << endl;
    cout << endl << "Notes:" << endl;
    cout << "\t- Accepted read file formats are fasta, fastq or a list of sequences, one per line." << endl;
    cout << "\t- Input read libraries can be either paired or single." << endl;
    cout << "\t- Each paired read library can be input as either two separated files or one interleaved file." << endl;
    cout << "\t- Each line of the input text file is expected in one of the following forms:" << endl;
    cout << "\t\tpaired [separated_pair_file_1] [separated_pair_file_2]" << endl;
    cout << "\t\tpaired [interleaved_pair_file]" << endl;
    cout << "\t\tsingle [unpaired_file]" << endl;
    cout << "\t- At least one paired read library is required." << endl;
}
