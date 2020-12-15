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

#include "transform.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <cassert>
#include <algorithm>
//#include <chrono>
//#include <iomanip>

void Transform::load( PreprocessFiles* fns, vector< vector<ReadFile*> >& libs, uint8_t pairedLibCount, bool revComp )
{
    sort( libs.begin(), libs.end(), []( vector<ReadFile*> &a, vector<ReadFile*> &b ){
        return a.size() > b.size();
    } );
    
    // Set base read length
    uint8_t readLen = 0;
    for ( vector<ReadFile*> &lib : libs )
    {
        for ( ReadFile* readFile : lib )
        {
            readLen = max( readLen, readFile->readLen );
        }
    }
    uint8_t minLen = 45;
    assert( readLen >= 80 && readLen <= 255 );
    
    cout << "    Read length set to " << to_string( readLen ) << "." << endl;
    cout << "    Min length set to " << to_string( minLen ) << "." << endl;

    ofstream tmpSingles = fns->getWriteStream( fns->tmpSingles );
    ReadId readCount = 0, discardCount = 0;
    uint8_t fileCount = 0;
    double readStartTime = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
    
    BinaryWriter* binWrite = new BinaryWriter( fns, pairedLibCount, readLen, revComp );
    
    // Write binary sequence file and first transform cycle
    while ( !libs.empty() )
    {
        ReadId thisReadCount = 0;
        
        // Process paired libraries
        if ( libs[0].size() == 2 )
        {
            string reads[2];
            while ( libs[0][0]->getNext( reads[0] ) && libs[0][1]->getNext( reads[1] ) )
            {
                if ( reads[0].length() >= minLen && reads[1].length() >= minLen )
                {
                    binWrite->write( reads[0] );
                    binWrite->write( reads[1] );
                    thisReadCount += 2;
                }
                else
                {
                    discardCount += 2;
                    if ( reads[0].length() >= minLen )
                    {
                        tmpSingles << reads[0] << '\n';
                        --discardCount;
                    }
                    if ( reads[1].length() >= minLen )
                    {
                        tmpSingles << reads[1] << '\n';
                        --discardCount;
                    }
                }
            }
            
            binWrite->setNextLibrary();
            
            fileCount += libs[0][0] == libs[0][1] ? 1 : 2;
            delete libs[0][0];
            if ( libs[0][1] != libs[0][0] ) delete libs[0][1];
            
            if ( libs.size() == 1 || libs[1].size() == 1 )
            {
                fileCount += libs.size() - 1;
                tmpSingles.close();
                ReadFile* readFile = new ReadFile( fns->tmpSingles, readLen, 0 );
                vector<ReadFile*> lib = { readFile };
                libs.push_back( lib );
            }
            
            cout << "\tRead " << to_string( thisReadCount ) << " paired reads from library" << endl;
        }
        
        // Process singleton libraries
        else if ( libs[0].size() == 1 )
        {
            string read;
            while ( libs[0][0]->getNext( read ) )
            {
                if ( read.length() >= minLen )
                {
                    binWrite->write( read );
                }
                else
                {
                    discardCount++;
                }
                thisReadCount++;
            }
            delete libs[0][0];
            
            if ( thisReadCount )
            {
                cout << "\tRead " << to_string( thisReadCount ) << " single reads" << endl;
            }
        }
        
        libs.erase( libs.begin() );
        readCount += thisReadCount;
    }
    binWrite->close();
    
    cout << endl << "Reading inputs files... completed!" << endl << endl;
    cout << "Summary:" << endl;
    cout << "Read in " << to_string( readCount ) << " sequence reads and discarded " << to_string( discardCount ) << endl;
    cout << "Read from " << to_string( fileCount ) << " read files, including " << to_string( pairedLibCount ) << " paired libraries." << endl;
    cout << "Time taken: " << getDuration( readStartTime );
//    cout << "   " << std::fixed << std::setprecision(2) << ( clock() - readStartTime ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << endl << endl;
    cout << endl << endl;
    delete binWrite;
}

void Transform::run( PreprocessFiles* fns )
{
    cout << "Preprocessing step 2 of 3: transforming sequence data..." << endl << endl;
    
    BinaryReader* bin = new BinaryReader( fns );
    BwtCycler* cycler = new BwtCycler( fns );
    double totalStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
    
    while ( bin->cycle < bin->readLen )
    {
        double cycleStart = clock();
        cout << "    Cycle " << to_string( bin->cycle ) << " of " << to_string( bin->readLen ) << "... " << flush;
        
        bin->read();
        cycler->run( bin->chars, ( bin->anyEnds ? bin->ends : NULL ), bin->cycle );
        bin->update();
        
        cout << " completed in " << getDuration( cycleStart ) << endl;
    }
    
    double finalStart = clock();
    cout << "    Cycle " << to_string( bin->cycle ) << " of " << to_string( bin->readLen ) << "... " << flush;
    cycler->finish( bin->cycle + 1 );
    cout << " completed in " << getDuration( finalStart ) << endl;
    bin->finish();
    fns->clean();
    delete bin;
    delete cycler;
    
    cout << endl << "Transforming sequence data... completed!" << endl;
    cout << "Time taken: " << getDuration( totalStart );
//    cout << "   " << std::fixed << std::setprecision(2) << ( clock() - totalStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << endl << endl;
    cout << endl << endl;
}
