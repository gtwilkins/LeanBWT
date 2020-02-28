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

#include "transform_structs.h"
#include <cassert>
#include <iostream>

ReadFile::ReadFile( string filename, int baseReadLen, int minScore )
: readLen( baseReadLen ), minPhred( minScore )
{
    fh.open( filename );
    
    if ( !fh.good() || !fh.is_open() )
    {
        cerr << "Error: could not open file \"" << filename << "\"" << endl;
        exit( EXIT_FAILURE );
    }
    
    string line;
    if ( getline( fh, line ) && !line.empty() )
    {
        if ( line[0] == '>' )
        {
            fileType = 1;
        }
        else if ( line[0] == '@' )
        {
            fileType = 3;
        }
        else if ( charToInt[ line[0] ] < 6 )
        {
            fileType = 0;
            fh.seekg( 0 );
        }
        else
        {
            cerr << "Error: Unrecognised sequence file type." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    setReadLen();
}

bool ReadFile::getNext( string &seq )
{
    if ( getline( fh, seq ) )
    {
        if ( fileType ) getline( fh, line );
        if ( fileType == 3 )
        {
            getline( fh, line );
            for ( int i( 0 ); i < line.size(); i++ )
            {
                if ( line[i] < minPhred ) seq[i] = 'N';
            }
            getline( fh, line );
        }
        
        trimSeq( seq );
        
        return true;
    }
    
    return false;
}

void ReadFile::setReadLen()
{
    string seq;
    size_t curr = fh.tellg();
    fh.seekg( 0 );
    int i = 0, j = 0;
    while ( getline( fh, seq ) && !seq.empty() && i < 1000 )
    {
        if ( !fileType || j == 1 )
        {
            if ( seq.length() > 255 )
            {
                cerr << "Error: Read length of " << seq.length() << " detected. Maximum length of 255 is supported." << endl;
                exit( EXIT_FAILURE );
            }
            readLen = max( readLen, (uint8_t)seq.length() );
        }
        if ( j++ == fileType )
        {
            j = 0;
            i++;
        }
    }
    
    if ( readLen < 80 )
    {
        cerr << "Error: Read length of " << seq.length() << " detected. Minimum length of 80 is supported." << endl;
        exit( EXIT_FAILURE );
    }
    fh.seekg( curr );
    fh.clear();
}

void ReadFile::trimSeq( string &seq )
{
    int iBest = 0;
    int bestLen = -1;
    int iCurr = 0;
    for ( int i( 0 ); i < seq.length(); i++ )
    {
        if ( charToInt[ seq[i] ] >= 4 )
        {
            if ( charToInt[ seq[i] ] == 6 )
            {
                cerr << "Error: Unrecognised character \"" << seq[i] << "\" in sequence file." << endl;
                exit( EXIT_FAILURE );
            }
            int len = i - iCurr;
            if ( len > bestLen )
            {
                iBest = iCurr;
                bestLen = len;
            }
            iCurr = i + 1;
        }
    }
    
    if ( bestLen > -1 )
    {
        int currLen = seq.length() - iCurr;
        if ( currLen > bestLen )
        {
            iBest = iCurr;
            bestLen = currLen;
        }
        seq = seq.substr( iBest, bestLen );
    }
}
