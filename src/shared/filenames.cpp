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

#include "filenames.h"
#include <iostream>
#include <cassert>
#include <sys/stat.h>

Filenames::Filenames( string inPrefix )
: prefix( inPrefix )
{
    string folder = inPrefix.substr( 0, inPrefix.find_last_of( '/' ) );
    makeFolder( folder );
    
    bin = prefix + "-bin.dat";
    bwt = prefix + "-bwt.dat";
    ids = prefix + "-ids.dat";
    idx = prefix + "-idx.dat";
    mer = prefix + "-mer.dat";
}

bool Filenames::exists( string &filename )
{
    if ( filename.empty() ) return false;
    ifstream ifs( filename );
    bool result = ifs.good();
    ifs.close();
    return result;
}

FILE* Filenames::getReadPointer( string &filename, bool doEdit, bool allowFail )
{
    FILE* fp = fopen( filename.c_str(), ( doEdit ? "rb+" : "rb" ) );
    if ( fp == NULL && !allowFail )
    {
        cerr << "Error opening file \"" << filename << "\"." << endl;
        exit( EXIT_FAILURE );
    }
    return fp;
}

FILE* Filenames::getWritePointer( string &filename )
{
    FILE* fp = fopen( filename.c_str(), "wb" );
    if ( fp == NULL )
    {
        cerr << "Error opening file \"" << filename << "\"." << endl;
        exit( EXIT_FAILURE );
    }
    return fp;
}

ifstream Filenames::getReadStream( string &filename )
{
    ifstream fp( filename );
    
    if ( !fp.good() || !fp.is_open() )
    {
        cerr << "Error opening file \"" << filename << "\"" << endl;
        exit( EXIT_FAILURE );
    }
    
    return fp;
}

ofstream Filenames::getWriteStream( string &filename )
{
    ofstream fp( filename );
    if ( !fp.good() || !fp.is_open() )
    {
        cerr << "Error creating file \"" << filename << "\"" << endl;
        exit( EXIT_FAILURE );
    }
    return fp;
}

FILE* Filenames::getBinary( bool doRead, bool doEdit )
{
    return ( doRead ? getReadPointer( bin, doEdit ) : getWritePointer( bin )  );
}

bool Filenames::isFolder( string folder )
{
    struct stat st;
    return ( stat( folder.c_str(), &st ) == 0 && ( S_ISDIR( st.st_mode ) ) );
}

void Filenames::makeFolder( string folder )
{
    if ( isFolder( folder ) ) return;
    size_t it = folder.find_last_of( '/' );
    if ( it != folder.npos )
    {
        string parent = folder.substr( 0, it );
        if ( !isFolder( parent ) ) makeFolder( parent );
    }
    
    const int dir_err = mkdir( folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    if ( dir_err == -1 )
    {
        cerr << "Error creating directory \"" << folder << "\"!" << endl;
        exit(1);
    }
}

void Filenames::removeFile( string &filename )
{
    if ( remove( filename.c_str() ) )
    {
        cerr << "Warning: could not remove file \"" << filename << "\"" << endl;
    }
}

void Filenames::setIndex( FILE* &inBin, FILE* &inBwt, FILE* &inIdx, FILE* &inMer )
{
    inBin = getReadPointer( bin, false );
    inBwt = getReadPointer( bwt, false );
    inIdx = getReadPointer( idx, false );
    inMer = getReadPointer( mer, false, true );
}

PreprocessFiles::PreprocessFiles( string inPrefix, bool overwrite )
: Filenames( inPrefix )
{
    tmpSingles = prefix + "-tmpSingles.seq";
    for ( int i( 0 ); i < 2; i++ )
    {
        tmpBwt[i] = prefix + "-bwt-tmp" + to_string( i + 1 );
        tmpEnd[i] = prefix + "-end-tmp" + to_string( i + 1 );
        for ( int j( 0 ); j < 4; j++ )
        {
            tmpIns[i][j] = prefix + "-ins-" + to_string( j + 1 ) + "-tmp" + to_string( i + 1 );
            for ( int k( 0 ); k < 5; k++ )
            {
                tmpIds[i][j][k] = prefix + "-ids-" + to_string( j + 1 ) + to_string( k + 1 ) + "-tmp" + to_string( i + 1 );
            }
        }
    }
    
    for ( string const &fn : { bwt, bin, ids, idx, mer } )
    {
        if ( ifstream( fn ) && !overwrite )
        {
            cerr << "Error: file \"" << fn << "\" already exists. Either remove it, rename it, or choose a different file prefix." << endl;
            exit( EXIT_FAILURE );
        }
    }
}

void PreprocessFiles::clean()
{
    removeFile( tmpSingles );
    for ( int i( 0 ); i < 2; i++ )
    {
        removeFile( tmpBwt[i] );
        removeFile( tmpEnd[i] );
        for ( int j( 0 ); j < 4; j++ )
        {
            removeFile( tmpIns[i][j] );
            for ( int k( 0 ); k < 5; k++ )
            {
                removeFile( tmpIds[i][j][k] );
            }
        }
    }
}

void PreprocessFiles::setBinaryWrite( FILE* &outBin, FILE* &outBwt, FILE* &outEnd, FILE* (&outIns)[4], FILE* (&outIds)[4][4] )
{
    outBin = getWritePointer( bin );
    outBwt = getWritePointer( tmpBwt[0] );
    outEnd = getWritePointer( tmpEnd[0] );
    for ( int i ( 0 ); i < 4; i++ )
    {
        outIns[i] = getWritePointer( tmpIns[0][i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            outIds[i][j] = getWritePointer( tmpIds[0][i][j] );
        }
    }
}

void PreprocessFiles::setCycler( FILE* &inBwt, FILE* &outBwt, FILE* &inEnd, FILE* &outEnd, FILE* (&outIns)[4], FILE* (&outIds)[4][5], uint8_t cycle )
{
    uint8_t iIn = cycle & 1;
    
    inBwt = getReadPointer( tmpBwt[iIn], false );
    outBwt = getWritePointer( tmpBwt[!iIn] );
    inEnd = getReadPointer( tmpEnd[iIn], false );
    outEnd = getWritePointer( tmpEnd[!iIn] );
    for ( int i ( 0 ); i < 4; i++ )
    {
        outIns[i] = getWritePointer( tmpIns[!iIn][i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            outIds[i][j] = getWritePointer( tmpIds[!iIn][i][j] );
        }
    }
}

void PreprocessFiles::setCyclerIter( FILE* &inIns, FILE* (&inIds)[5], uint8_t cycle, uint8_t i )
{
    uint8_t iIn = cycle & 1;
    
    inIns = getReadPointer( tmpIns[iIn][i], false );
    for ( int j ( 0 ); j < 5; j++ )
    {
        inIds[j] = getReadPointer( tmpIds[iIn][i][j], false );
    }
}

void PreprocessFiles::setCyclerFinal( FILE* &inBwt, FILE* &outBwt, FILE* &inEnd, FILE* &outEnd, uint8_t cycle )
{
    uint8_t iIn = cycle & 1;
    
    inBwt = getReadPointer( tmpBwt[iIn], false );
    outBwt = getWritePointer( bwt );
    inEnd = getReadPointer( tmpEnd[iIn], false );
    outEnd = getWritePointer( ids );
}

void PreprocessFiles::setCyclerFinalIter( FILE* &inIns, FILE* &inIds, uint8_t cycle, uint8_t i )
{
    uint8_t iIn = cycle & 1;
    
    inIns = getReadPointer( tmpIns[iIn][i], false );
    inIds = getReadPointer( tmpIds[iIn][i][4], false );
}

void PreprocessFiles::setCyclerUpdate( FILE* &outBwt, FILE* &outEnd, FILE* (&outIns)[4], FILE* (&outIds)[4][5], uint8_t cycle )
{
    uint8_t iIn = cycle & 1;
    
    outBwt = getReadPointer( tmpBwt[!iIn], true );
    outEnd = getReadPointer( tmpEnd[!iIn], true );
    for ( int i ( 0 ); i < 4; i++ )
    {
        outIns[i] = getReadPointer( tmpIns[!iIn][i], true );
        for ( int j ( 0 ); j < 5; j++ )
        {
            outIds[i][j] = getReadPointer( tmpIds[!iIn][i][j], true );
        }
    }
}

void PreprocessFiles::setIndexWrite( FILE* &inBwt, FILE* &outIdx )
{
    inBwt = getReadPointer( bwt, false );
    outIdx = getWritePointer( idx );
}

void PreprocessFiles::setMersWrite( FILE* &outMer )
{
    outMer = getWritePointer( mer );
}
