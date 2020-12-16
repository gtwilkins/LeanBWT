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

#include "shared_functions.h"
#include <algorithm>
#include <cassert>
#include <string.h>

char getComp( char c )
{
    if ( c == 'A' ) return 'T';
    if ( c == 'C' ) return 'G';
    if ( c == 'G' ) return 'C';
    if ( c == 'T' ) return 'A';
    if ( c == 'N' ) return 'N';
    assert( false );
    return 'N';
}

int getEndTrim( string &q, string trim, bool drxn )
{
    if ( trim.empty() ) return 0;
    int len = min( q.length(), trim.length() );
    if ( drxn )
    {
        for ( int i = q.length() - len; i < q.length(); i++ )
        {
            int j = 0;
            while ( q[i+j] == trim[j] && j < len ) j++;
            if ( i + j == q.length() ) return j;
        }
    }
    else
    {
        assert( false );
        for ( int i = len; --i >= 0; )
        {
            int j = 0;
            while ( q[i-j] == trim[ trim.length() - j ] && i >= j ) j++;
            if ( j > i ) return j;
        }
    }
    
    return 0;
}

int getHomopolymerLen( string &s, bool drxn )
{
    int i = 0;
    while ( i+1 < s.length() &&
            ( drxn ? s[s.length()-i-1] == s[s.length()-i-2] : s[i] == s[i+1] ) ) i++;
    if ( i+2 < s.length() && ( drxn ? s[s.length()-i-2] == s[s.length()-i-3] : s[i+1] == s[i+2] ) )
    {
        int j = ++i;
        while ( i+1 < s.length() &&
                ( drxn ? s[s.length()-j-1] == s[s.length()-i-2] : s[j] == s[i+1] ) ) i++;
    }
    if ( i ) i += 2;
    return i;
}

int getHomopolymerScore( string &s )
{
    int score = 0;
    int i = 0, j = 0;
    while ( j++ < s.length() )
    {
        if ( s[i] != s[j] )
        {
            if ( j > i + 2 )
            {
                score += j - i - 2;
            }
            i = j;
        }
    }
    
    score += j - i - 1;
    return score;
}

bool getSeq( ifstream& ifs, string& header, string& seq )
{
    string line;
    if ( !getline( ifs, line ) || line.empty() ) return false;
    if ( line[0] == '>' )
    {
        header = line.substr( 1 );
        if ( !getline( ifs, seq ) || !isSequence( seq ) ) return false;
    }
    else if ( line[0] == '@' )
    {
        header = line.substr( 1 );
        if ( !getline( ifs, seq ) || !isSequence( seq ) ) return false;
        assert( getline( ifs, line ) && line[0] == '+' );
        assert( getline( ifs, line ) && line.size() == seq.size() );
    }
    else if ( !isSequence( line ) ) return false;
    else seq = line;
    
    return true;
}

bool isSequence( string &s )
{
    for ( char c : s )  if ( !strchr( "ACGTN", c ) ) return false;
    return true;
}

bool mapSeq( string &q, string &t, int* coords, int minLen )
{
    if ( q.find( t ) != string::npos )
    {
        coords[0] = 0;
        coords[1] = t.size();
        return true;
    }
    size_t it = t.find( q );
    if ( it != string::npos )
    {
        coords[0] = it;
        coords[1] = coords[0] + q.size();
        return true;
    }
    if ( minLen > q.size() ) return false;
    int ols[2] = { mapSeqOverlap( q, t, minLen ), mapSeqOverlap( t, q, minLen ) };
    if ( ols[0] && ols[0] >= ols[1] )
    {
        coords[0] = 0;
        coords[1] = ols[0];
        return true;
    }
    if ( ols[1] )
    {
        coords[0] = t.size() - ols[1];
        coords[1] = t.size();
        return true;
    }
    return false;
}

int mapCongruence( string &left, string &right, int len )
{
    if ( len > min( left.size(), right.size() ) ) return 0;
    if ( left.size() < right.size() )
    {
        size_t i = left.size() - len;
        size_t j = right.find( left.substr( i ) );
        if ( j == right.npos ) return 0;
        while ( i > 0 && j > 0 )
        {
            if ( left[--i] != right[--j] ) return 0;
            len++;
        }
    }
    else
    {
        size_t i = left.find( right.substr( 0, len ) );
        size_t j = len;
        if ( i == left.npos ) return 0;
        i += len;
        while ( i < left.size() && j < right.size() )
        {
            if ( left[i++] != right[j++] ) return 0;
            len++;
        }
    }
    
    return len;
}

bool mapSeqCongruent( string &left, string &right, int offset )
{
    int i = 0;
    while ( i + offset < left.length() && i < right.length() && left[i+offset] == right[i++] );
    return i + offset == left.length();
}

bool mapSeqEnd( string &q, string &t, int minLen, int32_t* coords, bool drxn )
{
    if ( minLen > q.length() ) return false;
    string qSub = drxn ? q.substr( q.length() - minLen ) : q.substr( 0, minLen );
    size_t it = t.find( qSub );
    bool didMap = false;
    while ( it != t.npos )
    {
        didMap = true;
        coords[0] = it;
        coords[1] = it + qSub.length();
        int i = drxn ? q.length() - qSub.length(): qSub.length();
        if ( drxn )
        {
            while ( coords[0]-1 >= 0 && i-1 >= 0 && q[i-1] == t[coords[0]-1] )
            {
                coords[0]--;
                i--;
            }
        }
        else
        {
            while ( i < q.length() && coords[1] < t.length() && q[i] == t[coords[1]] )
            {
                coords[1]++;
                i++;
            }
        }
        minLen = coords[1] - coords[0] + 1;
        if ( minLen > q.length() ) break;
        qSub = drxn ? q.substr( q.length() - minLen ) : q.substr( 0, minLen );
        it = t.find( qSub );
    }
    
    return didMap;
}

int mapSeqOverlap( string &left, string &right, int minLen )
{
    if ( minLen <= right.length() )
    {
        string q = right.substr( 0, minLen );
        size_t it = left.find( q );
        while ( it != left.npos )
        {
            int i = it + minLen;
            int j = minLen;
            while ( i < left.length() && j < right.length() && left[i] == right[j] )
            {
                i++;
                j++;
            }
            if ( i == left.length() ) return j;
            it = left.find( q, it + 1 );
        }
    }
    return 0;
}

int mapSeqOverlap( string &q, string &t, int minLen, bool drxn )
{
    return mapSeqOverlap( drxn ? q : t, drxn ? t : q, minLen );
}

void revComp( string &seq )
{
    reverse( seq.begin(), seq.end() );
    for ( int i = 0; i < seq.length(); i++ )
    {
        if ( seq[i] == 'A' ) seq[i] = 'T';
        else if ( seq[i] == 'T' ) seq[i] = 'A';
        else if ( seq[i] == 'C' ) seq[i] = 'G';
        else if ( seq[i] == 'G' ) seq[i] = 'C';
        else if ( seq[i] == 'a' ) seq[i] = 't';
        else if ( seq[i] == 't' ) seq[i] = 'a';
        else if ( seq[i] == 'c' ) seq[i] = 'g';
        else if ( seq[i] == 'g' ) seq[i] = 'c';
        else seq[i] = 'N';
    }
    
}

string revCompNew( string &seq )
{
    string rev;
    for ( auto it = seq.rbegin(); it != seq.rend(); it++ )
    {
        if ( *it == 'A' ) rev += 'T';
        else if ( *it == 'T' ) rev += 'A';
        else if ( *it == 'C' ) rev += 'G';
        else if ( *it == 'G' ) rev += 'C';
        else if ( *it == 'a' ) rev += 't';
        else if ( *it == 't' ) rev += 'a';
        else if ( *it == 'c' ) rev += 'g';
        else if ( *it == 'g' ) rev += 'c';
        else rev += 'N';
    }
    return rev;
}
