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

#ifndef LOCAL_ALIGNMENT_H
#define LOCAL_ALIGNMENT_H

#include <vector>
#include <string>

class LocalAlignment
{
public:
    LocalAlignment( std::string &a, std::string &b, bool glocal, bool freePolymer );
    
    static int isGapPoly( std::string (&a)[2], int d, int i, int gap );
    void print( int i, int j );
    void realign( std::string &a, std::string &b, bool conform, bool bluntStart=false, bool bluntEnd=false, bool trimEnd=false );
private:
    void blunt( bool start, bool end );
    int score( int i, int j, char &c, int* coords );
    void setRuns( int* run, int i, int j );
    
    std::string a_, b_;
    std::vector< std::vector<int> > m_;
    bool freePolymer_;
};



#endif /* LOCAL_ALIGNMENT_H */

