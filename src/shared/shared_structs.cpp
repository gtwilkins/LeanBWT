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

#include "shared_structs.h"
#include <algorithm>

Read::Read( std::string seq, ReadId id, int i, int j )
: seq_( seq ), id_( id ), coords_( Coords( i, j ) )
{
    
}

void Read::sort( vector<Read>& reads, bool ascending, bool coordDrxn )
{
    if ( ascending ) std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] < b.coords_[coordDrxn]; } );
    else std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] > b.coords_[coordDrxn]; } );
}