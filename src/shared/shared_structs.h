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

#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H

#include <vector>
#include "types.h"

struct Coords
{
    Coords( int32_t i, int32_t j ){ coords_[0] = i; coords_[1] = j; };
    int32_t& operator[]( int i ){ return coords_[i]; };
    int32_t coords_[2];
};

struct Read
{
    Read( string seq, ReadId id, int i, int j );
    static void sort( vector<Read>& reads, bool ascending, bool coordDrxn );
    string seq_;
    ReadId id_;
    Coords coords_;
};


#endif /* SHARED_STRUCTS_H */

