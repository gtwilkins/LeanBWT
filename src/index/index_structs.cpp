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

#include "index_structs.h"

void CharCount::clear()
{
    counts[0] = counts[1] = counts[2] = counts[3] = 0;
    endCounts = 0;
}

int CharCount::getBranchCount()
{
    return bool(counts[0]) + bool(counts[1]) + bool(counts[2]) + bool(counts[3]);
}

int CharCount::getMaxBranch()
{
    int best = counts[0] > 0 ? 0 : ( counts[1] > 0 ? 1 : ( counts[2] > 0 ? 2 : ( counts[3] > 0 ? 3 : 4 ) ) );
    bool bad = false;
    for ( int i = best + 1; i < 4; i++ )
    {
        if ( !counts[i] ) continue;
        if ( counts[i] == counts[best] ) bad = true;
        else if ( counts[i] > counts[best] )
        {
            best = i;
            bad = false;
        }
    }
    return bad ? 4 : best;
}

void CharCount::operator -=( CharCount &rhs )
{
    for ( int i = 0; i < 4; i++ ) counts[i] -= rhs.counts[i];
    endCounts -= rhs.endCounts;
}