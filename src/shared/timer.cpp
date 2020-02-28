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

#include "timer.h"
#include <ctime>

string getDuration( double startTime )
{
    double totalDuration = ( clock() - startTime ) / CLOCKS_PER_SEC;
    int hrs = totalDuration / 3600;
    int mins = int( totalDuration / 60 ) % 60;
    int secs = int(max( 1.0, totalDuration ) ) % 60;
    return ( ( hrs > 0 ? to_string( hrs ) + " hr " : "" )
            + ( mins > 0 ? to_string( mins ) + " min " : "" )
            + to_string( secs ) + " sec" );
}
