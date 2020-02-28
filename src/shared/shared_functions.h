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

#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include "types.h"

char getComp( char c );
int getEndTrim( string &q, string trim, bool drxn );
int getHomopolymerLen( string &s, bool drxn );
int getHomopolymerScore( string &s );
bool isSequence( string &s );
bool mapSeq( string &q, string &t, int* coords, int minLen );
int mapCongruence( string &left, string &right, int len );
bool mapSeqCongruent( string &left, string &right, int offset );
bool mapSeqEnd( string &q, string &t, int minLen, int32_t* coords, bool drxn );
int mapSeqOverlap( string &left, string &right, int minLen );
int mapSeqOverlap( string &q, string &t, int minLen, bool drxn );
void revComp( string &seq );
string revCompNew( string &seq );

#endif /* SHARED_FUNCTIONS_H */

