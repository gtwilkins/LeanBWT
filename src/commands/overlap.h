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

#ifndef OVERLAP_H
#define OVERLAP_H

#include "types.h"
#include "index_reader.h"

class Overlap
{
public:
    Overlap( int argc, char** argv );
    
private:
    void overlap( string& seq, ofstream* ofs, int errors );
    
    void test( IndexReader* ir );
    void test( IndexReader* ir, CharId rank, CharId edge, CharId count, CharId& reads, uint8_t c );
    void printUsage();
};

#endif /* OVERLAP_H */

