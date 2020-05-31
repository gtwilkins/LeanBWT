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

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "types.h"
#include "transform.h"
#include "index_writer.h"

class Index
{
public:
    Index( int argc, char** argv );
    
    void newTransform( PreprocessFiles* fns, int minScore, ifstream &infile, bool revComp );
    void resumeTransform( PreprocessFiles* fns );
    
    void printUsage();
private:
};

#endif /* PREPROCESS_H */

