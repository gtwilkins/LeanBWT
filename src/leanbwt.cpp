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

#include <cstdlib>
#include <cassert>
#include <string.h>
#include <iostream>
#include "shared/types.h"
#include "commands/index.h"
#include "shared_functions.h"
#include "parameters.h"

Parameters params;

void printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tleanbwt <command> [args]" << endl;
    cout << endl << "Accepted commands:" << endl;
    cout << "\tindex\t\tTransforms and indexes raw sequence data." << endl;
    cout << endl << "To view the usage and accepted arguments for a given command type:" << endl;
    cout << "\tleanbwt <command> -h" << endl;
}

int main( int argc, char** argv )
{
    if ( argc > 1 )
    {
        if ( !strcmp( argv[1], "-h" ) || !strcmp( argv[1], "--help" ) || !strcmp( argv[1], "-help" ) )
        {
            printUsage();
        }
        else if ( !strcmp( argv[1], "index" ) )
        {
            Index( argc, argv );
        }
//        else if ( !strcmp( argv[1], "match" ) )
//        {
//            Match( argc, argv );
//        }
//        else if ( !strcmp( argv[1], "overlap" ) )
//        {
//            Overlap( argc, argv );
//        }
        else
        {
            cerr << "Unrecognised command: \"" << argv[1] << "\"" << endl << endl;
            printUsage();
        }
    }
    else
    {
        printUsage();
    }
    return 0;
}

