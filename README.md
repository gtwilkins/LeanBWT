# LeanBWT

## Overview
LeanBWT is an application for the assembly memory-efficient Burrows-Wheeler Transform and FM-index construction.

## Limitations
LeanBWT is designed for use with short read sequence data, namely those produced by Illumina

## Requirements
* gcc

## Installation
The install directory can be specified with the following command (if this omitted, LeanBWT is installed to /usr/local/bin/):

	./configure --prefix=/path/to/directory/

LeanBWT can be installed with the following commands:

	make
	sudo make install

## Use
leanbwt -h
