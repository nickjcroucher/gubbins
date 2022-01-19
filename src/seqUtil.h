#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "string_cat.h"

/*
 *  Yu-Wei Wu  http://yuweibioinfo.blogspot.com/2008/10/newick-tree-parser-in-c-make-use-of.html
 *  Copyright (C) 2011  Yu-Wei Wu
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#define APPEND_LEN	256

#ifdef __SEQUTIL_C__
    void seqMemInit(void);
	void* seqMalloc(int size);
    void seqFree(void*);
    void seqFreeAll(void);
	void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
#else
    extern void seqMemInit(void);
	extern void* seqMalloc(int size);
    extern void seqFreeAll(void);
    extern void seqFree(void*);
	extern void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
#endif



#endif

