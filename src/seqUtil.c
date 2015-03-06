#define __SEQUTIL_C__

#include "seqUtil.h"
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

typedef struct seqMem
{
	void *pos;
	struct seqMem *next;
} seqMem;

seqMem *start;
seqMem *current;

void seqMemInit()
{
	start = NULL;
}

void* seqMalloc(int size)
{
	if (start == NULL)
	{
		start = malloc(sizeof(seqMem));
		memset(start, '\0', sizeof(seqMem));
		current = start;
	}
	else
	{
		current->next = malloc(sizeof(seqMem));
		memset(current->next, '\0', sizeof(seqMem));
		current = current->next;
	}
	current->pos = malloc(size);
	memset(current->pos, '\0', size);
	return(current->pos);
}

void seqFreeAll()
{
	while (start != NULL)
	{
		current = start->next;
		free(start->pos);
		free(start);
		start = current;
	}

	start = NULL;
}

void seqFree(void* pos)
{
	seqMem *node, *prenode;
	node = start;
	prenode = start;
	while (node != NULL)
	{
		if (node->pos == pos)
		{
			free(node->pos);
			if (node == start)
			{
				start = node->next;
			}
			else if (node->next == NULL)
			{
				current = prenode;
				prenode->next = NULL;
			}
			else
			{
				prenode->next = node->next;
			}
			free(node);
			break;
		}

		prenode = node;
		node = node->next;
	}
}

void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen)
{
	int inputLen;
	char *temp;
	inputLen = size_of_string(input);
	if (inputLen == 0)
	{
		return;
	}
	while (*iMaxLen < (*iLen + inputLen) + 1)
	{
		*iMaxLen = *iMaxLen + APPEND_LEN;
	}
	temp = seqMalloc(*iMaxLen);
	if (*ppcStr == NULL)
	{
		memcpy(temp, input, inputLen);
	}
	else
	{
		memcpy(temp, *ppcStr, *iLen);
		strcat(temp, input);
	}
	*iLen = *iLen + inputLen;
	if (*ppcStr != NULL)
	{
		seqFree(*ppcStr);
	}
	*ppcStr = temp;
}

