#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <limits.h>
#include <time.h>


#define PAUSE fgets(stp,sizeof(stp),stdin) // shortcut for pausing

#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

#define NEW(p,type)     if ((p=(type *) malloc (sizeof(type))) == NULL) {\
                                printf ("Out of Memory!\n");\
                                exit(0);\
                        }

#define FREE(p)         if (p) { free ((char *) p); p = NULL; }


#define ADD( head, p )  if ( head )  { \
                                p->next = head; \
                                p->prev = head->prev; \
                                head->prev = p; \
                                p->prev->next = p; \
                        } \
                        else { \
                                head = p; \
                                head->next = head->prev = p; \
                        }

#endif
