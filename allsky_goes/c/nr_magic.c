/* based on Numerical Recipes Utilities  */
/* providing byte-type and short int vectors and matrices       */
/* Ch. Reise, July 1993 , R. Mueller Apr. 2009      */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nr_magic.h"
#define NR_END 1
#define FREE_ARG char*

byte *bvector(long nl, long nh)
/* allocate a byte vector with subscript range v[nl..nh] */
{
    byte *v;

    v=(byte *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(byte)));
    if (!v) nrerror("allocation failure in bvector()");
    return v-nl+NR_END;
}

short int *svector(long nl, long nh)
/* allocate a byte vector with subscript range v[nl..nh] */
{
    short int *v;

    v=(short int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(short int)));
    if (!v) nrerror("allocation failure in svector()");
    return v-nl+NR_END;
}


byte **bmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a byte matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    byte **m;

    /* allocate pointers to rows */
    m=(byte **) malloc((size_t) ((nrow+NR_END)*sizeof(byte*)));
    if (!m) nrerror("allocation failure 1 in bmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(byte *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(byte)));
    if (!m[nrl]) nrerror("allocation failure 2 in bmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= nrl;

    for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

short int **smatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a byte matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    short int **m;

    /* allocate pointers to rows */
    m=(short int **) malloc((size_t) ((nrow+NR_END)*sizeof(short int*)));
    if (!m) nrerror("allocation failure 1 in bmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(short int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short int)));
    if (!m[nrl]) nrerror("allocation failure 2 in bmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= nrl;

    for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_bvector(byte *v, long nl, long nh)
/* free a byte vector allocated with bvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_svector(short int *v, long nl, long nh)
/* free a byte vector allocated with bvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_bmatrix(byte **m, long nrl, long nrh, long ncl, long nch)
/* free a byte matrix allocated by bmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

void free_smatrix(short int **m, long nrl, long nrh, long ncl, long nch)
/* free a byte matrix allocated by bmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
