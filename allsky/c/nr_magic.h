#ifndef _NR_MAGIC_H_
#define _NR_MAGIC_H_

typedef unsigned char byte;

byte *bvector(long nl, long nh);
byte **bmatrix(long nrl, long nrh, long ncl, long nch);
short int *svector(long nl, long nh);
short int **smatrix(long nrl, long nrh, long ncl, long nch);

void free_bvector(byte *v, long nl, long nh);
void free_bmatrix(byte **m, long nrl, long nrh, long ncl, long nch);
void free_svector(short int *v, long nl, long nh);
void free_smatrix(short int **m, long nrl, long nrh, long ncl, long nch);


#endif 
