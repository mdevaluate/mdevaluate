#include <stdio.h>
#include <memory.h>
#include <stdlib.h>


extern void test(int *b) {
	b[0] = 1;
	b[1] = 0;
	b[2] = 2;
}

static const int magicints[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, 
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 
};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

/* Internal support routines for reading/writing compressed coordinates
 * sizeofint - calculate smallest number of bits necessary
 * to represent a certain integer.
 */
static int
sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32)
    {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

/*
 * sizeofints - calculate 'bitsize' of compressed ints
 *
 * given a number of small unsigned integers and the maximum value
 * return the number of bits needed to read or write them with the
 * routines encodeints/decodeints. You need this parameter when
 * calling those routines.
 * (However, in some cases we can just use the variable 'smallidx'
 * which is the exact number of bits, and them we dont need to call
 * this routine).
 */
static int
sizeofints(int num_of_ints, unsigned int sizes[])
{
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++)
    {
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num)
    {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
    
}

/*
 * decodebits - decode number from buf using specified number of bits
 *
 * extract the number of bits from the array buf and construct an integer
 * from it. Return that value.
 *
 */

static int
decodebits(int state[], int _buf[], int num_of_bits)
{
    
    int cnt, num;
    unsigned int lastbits, lastbyte;
    int mask = (1 << num_of_bits) -1;
    unsigned char * cbuf = (unsigned char *)_buf;
    cnt = state[0];
    lastbits = (unsigned int) state[1];
    lastbyte = (unsigned int) state[2];
    
    num = 0;
    while (num_of_bits >= 8)
    {
        lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
        num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -=8;
    }
    if (num_of_bits > 0)
    {
        if (lastbits < num_of_bits)
        {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    state[0] = cnt;
    state[1] = lastbits;
    state[2] = lastbyte;
    return num;
}

/*
 * decodeints - decode 'small' integers from the buf array
 *
 * this routine is the inverse from encodeints() and decodes the small integers
 * written to buf by calculating the remainder and doing divisions with
 * the given sizes[]. You need to specify the total number of bits to be
 * used from buf in num_of_bits.
 *
 */

static void
decodeints(int state[], int buf[], int num_of_ints, int num_of_bits,
           unsigned int sizes[], int nums[])
{
    
    int bytes[32];
    int i, j, num_of_bytes, p, num;
    
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8)
    {
        bytes[num_of_bytes++] = decodebits(state, buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        bytes[num_of_bytes++] = decodebits(state, buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--)
    {
        num = 0;
        for (j = num_of_bytes-1; j >=0; j--)
        {
            num = (num << 8) | bytes[j];
            p = num / sizes[i];
            bytes[j] = p;
            num = num - p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

extern int
xdrfile_decompress_coord_float(float     *coordinates,
                               int       size,
                               float     precision,
                               int       minint[3],
                               int       maxint[3],
                               int       smallidx,
                               char*     compressed_blob_,
                               size_t    blob_len,
                               size_t *     readed_len)
{
    int *compressed_blob = (int*)compressed_blob_;
    int *lip;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
    int k, *buf1 , flag;
    int smallnum, smaller, i, is_smaller, run;
    float *lfp, inv_precision;
    int tmp, *thiscoord,  prevcoord[3];
    unsigned int bitsize;
    
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;
    
    
    size3 = size * 3;
    
    if((buf1=(int *)malloc(sizeof(int)*size3))==NULL)
    {
        fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
        return -1;
    }
    
    /* Dont bother with compression for three atoms or less */
    if(size<=9)
    {
        
        // TODO...
        /* return number of coords, not floats */
    }
    
    
    sizeint[0] = maxint[0] - minint[0]+1;
    sizeint[1] = maxint[1] - minint[1]+1;
    sizeint[2] = maxint[2] - minint[2]+1;
    
    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
    {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        bitsize = 0; /* flag the use of large sizes */
    }
    else
    {
        bitsize = sizeofints(3, sizeint);
    }
    
    
    tmp = smallidx+8;
    tmp = smallidx-1;
    tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
    smaller = magicints[tmp] / 2;
    smallnum = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    
    
    lfp = coordinates;
    inv_precision = 1.0 / precision;
    run = 0;
    i = 0;
    lip = buf1;
    
    int state[3] = {0,0,0};
    
    while ( i < size )
    {
        thiscoord = (int *)(lip) + i * 3;
        
        if (bitsize == 0)
        {
            thiscoord[0] = decodebits(state, compressed_blob, bitsizeint[0]);
            thiscoord[1] = decodebits(state, compressed_blob, bitsizeint[1]);
            thiscoord[2] = decodebits(state, compressed_blob, bitsizeint[2]);
        }
        else
        {
            decodeints(state, compressed_blob, 3, bitsize, sizeint, thiscoord);
        }
        
        i++;
        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];
        
        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];
        
        flag = decodebits(state, compressed_blob, 1);
        is_smaller = 0;
        if (flag == 1)
        {
            run = decodebits(state, compressed_blob, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }
        if (run > 0)
        {
            thiscoord += 3;
            for (k = 0; k < run; k+=3)
            {
                decodeints(state, compressed_blob, 3, smallidx, sizesmall, thiscoord);
                i++;
                thiscoord[0] += prevcoord[0] - smallnum;
                thiscoord[1] += prevcoord[1] - smallnum;
                thiscoord[2] += prevcoord[2] - smallnum;
                if (k == 0) {
                    /* interchange first with second atom for better
                     * compression of water molecules
                     */
                    tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
                    prevcoord[0] = tmp;
                    tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
                    prevcoord[1] = tmp;
                    tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
                    prevcoord[2] = tmp;
                    *lfp++ = prevcoord[0] * inv_precision;
                    *lfp++ = prevcoord[1] * inv_precision;
                    *lfp++ = prevcoord[2] * inv_precision;
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
        } 
        else
        {
            *lfp++ = thiscoord[0] * inv_precision;
            *lfp++ = thiscoord[1] * inv_precision;
            *lfp++ = thiscoord[2] * inv_precision;		
        }
        smallidx += is_smaller;
        if (is_smaller < 0) 
        {
            smallnum = smaller;
            
            if (smallidx > FIRSTIDX) 
            {
                smaller = magicints[smallidx - 1] /2;
            } 
            else 
            {
                smaller = 0;
            }
        } 
        else if (is_smaller > 0)
        {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    }
    
    free((void*)buf1);
    *readed_len = state[0];
    return 0;
}
