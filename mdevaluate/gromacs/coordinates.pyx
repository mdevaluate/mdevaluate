# distutils: sources = mdevaluate/gromacs/compression.c

from cython cimport view
from array import array
import numpy as np
cimport numpy as np


cdef extern from "compression.h":
    int xdrfile_decompress_coord_float(float     *coordinates,
                                   int       size,
                                   float     precision,
                                   int       minint[3],
                                   int       maxint[3],
                                   int       smallidx,
                                   char *    compressed_blob,
                                   size_t    blob_len,
                                   size_t *  readed_len)


def decompress(int size, float precision, minint, maxint, int smallidx, bytes blob):
     
    assert len(minint) == 3
    assert len(maxint) == 3
    
    cdef int minints[3]
    cdef int maxints[3]
    
    for i in range(3):
        minints[i] = minint[i]
        maxints[i] = maxint[i]
    
    
    cdef np.ndarray[float, ndim=2] coordinates = np.empty((size,3), dtype=np.float32)
    cdef size_t readed
    
    xdrfile_decompress_coord_float(
        <float *>coordinates.data, 
        size, 
        precision,
        maxints,
        minints, 
        smallidx, 
        blob,
        len(blob),
        &readed)
    return coordinates
