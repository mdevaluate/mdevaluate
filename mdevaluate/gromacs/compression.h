#include <stdio.h>

int
xdrfile_decompress_coord_float(float     *coordinates,
                               int       size,
                               float     precision,
                               int       minint[3],
                               int       maxint[3],
                               int       smallidx,
                               int *     compressed_blob,
                               size_t    blob_len,
                               size_t *  readed_len);
