from .reader import XTCReader, index_filename_for_xtc, INDEX_MAGIC
from xdrlib import Packer
import os
import sys
import logging


def index_xtcfile(filename):

    index_filename = index_filename_for_xtc(filename)
    packer = Packer()
    xtc_stat = os.stat(filename)
    c_time = int(xtc_stat.st_ctime)
    m_time = int(xtc_stat.st_mtime)
    size = xtc_stat.st_size

    with XTCReader(filename, load_cache_file=False) as reader, open(index_filename, 'wb') as idx_fd:
        packer.pack_hyper(INDEX_MAGIC)
        packer.pack_hyper(c_time)
        packer.pack_hyper(m_time)
        packer.pack_hyper(size)
        count = 0
        try:
            while True:
                if count % 300 == 0:
                    print("Frame {}".format(count), end="\r")
                position = reader.get_position()
                frame = reader.decode_frame()
                packer.pack_hyper(position)
                packer.pack_float(frame.time)

                count += 1

        except EOFError:
            pass
        idx_fd.write(packer.get_buffer())
        packer.reset()

def main():
    if len(sys.argv) < 2:
        print("Usage {} xtc-file [xtc-file]".format(sys.argv[0]))
        return -1

    for filename in sys.argv[1:]:
        index_xtcfile(filename)



if __name__ == '__main__':
    sys.exit(main())