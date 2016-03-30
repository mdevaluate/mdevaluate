import xdrlib
import io
import os
from .coordinates import decompress
from ..utils import hash_anything, merge_hashes
from functools import partial
import numpy as np

TRR_MAGIC = 1993
XTC_MAGIC = 1995
INDEX_MAGIC = 2015

def index_filename_for_xtc(xtcfile):
    """
    Get the filename for the index file of a xtc file.

    In general, the index file is located in the same directory as the xtc file.
    If the xtc file is located in a read-only directory (for the current user)
    and does not exist, the index file will be instead located in a subdirectory
    of ~/.xtcindex, of the current user.

    The directory for user index files can be changed by setting the environment
    variable 'XTCINDEXDIR'.

    Example:
        xtcfile = '/data/niels/test.xtc'

        # case 1: index file exists, or current user is niels
        index_filename_for_xtc(xtcfile) == '/data/niels/.test.xtc.xtcindex'

        # case 2: index file doesn't exist, and nor write permission in /data/niels
        index_filename_for_xtc(xtcfile) == '~/.xtcindex/data/niels/.test.xtc.xtcindex'

    Warning:
        At this point, the index file is not checked for validity. If an invalid
        index file is present in the xtc files directory, this will be used and an
        error will be risen by XTCReader.

    Warning:
        On most systems, the home directory is on the local drive so that the indexing
        has probably to do be done on every system if it can not be saved in the directory
        of the xtc file.
    """
    base = os.path.basename(xtcfile)
    filename = "." + base + ".xtcindex"

    dir = os.path.abspath(os.path.dirname(xtcfile))

    if not os.path.exists(os.path.join(dir, filename)) and  not os.access(dir, os.W_OK):
        if 'XTCINDEXDIR' in os.environ:
            index_dir = os.environ['XTCINDEXDIR']
        else:
            index_dir = os.path.join(os.environ['HOME'], '.xtcindex')
        dir = os.path.join(index_dir, dir.lstrip('/'))
        os.makedirs(dir, exist_ok=True)

    return os.path.join(dir, filename)


class NumpyUnpacker(xdrlib.Unpacker):

    def unpack_float_array(self, n_items):
        i = self.get_position()
        j = i + 4 * n_items
        self.set_position(j)
        data = self.get_buffer()[i:j]

        if len(data) < 4:
            raise EOFError

        ret = np.frombuffer(data, '>f')
        if len(ret) != n_items:
            raise EOFError
        return ret

    def unpack_double_array(self, n_items):
        i = self.get_position()
        j = i + 8 * n_items
        self.set_position(j)
        data = self.get_buffer()[i:j]

        if len(data) < 8:
            raise EOFError

        ret = np.frombuffer(data, '>d')
        if len(ret) != n_items:
            raise EOFError
        return ret


class InvalidMagicException(Exception):
    pass


class InvalidIndexException(Exception):
    pass


class UnknownLenError(Exception):
    pass


class SubscriptableReader:

    def __init__(self, fd):
        self.fd = fd
        self.len = os.stat(self.fd.fileno()).st_size

    def __getitem__(self, r):
        if isinstance(r, slice):
            if r.start is not None:
                if r.start < 0:
                    self.fd.seek(r.start, io.SEEK_END)
                else:
                    self.fd.seek(r.start, io.SEEK_SET)

            if r.step is not None:
                raise NotImplementedError

            return self.fd.read(r.stop - r.start)
        else:
            self.fd.seek(r, io.SEEK_SET)
            return self.fd.read(1)

    def __len__(self):
        return self.len


class XTCFrame:
    __slots__ = ['_coordinates', 'index', 'time', 'box', '_compressed_coordinates']

    def __init__(self):
        self._coordinates = None

    @property
    def coordinates(self):
        if self._coordinates is None:
            self._coordinates = self._compressed_coordinates()
        return self._coordinates


class BaseReader:

    def __init__(self, file):
        self.fd = open(file, 'rb')
        self.filename = file

        self.reader = SubscriptableReader(self.fd)
        self.xdr_reader = NumpyUnpacker(self.reader)

    def __exit__(self, *_):
        self.fd.close()

    def __enter__(self, *_):
        return self

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['fd']
        del state['reader']
        del state['xdr_reader']
        state['filepos'] = self.get_position()
        return state

    def __setstate__(self, state):
        filepos = state.pop('filepos')
        self.__dict__.update(state)
        self.fd = open(self.filename, 'rb')
        self.reader = SubscriptableReader(self.fd)
        self.xdr_reader = NumpyUnpacker(self.reader)
        self.set_position(filepos)

    def get_position(self):
        return self.fd.tell()

    def set_position(self, pos):
        return self.fd.seek(pos)

    def skip_frames(self, num):
        for i in range(num):
            self.skip_frame()

    def skip_frame(self):
        self._read_header()
        self._read_frame()

    def dump_frame(self):
        raise NotImplemented

    def _read_header(self):
        raise NotImplemented

    def _read_frame(self):
        raise NotImplemented

XTC_HEADER_SIZE = 4 * 4
XTC_FRAME_HEADER_SIZE = 9 * 4 + 10 * 4


class XTCReader(BaseReader):

    len_available = False

    def __init__(self, xtcfile, load_cache_file=True):
        super().__init__(xtcfile)
        self._cache = [0]
        self._times = None
        self.current_frame = 0

        index_file_name = index_filename_for_xtc(xtcfile)
        if load_cache_file:
            self.load_index(index_file_name)

    def load_index(self, filename):
        xtc_stat = os.stat(self.filename)
        c_time = int(xtc_stat.st_ctime)
        m_time = int(xtc_stat.st_mtime)
        size = xtc_stat.st_size

        with open(filename, 'rb') as index_file_fd:
            # TODO: Is NumpyUnpacker even necessary at this point?
            #       Seems like xdrlib.Unpacker would be sufficient here ...
            reader = NumpyUnpacker(SubscriptableReader(index_file_fd))

            if reader.unpack_hyper() != INDEX_MAGIC:
                raise InvalidMagicException
            if reader.unpack_hyper() != c_time:
                raise InvalidIndexException
            if reader.unpack_hyper() != m_time:
                raise InvalidIndexException
            if reader.unpack_hyper() != size:
                raise InvalidIndexException

            self._cache = []
            self._times = []
            try:
                while True:
                    self._cache.append(reader.unpack_hyper())
                    self._times.append(reader.unpack_float())
            except EOFError:
                pass
        self.len_available = True

    def _raw_header(self):
        if len(self._cache) == self.current_frame:
            self._cache.append(self.get_position())
        return self.fd.read(XTC_HEADER_SIZE)

    def _raw_frame(self):
        frame_header = self.fd.read(XTC_FRAME_HEADER_SIZE)
        blob_size = xdrlib.Unpacker(frame_header[-4:]).unpack_int()
        blob_size = (blob_size + 3) // 4 * 4   # Padding to 4 bytes
        frame_blob = self.fd.read(blob_size)

        self.current_frame += 1
        return frame_header + frame_blob

    def _unpack_header(self, raw_header):
        unpacker = xdrlib.Unpacker(raw_header)

        magic = unpacker.unpack_int()

        if magic != XTC_MAGIC:
            raise InvalidMagicException

        n_atoms = unpacker.unpack_int()
        step = unpacker.unpack_int()
        time = unpacker.unpack_float()

        return n_atoms, step, time

    def _read_header(self):
        raw_header = self._raw_header()
        return self._unpack_header(raw_header)

    def _unpack_frame(self, raw_frame):
        unpacker = xdrlib.Unpacker(raw_frame)

        raw_box = unpacker.unpack_farray(9, unpacker.unpack_float)
        box = np.array(raw_box).reshape(3, 3)
        num_coords = unpacker.unpack_int()
        precision = unpacker.unpack_float()
        maxint = unpacker.unpack_farray(3, unpacker.unpack_int)
        minint = unpacker.unpack_farray(3, unpacker.unpack_int)
        smallindex = unpacker.unpack_int()

        blob_len = unpacker.unpack_int()
        blob = unpacker.unpack_fopaque(blob_len)

        return box, precision, num_coords, minint, maxint, smallindex, blob

    def _read_frame(self):
        raw_frame = self._raw_frame()
        return self._unpack_frame(raw_frame)

    def decode_frame(self, header=None, body=None):
        n_atoms, step, time = header or self._read_header()
        box, precision, num_coords, minint, maxint, smallindex, blob = body or self._read_frame()
        coordinates = partial(decompress, num_coords, precision, minint, maxint, smallindex, blob)

        frame = XTCFrame()
        frame._compressed_coordinates = coordinates
        frame.index = step
        frame.time = time
        frame.box = box
        return frame

    def dump_frame(self):
        """
        :return: Tuple: The binary data of the frame, the frame itself
        """

        raw_header = self._raw_header()
        header = self._unpack_header(raw_header)

        raw_body = self._raw_frame()
        body = self._unpack_frame(raw_body)

        return (raw_header + raw_body, self.decode_frame(header, body))

    def cached_position(self, item):
        try:
            return self._cache[item]
        except IndexError:
            return None

    def __getitem__(self, item):
        position = self.cached_position(item)

        if position is not None:
            self.set_position(position)
            self.current_frame = item
            return self.decode_frame()
        # TODO: Use elif and one single return at the end
        # Also: skip_frames is in fact not implemented at all!
        # TODO: Catch EOFError and raise IndexError
        if self.current_frame <= item:
            self.skip_frames(item - self.current_frame)
            return self.decode_frame()
        else:
            self.set_position(0)
            self.current_frame = 0
            self.skip_frames(item)
            return self.decode_frame()

    def __len__(self):
        if self.len_available:
            return len(self._cache)
        raise UnknownLenError

    def times_of_indices(self, indices):
        return [self[i].time for i in indices]

    def __hash__(self):
        return merge_hashes(hash_anything(self.filename), hash_anything(self._cache))


class TRRHeader:
    __slots__ = ['ir_size', 'e_size', 'box_size', 'vir_size', 'pres_size', 'top_size', 'sym_size', 'x_size', 'v_size',
                 'f_size', 'n_atoms', 'step', 'nre', 't', '_lambda', 'is_double']

    @property
    def frame_size(self):
        return self.ir_size + \
            self.e_size + \
            self.box_size + \
            self.vir_size + \
            self.top_size + \
            self.sym_size + \
            self.x_size + \
            self.v_size + \
            self.f_size


class TRRFrame:
    __slots__ = ['header', 'box', 'x', 'v', 'f']


class TRRReader(BaseReader):

    def __iter__(self):
        return TRRIterator(self.filename)

    def _unpack_float(self, is_double):
        if is_double:
            return self.xdr_reader.unpack_double()
        else:
            return self.xdr_reader.unpack_float()

    def _unpack_np_float(self, is_double, *dim):
        total = np.product(dim)

        if is_double:
            box = self.xdr_reader.unpack_double_array(total)
        else:
            box = self.xdr_reader.unpack_float_array(total)

        box = box.reshape(*dim)
        return box

    def _read_header(self):

        data = self.fd.read(8)  # 4 Magic + 4 VersionStringLen
        unpacker = NumpyUnpacker(data)

        magic = unpacker.unpack_int()

        if magic != TRR_MAGIC:
            raise InvalidMagicException

        version_string_len = unpacker.unpack_int()
        data = self.fd.read(version_string_len + 55)  # 4 Magic + 4 VersionStringLen
        unpacker.reset(data)

        version = unpacker.unpack_string()
        assert version == b'GMX_trn_file'

        header = TRRHeader()

        header.ir_size = unpacker.unpack_int()
        header.e_size = unpacker.unpack_int()
        header.box_size = unpacker.unpack_int()
        header.vir_size = unpacker.unpack_int()
        header.pres_size = unpacker.unpack_int()
        header.top_size = unpacker.unpack_int()
        header.sym_size = unpacker.unpack_int()
        header.x_size = unpacker.unpack_int()
        header.v_size = unpacker.unpack_int()
        header.f_size = unpacker.unpack_int()
        header.n_atoms = unpacker.unpack_int()
        header.step = unpacker.unpack_int()
        header.nre = unpacker.unpack_int()

        if header.x_size is not None:
            header.is_double = (header.x_size // header.n_atoms) == 8
        elif header.v_size is not None:
            header.is_double = (header.v_size // header.n_atoms) == 8
        elif header.f_size is not None:
            header.is_double = (header.f_size // header.n_atoms) == 8

        if header.is_double:
            data = self.fd.read(16)
        else:
            data = self.fd.read(8)
        unpacker.reset(data)
        self.xdr_reader = unpacker

        header.t = self._unpack_float(header.is_double)
        header._lambda = self._unpack_float(header.is_double)

        return header

    def _read_frame(self, header):
        frame = TRRFrame()
        frame.header = header

        data = self.fd.read(header.frame_size)
        self.xdr_reader = NumpyUnpacker(data)

        frame.box = None
        if header.box_size:
            frame.box = self._unpack_np_float(header.is_double, 3, 3)

        if header.vir_size:
            for i in range(9):
                self._unpack_float(header.is_double)

        if header.pres_size:
            for i in range(9):
                self._unpack_float(header.is_double)

        frame.x = None
        if header.x_size:
            frame.x = self._unpack_np_float(header.is_double, header.n_atoms, 3)

        frame.v = None
        if header.v_size:
            frame.v = self._unpack_np_float(header.is_double, header.n_atoms, 3)

        frame.f = None
        if header.f_size:
            frame.f = self._unpack_np_float(header.is_double, header.n_atoms, 3)

        return frame

    def decode_frame(self):
        h = self._read_header()
        return self._read_frame(h)

    def get_position(self):
        return self.fd.tell()

    def set_position(self, pos):
        self.fd.seek(pos)


class TRRIterator():

    def __init__(self, file):
        self.__reader = TRRReader(file)
        self.__reader.__enter__()

    def __next__(self):
        try:
            return self.__reader.decode_frame()
        except EOFError as err:
            self.__reader.__exit__()
            raise StopIteration

    def __iter__(self):
        return self
