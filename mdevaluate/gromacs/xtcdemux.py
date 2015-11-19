import sys
import os
from .logarithmic import is_log_step
from .reader import XTCReader

def next_log_step(step, per_decade):
    step += 1
    while not is_log_step(step, per_decade): step += 1
    return step

def next_lin_step(step, frequency):
    return step+frequency

class LogWriter:
    current_offset = 0
    log_step = 0
    fd = None

    def __init__(self, file, offset, log_freq):
        self.fd = open(file,'wb')
        self.offset = offset
        self.log_freq = log_freq

    def consume_frame(self, step, raw_frame):
        if step > self.log_step+self.offset:
            raise Exception("Missing frame {}".format(step))

        if step == self.log_step+self.offset:
            self.fd.write(raw_frame)
            self.log_step = next_log_step(self.log_step, self.log_freq)


def main():
    filename = sys.argv[1]
    reader = XTCReader(filename, load_cache_file=False)
    lin_freq = int(sys.argv[2])
    log_freq = int(sys.argv[3])
    log_restart = int(sys.argv[4])

    base,ext = os.path.splitext(filename)

    lin_file = base + '.lin' + ext
    lin_step = 0

    log_writers = []

    with open(lin_file,'wb') as lin_fd:
        try:
            while True:
                raw_frame, frame = reader.dump_frame()

                if frame.index % log_restart == 0:
                    print('Starting new log band: {}'.format(frame.index))
                    log_filename = '{}.{}.log{}'.format(base, len(log_writers), ext)
                    log_writers.append(LogWriter(log_filename, frame.index, log_freq))

                if frame.index > lin_step:
                    raise Exception("Missing frames")

                if frame.index == lin_step:
                    lin_fd.write(raw_frame)
                    lin_step = next_lin_step(lin_step, lin_freq)

                for writer in log_writers:
                    writer.consume_frame(frame.index, raw_frame)

        except EOFError:
            pass

if __name__ == '__main__':
    main()