import argparse
from .logging import logger

import mdevaluate as md


def run(*args, **kwargs):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'xtcfile',
        help='The xtc file to index.',
    )
    parser.add_argument(
        '--tpr',
        help='The tprfile of the trajectory.',
        dest='tpr', default=None
    )
    parser.add_argument(
        '--nojump',
        help='Generate Nojump Matrices, requires a tpr file.',
        dest='nojump', action='store_true', default=False
    )
    parser.add_argument(
        '--debug',
        help='Set logging level to debug.',
        dest='debug', action='store_true', default=False
    )
    args = parser.parse_args()
    level = logger.DEBUG if args.debug else logging.INFO
    logger.basicConfig(level=level)

    md.open('', trajectory=args.xtcfile, topology=args.tpr, nojump=args.nojump)


if __name__ == '__main__':
    run()
