"""Compatibility entry point; implementation is installed with KinBot."""
from kinbot import requested_orbital_recon as _implementation


def __getattr__(name):
    return getattr(_implementation, name)


if __name__ == '__main__':
    _implementation.main()
