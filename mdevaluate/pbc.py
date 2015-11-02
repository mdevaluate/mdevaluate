def pbc_diff(v1, v2, box):
    v = v1-v2
    if box is None:
        return v

    v -= (v > box/2)*box
    v += (v < -box/2) * box

    return v
