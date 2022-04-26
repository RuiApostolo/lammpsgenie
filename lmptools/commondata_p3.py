########
# module containing commonly used data reading and analysis tools
# Author: Rui Apóstolo
# email: ruiapostolo@gmail.com
# Inspired by previous work by Michael Doig
from .atoms import getNatoms


def getDumpTSRange(lines):
    """
    Reads dump file to get list of timestep numbers.

    Parameters
    ----------
    lines : list of str
        Read list of lines.

    Returns
    -------
    tsrange : list of int
        List of timestep numbers.

    Notes
    -----
    Fails hard if lines don't contain the target string.
    """

    tsrange = []
    for i, _ in enumerate(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            tsrange.append(int(lines[i+1]))
    return tsrange


def readTS(lines, tsnum, atomnames, header=9):
    """
    Reads dumpfile and returns the atomic information of a timestep.

    Parameters
    ----------
    lines : list of str
        Reads list of lines from dumpfile.

    tsnum : int
        Number of the frame to read.

    atomnames : dict
        Dictionary of atom types and names. Use getAtomType().

    header : int, default=9
        Number of lines in the dumpfile 'header', from 'ITEM: TIMESTEP' to
        'ITEM: ATOM', inclusive.
    Returns
    -------
    traj : dict
        A trajectory dict with more than one trajectory. Keys are the
        timestep labels.
    TODO: remove `ts`?
    """

    traj = {}
    natoms = getNatoms(lines)
    nlines = natoms + header
    firstline = ((tsnum-1)*nlines)
    lastline = firstline + nlines-1
    for line in range(firstline, lastline):
        if lines[line].startswith("ITEM: TIMESTEP"):
            ts = int(lines[line+1])
            traj[ts] = {}
            traj[ts]["atom"] = {}
        if any(s in lines[line] for s in ["BOX BOUNDS", "xy xz yz pp pp"]):
            boxx = list(map(float, lines[line + 1].split()))
            boxy = list(map(float, lines[line + 2].split()))
            boxz = list(map(float, lines[line + 3].split()))
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            traj[ts]["boxsize"] = (lx, ly, lz)
        if lines[line].startswith("ITEM: ATOMS"):
            labels = lines[line].split()
            for atomid in range(1, natoms + 1):
                atomvals = lines[line+atomid].split()
                traj[ts]["atom"][atomid] = {}
                for label_idx, k in enumerate(labels[3:], start=1):
                    if k == 'type':
                        traj[ts]["atom"][atomid][k] = int(atomvals[label_idx])
                    else:
                        traj[ts]["atom"][atomid][k] = float(
                            atomvals[label_idx])
                # Change from numerical atom type (from dump) to atom type
                # name (from data file)
                traj[ts]["atom"][atomid]["type"] = \
                    atomnames[traj[ts]["atom"][atomid]["type"]]
    return traj


def getTrajTSRange(traj, ts1, ts2):
    """
    Select a partial timestep range from a trajectory.

    Parameters
    ----------
    traj : dict
        A trajectory dict with more than one trajectory. Keys are the
        timestep labels.
    ts1 : int or 'first'
        The timestep label for the start of the range.
        If 'first' uses the first timestep of the trajectory, regardless
        of number.
    ts2 : int or 'last'
        The timestep label for the end of the range.
        If 'last' uses the last timestep of the trajectory, regardless
        of number.

    Returns
    -------
    tsrange : list of ints
        List of timestep labels.
    """

    sortedts = sorted(traj.keys())
    try:
        ts1 = 0 if ts1 == 'first' else sortedts.index(ts1)
        ts2 = None if ts2 == 'last' else sortedts.index(ts2)
        if ts2 is None:
            tsrange = sortedts[ts1:None]
        else:
            tsrange = sortedts[ts1:ts2 + 1]
        return tsrange
    except(ValueError):
        return []
