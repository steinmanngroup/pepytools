import sys
import numpy
import numpy.linalg

from .constants import AATOBOHR


def read_coordinates(file):
    coordinates = []
    labels = []
    lines = 0
    read = True
    factor = 1.0
    while read:
        line = file.readline()
        tokens = line.split()
        if lines == 0:
            lines = int(tokens[0])  # 'the title line'
            line = file.readline()
            if 'AA' in line.upper():
                factor = AATOBOHR
            continue

        lines = lines - 1
        coordinates.append(list(map(float, tokens[1:4])))
        labels.append(tokens[0])
        read = (lines > 0)
    return factor * numpy.array(coordinates), labels


def read_polarizabilites(file):
    tensors = []
    lines = 0
    read = True
    while read:
        line = file.readline()
        tokens = line.split()

        if lines == 0 and read:
            line = file.readline()
            tokens = line.split()
            lines = int(tokens[0])
            continue

        tensors.append(list(map(float, tokens[1:])))
        lines = lines - 1
        read = (lines > 0)
    return tensors


def read_multipoles(file, order=0):
    multipoles = []
    lines = 0
    read = True
    while read:
        line = file.readline()
        tokens = line.split()

        if lines == 0 and read:
            lines = int(tokens[0])
            continue

        data = list(map(float, tokens[1:]))
        multipoles.append(numpy.array(data))
        lines = lines - 1
        read = (lines > 0)

    return multipoles


def read_exclusionlists(file):
    exclusions = {}
    lines = 0
    read = True
    excl_counter = 0
    while read:
        line = file.readline()
        tokens = line.split()

        if lines == 0 and read:
            lines = int(tokens[0])
            continue

        data = list(map(int, tokens))
        exclusions[excl_counter] = numpy.array(data[1:]) - 1  # the first one is just the site, the rest are what they are not supposed to mix with. I guess.
        lines = lines - 1
        read = (lines > 0)
        excl_counter += 1

    return exclusions


def read_potential_from_file(filename):
    readers = {'COORD': read_coordinates,
               'MULT': read_multipoles,
               'POLAR': read_polarizabilites,
               'EXCL': read_exclusionlists}

    f = open(filename, 'r')

    line = "\n"
    reader = None
    coordinates = None
    atom_labels = None
    polarizabilites = None
    multipoles = {}
    exclusion_list = None
    while len(line) > 0:
        line = f.readline()
        #print "line:", line

        for key in readers.keys():
            if key in line:
                reader = readers[key]
                break

        if reader is not None and reader == read_coordinates:
            coordinates, atom_labels = reader(f)
            reader = None

        if reader is not None and reader == read_polarizabilites:
            polarizabilites = reader(f)
            reader = None

        if reader is not None and reader == read_multipoles:
            order_line = line.split()
            try:
                order = int(order_line[1])
            except IndexError:
                line = f.readline()
                order_line = line.split()

                # if the file has no polarizabilities it just ends here
                if len(order_line) == 0:
                    break
                order = int(order_line[1])
            multipoles[order] = reader(f, order)

        if reader is not None and reader == read_exclusionlists:
            exclusion_list = reader(f)
            reader = None

    f.close()

    return (coordinates, atom_labels, multipoles, polarizabilites, exclusion_list)
