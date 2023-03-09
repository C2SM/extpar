#!/usr/bin/env python3
import glob
import shutil
import os
import argparse


def filter_string_list(string_list, filter_strings):
    '''
    Keep all entries of a list of strings, that do not
    contain an entry of filter_strings
    '''
    filtered_string = []
    for line in string_list:
        if not any(filter_string in line.lower()
                   for filter_string in filter_strings):
            filtered_string.append(line)

    return filtered_string


def parse_for_data_files(namelist, comment, sep):
    '''
    Parse a namelist line-by-line and find input
    datafiles. For data in tiles create a single
    entry for each file
    '''

    datafiles_raw = []

    with open(namelist, 'r') as f:

        # read line by line
        for line in f:
            line = line.rstrip().lstrip()

            # line is commented
            if line.startswith(comment):
                print('*** Ignore commented line: '
                      '{}'.format(line))

            # valid entry in namelist
            else:

                if ".nc" in line or ".data" in line:
                    split = line.split(sep)

                    # return last element of split
                    raw_variable = split[-1].strip()

                    characters_to_strip = ["'", ",", '"']
                    for character in characters_to_strip:
                        raw_variable = raw_variable.strip(character)

                    # for entries with more than one data file i.e. Aster tiles
                    if " " in raw_variable:
                        tile_data = raw_variable.split(" ")
                        characters_to_strip = ["'", ",", '"']
                        for tile in tile_data:
                            for character in characters_to_strip:
                                tile = tile.strip(character)
                            if len(tile) != 0:
                                datafiles_raw.append(tile.strip("'"))
                    else:
                        datafiles_raw.append(str(raw_variable.rstrip("'")))

    return datafiles_raw


def extract_data_files_from_namelists(dir_glob, list_glob, sep, comment):

    files_from_namelists = []

    print('Extract data files contained in {}'.format(list_glob))

    test_folders = glob.glob(dir_glob)

    for test in test_folders:
        if 'intel' not in test:
            print('  - {}'.format(test))

            namelists_per_test = glob.glob(test + '/' + list_glob)

            for namelist in namelists_per_test:

                files_from_namelists += parse_for_data_files(
                    namelist, comment, sep)

    files_no_duplicates = list(dict.fromkeys(files_from_namelists))

    # lowercase words
    bad_words = ['buffer', 'icon', 'external', '@', '_tg','clim_tsea_corine','clim_t2m_corine']

    return filter_string_list(files_no_duplicates, bad_words)


if __name__ == "__main__":

    # parsing arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--dir_glob',
                        dest='dir_glob',
                        default='data/*/*',
                        type=str,
                        help='glob-pattern for directories')

    args = parser.parse_args()

    # get data filenames
    files_fortran_nml = extract_data_files_from_namelists(
        args.dir_glob, 'INPUT_*', '=', '!')

    files_python_nml = extract_data_files_from_namelists(
        args.dir_glob, 'namelist.py', ':', '#')

    files_nml = list(dict.fromkeys(files_fortran_nml + files_python_nml))

    # write to file
    with open('transfer.txt', 'w') as f:
        for file in files_nml:
            f.write(file)
            f.write('\n')

    print('')
    print('Files written to transfer.txt')
