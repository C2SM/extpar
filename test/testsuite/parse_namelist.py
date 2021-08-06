import glob
import shutil
import os

def filter_string_list(string_list, filter_strings):
    '''
    Keep all entries of a list of strings, that do not
    contain an entry of filter_strings
    '''
    filtered_string = []
    for line in string_list:
        if not any(filter_string in line.lower() for filter_string in filter_strings):
            filtered_string.append(line)

    return filtered_string


def extract_data_files_from_namelists(dir_glob,list_glob,sep,comment):

    datafiles_raw = []

    datafiles_no_duplicates = []

    print('Extract data files contained in {}'.format(list_glob))

    test_folders = glob.glob(dir_glob)

    for test in test_folders:
        if 'intel' not in test:
            print('  - {}'.format(test))

            namelists_per_test = glob.glob(test + '/' + list_glob)

            for namelist in namelists_per_test:

                if 'INPUT_CHECK' not in namelist:

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

                                if ".nc" in line:
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

    datafiles_no_duplicates = list(dict.fromkeys(datafiles_raw))

    bad_words = ['buffer','icon','external','@','corine']

    return filter_string_list(datafiles_no_duplicates, bad_words)


def main():

    data_dir = '/store/c2sm/extpar_raw_data/linked_data'
    dest_dir = 'input-data'
    os.makedirs(dest_dir,exist_ok=True)

    # get data filenames
    files_fortran_nml = extract_data_files_from_namelists('data/*/*','INPUT_*','=','!')

    files_python_nml = extract_data_files_from_namelists('data/*/*','namelist.py',':','#')

    files_nml = list(dict.fromkeys(files_fortran_nml + files_python_nml))

    with open('files_to_copy.txt', 'w') as f:
        for file in files_nml:
            f.write(os.path.join(data_dir,file))
            f.write('\n')



if __name__ == "__main__":
    main()
