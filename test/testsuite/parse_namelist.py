
import glob
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

testlist = glob.glob('data/*/*')

datafiles_raw = []
datafiles_no_duplicates = []

print('Extract data files contained in namelists:')
for test in testlist:
    if 'intel' not in test:
        print('  - {}'.format(test))

        namelists_per_test = glob.glob(test + '/INPUT_*')

        for namelist in namelists_per_test:
            if 'INPUT_CHECK' not in namelist:

                with open(namelist, 'r') as f: 

                    # read line by line
                    for line in f:
                        line = line.rstrip().lstrip() 

                        # line is commented
                        if line.startswith("!"):
                            print('*** Ignore commented line: '
                                            f'{line}')

                        # valid entry in namelist
                        else:

                            if ".nc" in line:
                                split = line.split('=')

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

print(len(datafiles_raw))
datafiles_no_duplicates = list(dict.fromkeys(datafiles_raw))
print(len(datafiles_no_duplicates))

bad_words = ['buffer','icon','external','@','corine']
datafiles_clean = filter_string_list(datafiles_no_duplicates, bad_words)
print(len(datafiles_clean))

for file in datafiles_clean:
    print(file)

