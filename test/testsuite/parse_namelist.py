
import glob
namelists = glob.glob('INPUT_*')
#namelists = ['INPUT_ORO','INPUT_AOT']
variable = 'test'
type_to_convert = str
clean_files = []
for namelist in namelists:
    with open(namelist, 'r') as f: 
        # read line by line
        for line in f:
            line = line.rstrip().lstrip() 

            # line is commented
            if line.startswith("!"):
                logging.warning('Ignore commented line: '
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


                    if " " in raw_variable:
                        split = raw_variable.split(" ")
                        characters_to_strip = ["'", ",", '"']
                        for file in split:
                            for character in characters_to_strip:
                                file = file.strip(character)
                            if len(file) != 0:
                                    clean_files.append(file)
                    else:
                        if not "BUFFER" in raw_variable.upper():
                            if not "ICON" in raw_variable.upper():
                                clean_files.append(str(raw_variable))

for file in clean_files:
    print(file)

