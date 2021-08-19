import re
filename = 'INPUT_AOT'
pattern = re.compile(r'nc')
with open(filename) as f:
    for line in f:
            match=pattern.search(line)
            if match:
                print(match.group(1))
            pattern.findall(line)
