'''Convert CSV catalyst list to plaintext catalyst list.'''

import sys
import csv
import re

namesToIgnore=[] # "cat *", "unnamed", "two blocks *", "two fishhooks *"]

assert(len(sys.argv) == 2 or len(sys.argv) == 3), "Too many or to few arguments.\nPlease specify the name of the CSV-style catalyst list."
includeNames = (len(sys.argv) == 3 and sys.argv[2][0].lower() == 'n')
lastPeriod = -1
with open(sys.argv[1], 'r', newline='', encoding = 'utf-8') as f:
    try:
        lineNum = 0
        reader = csv.DictReader(f)
        for data in reader:
            keys = data.keys()
            if includeNames:
                if 'period' in data and data['period'] != lastPeriod:
                    print(f"\n# ============ period {data['period']} ============")
                if all([re.search(regexPat,data['name']) is None for regexPat in namesToIgnore]):
                    print(f"# ==={data['name']}===")
            if 'stable' in keys:
                stableInterval = data['stable']
            else:
                stableInterval = data['absence']
            symChar = data['symType'] if 'symType' in keys else data['symmetry']
            print(f"cat {data['rle']} {stableInterval} {data['dx']} {data['dy']} {data['symType']}", end = '')
            if 'required' in keys or 'req rle' in keys:
                reqRLE = data['required'] if 'required' in keys else data['req rle']
                if reqRLE != '':
                    print(f" required {reqRLE} {data['req dx']} {data['req dy']}", end = '')
            if 'antirequired' in keys or 'antireq rle' in keys:
                antireqRLE = data['antirequired'] if 'antirequired' in keys else data['antireq rle']
                if antireqRLE != '':
                    print(f" antirequired {antireqRLE} {data['antireq dx']} {data['antireq dy']}", end = '')
            if 'locus' in keys and 'locus' in data and 'o' in data['locus']:
                print(f" locus {data['locus']} {data['locus dx']} {data['locus dy']}", end = '')
            n = 1
            while f'forbidden {n}' in keys and data[f'forbidden {n}'] != '':
                print(f" forbidden {data[f'forbidden {n}']} {data[f'forbid {n} dx']} {data[f'forbid {n} dy']}", end = '')
                n += 1
            if 'period' in keys:
                print(f" period {data['period']}", end = '')
                lastPeriod = data['period']
            else:
                lastPeriod = 1
            print('\n', end = '')
            lineNum += 1
    except Exception as e:
        print(repr(e))
        print(f"check line number {lineNum}")