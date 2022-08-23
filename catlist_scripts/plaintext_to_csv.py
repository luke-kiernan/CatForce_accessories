'''Catlist text to csv.'''

import sys
import csv

assert(len(sys.argv) == 3), "Requires 2 arguments: [catalyst plaintext file] [path to save catalyst csv file]"

def CheckFormatting(line, subarray, entryType):
    assert('o' in subarray[0]), f"RLE not found for entry {subarray[0]} in line: \n\t \'{line}\'"
    assert(subarray[1].lstrip('-').isnumeric() and subarray[2].lstrip('-').isnumeric()), f"x,y position not found for entry {entryType} in line: \n\t \'{line}\'"

maxForbidden = 0
csvEntries = []
periodic = False
with open(sys.argv[1]) as f:
    
    for line in f:
        data = line.rstrip().split(' ')
        if data[0] == 'cat':
            assert(data[2].isnumeric()), f"Stable interval not found in line: \n\t \'{line}\'"
            csvEntries.append(["[name]", data[2]])
            # cat rle stable x y sym 
            # 0   1    2     3 4  5
            posData = [data[1], data[3], data[4]]
            CheckFormatting(line, posData, "catalyst")
            csvEntries[-1] += posData
            assert(data[5] in ['@', '*', '.','/', 'x', '+', '|', '-']),f"Symmetry character wrong in line: \n\t \'{line}\'"
            csvEntries[-1].append(data[5])
            if 'period' in data:
                csvEntries[-1].append(data[data.index('period')+1])
                periodic = True
            for keyword in ["required"]: #, "locus"]:
                if keyword in data:
                    i = data.index(keyword)
                    CheckFormatting(line, data[i+1:i+4], keyword)
                    csvEntries[-1] += data[i+1:i+4]
                else:
                    csvEntries[-1] += ['','','']  
            numForbidden = 0
            for i in range(len(data)):
                if data[i] == "forbidden":
                    numForbidden += 1
                    CheckFormatting(line, data[i+1:i+4], "forbidden")
                    csvEntries[-1] += data[i+1:i+4]
                    i += 4
            maxForbidden = max(numForbidden, maxForbidden)
for row in csvEntries:
    while(len(row) < 6*6*maxForbidden):
        row.append('')

headerList = ['name', 'absence', 'rle', 'dx', 'dy', 'symType']
if periodic:
    headerList.append('period')
headerList +=['required', 'req dx', 'req dy']#, 'locus rle', 'locus x', 'locus y']
for i in range(1,maxForbidden+1):
    headerList += [f'forbidden {i}', f'forbid {i} dx', f'forbid {i} dy']
with open(sys.argv[2], 'w') as file:
    dw = csv.writer(file, delimiter=',')
    dw.writerow(headerList)
    for row in csvEntries:
        dw.writerow(row)
