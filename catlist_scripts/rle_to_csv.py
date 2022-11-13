from math import ceil
from copy import copy
import golly as g
import os
import csv
import sys
from itertools import product
MAX_PERIOD = 25

def TrimCellList(cellList):
    if len(cellList) % 6 == 1: cellList.pop()
    if len(cellList) % 3 != 0:
        g.warn(str(cellList))
        g.warn(str(len(cellList)))
    assert(len(cellList) % 3 == 0)

def AddPadding(cellList):
    if len(cellList) % 2 == 0:
        cellList.append(0)

def LifeHistoryToLife(cellList, states = None):
    if states is None:
        states = [1,3,5]
    if len(states) == 0:
        return []
    TrimCellList(cellList)
    newList = []
    for i in range(0,len(cellList),3):
        if cellList[i+2] in states:
            newList += cellList[i:i+2]
    return newList

def GetX0Y0(rectList, state = None):
    cellCoords = ExtractState(rectList, state)
    return (min([t[0] for t in cellCoords]), min([t[1] for t in cellCoords]))

def LifeHistoryToRLE(cellList, states=None):
    if states is None:
        states = [1,3,5]
    TrimCellList(cellList)
    relevantCells = []
    for i in range(0, len(cellList), 3):
        if cellList[i+2] in states:
            relevantCells += cellList[i:i+2]
    g.store(relevantCells, 'temp.rle')
    # it really seems like there ought to be a way to do this
    # without creating a temporary file, but this works, so whatever.
    with open('temp.rle', 'r') as f:
        rle = ''
        for line in f:
            if line.startswith('x ='):
                data = line.split(', ')
                w = int(data[0].removeprefix('x = '))
                h = int(data[1].removeprefix('y = '))
            if not (line.startswith('#') or line.startswith('x =')):
                rle += line.rstrip()
    os.remove('temp.rle')
    for char in ['A', 'B', 'C', 'D', 'E']:
        rle = rle.replace(char, 'o')
    rle = rle.replace('.', 'b')
    return rle #, w, h

def ExtractState(rectList, states = None, sortList = False):
    '''Return a list of the coordinates of all cells with a given state [or, 1/3/5 if None].'''
    if states is None:
        states = [1,3,5]
    cellList = g.getcells(rectList)
    TrimCellList(cellList)
    coords = []
    for i in range(0, len(cellList), 3):
        if cellList[i+2] in states:
            coords.append(tuple(cellList[i:i+2]))
    if sortList: coords.sort(key = lambda t: tuple(reversed(t)))
    return coords


def GetMinRect(rectList):
    g.select(rectList)
    g.strink()
    return g.getselrect()

def SplitIntoBlocksByRepeatedZeros(stringOfZerosAndOnes, howManyZeros):
    blocks = []
    blockStart = stringOfZerosAndOnes.find('1')
    blockEnd = stringOfZerosAndOnes.find(howManyZeros*'0', blockStart)
    while(blockEnd != -1):
        blocks.append((blockStart, blockEnd))
        blockStart = stringOfZerosAndOnes.find('1', blockEnd)
        blockEnd = stringOfZerosAndOnes.find(howManyZeros*'0', blockStart)
    if blockStart != -1:
        blocks.append((blockStart, len(stringOfZerosAndOnes)-1))
    return blocks

inGolly = True
try: # compatibility with run_in_golly.sh
    g.open(pseudo_argv[1])
    inGolly = False
    fname = pseudo_argv[2]
    if len(pseudo_argv) > 3:
        periodic = pseudo_argv[3] == 'periodic'
    if len(pseudo_argv) > 4:
        minRowGap = int(pseudo_argv[4])
        minEntryGap = int(pseudo_argv[4])
    else:
        minRowGap = 4
        minEntryGap = 4

except NameError:

    g.note('Program for producing lists for CatForce.')
    g.note('1 catalyst per row. \n First entry: sparker/catalyst, with state 4/5 for required cells. \n Optionally: locus as second entry, with state 4/5 for locus. \n Rest are forbidden, only state 3 cells.')
    #g.note(f'Include at least {minRowGap} cells between rows, and {minEntryGap}  cells between entries in the in the row.')
    g.note('Patterns don\'t need to exactly line up in each row, but they shouldn\'t be off by more than 7.')
    g.note('In the first entry, draw a single white cell at the spark--the catalyst will be centered there.')
    g.note('Symmetry and max absence time must be entered by hand.')

    periodic = g.getstring("periodic catalyst list? y/n \n(Needed for required vs antirequired)") == 'y'
    minGap = int(g.getstring("Please enter minimum gap between entries and rows:"))

    minRowGap = int(minGap)
    minEntryGap = int(minGap)
    fname = g.savedialog("Save csv catalyst list", "*.csv")
    if fname == '':
        g.exit("File save canceled.")
# list of list of rectangle lists
# each list is the entries in 1 row.
rectanglesInRow = []
[xMin, yMin , w, h ] = g.getrect()
# figure out spacing, break up rle into rows and entries in each row
rowsNonempty = ''
for y in range(yMin, yMin+h+1):
    rowsNonempty += '1' if len(g.getcells([xMin,y,w,1])) > 1 else '0' 
for rowBlock in SplitIntoBlocksByRepeatedZeros(rowsNonempty, minRowGap):
    rectanglesInRow.append([])
    rowHeight = rowBlock[1]-rowBlock[0]+1
    rowStartY = yMin+rowBlock[0]
    colsNonempty = ''
    for x in range(xMin, xMin+w+1):
        colsNonempty += '1' if len(g.getcells([x,rowStartY,1,rowHeight])) > 1 else '0' 
    for colBlock in SplitIntoBlocksByRepeatedZeros(colsNonempty, minEntryGap):
        rectanglesInRow[-1].append([xMin+colBlock[0], rowStartY, colBlock[1]-colBlock[0]+1, rowHeight])

# done figuring out spacing
maxForbidden = 0
allCSVRows = []

catNumber = 1


for entryBboxes in rectanglesInRow:
    dataDict = {}
    
    stateRLE = LifeHistoryToRLE(g.getcells(entryBboxes[0]))
    
    patCells = LifeHistoryToLife(g.getcells(entryBboxes[0]))
    if len(patCells) == 0:
        g.select(entryBboxes[0])
        g.fitsel()
        g.warn("problem with selected box")
    x0Pat, y0Pat = min(patCells[::2]), min(patCells[1::2])
    w, h = max(patCells[::2])-x0Pat+1,max(patCells[1::2])-y0Pat+1


    # center = (ceil(x0Pat+w/2)//1, ceil(y0Pat+h/2)//1)
    # catForce universe [post-(dx,dy)] catalyst
    #          = golly universe catalyst, shifted by -center.

    # improvement: white = center.
    
    if len(ExtractState(entryBboxes[0], [3])) == 1:
        center = ExtractState(entryBboxes[0], [3])[0]
    else:
        center = (ceil(x0Pat+w/2)//1, ceil(y0Pat+h/2)//1)


    dataDict['rle'] = stateRLE
    dataDict['dx'] = -(center[0]-x0Pat)
    dataDict['dy'] = -(center[1]-y0Pat)
    if periodic:
        reqStates = [4,5]
        antiReqStates = []
    else:
        reqStates = [5]
        antiReqStates = [4]
    reqCells = LifeHistoryToLife(g.getcells(entryBboxes[0]), reqStates)
    antiReqCells = LifeHistoryToLife(g.getcells(entryBboxes[0]), antiReqStates)
    
    cellStates = [reqStates,antiReqStates]
    cells = [reqCells, antiReqCells]
    names = ['required', 'antirequired']
    shorthands = ['req', 'antireq']
    for i in range(2):
        if len(cells[i]) >= 1:
            dataDict[names[i]] = LifeHistoryToRLE(g.getcells(entryBboxes[0]), cellStates[i])
            x0, y0= min(cells[i][::2]), min(cells[i][1::2])
            dataDict[shorthands[i]+' dx'] = x0 - center[0]
            dataDict[shorthands[i]+' dy'] = y0 - center[1]
        if periodic:
            break # no antirequired


    # want (dxReq, dyReq) such that
    # default position required shifted by (dxReq, dyReq)
    #  equals golly universe required, shifted by -center.
    # that is, (x0Req, y0Req) - center = (dxReq, dyReq)

    patCellsCentered = g.transform(g.getcells(entryBboxes[0]), -center[0], -center[1])
    patLifeCellsCentered = LifeHistoryToLife(patCellsCentered)
    patCoordsCentered = [tuple(patLifeCellsCentered[i:i+2]) 
                        for i in range(0, len(patLifeCellsCentered), 2)]

    # make sure it really is periodic or stable.
    checkUntil = 1 if not periodic else MAX_PERIOD
    g.putcells(g.getcells(entryBboxes[0]), (xMin-100)-x0Pat, (yMin-100)-y0Pat)
    sortedCellsAtStart = ExtractState([xMin-100, yMin-100, 64, 64], [1,3,5], sortList = True)
    g.select([xMin-164, yMin-164, 128, 128])
    for n in range(checkUntil):
        g.advance(0,1)
        if ExtractState([xMin-164, yMin-164, 128, 128], None, sortList = True) == sortedCellsAtStart:
            if periodic:
                dataDict['period'] = str(n+1)
            g.select([xMin-164, yMin-164, 128, 128])
            g.clear(0)
            break
        if n + 1 == checkUntil:
            if inGolly:
                g.select(entryBboxes[0])
                g.fitsel()
                g.exit(f'A catalyst failed to be periodic (of period < {MAX_PERIOD}) or stable. View is centered on that catalyst.')
            else:
                sys.stderr.write(f'A catalyst failed to be periodic (of period < {MAX_PERIOD}) or stable.\n')
                sys.stderr.write(f'Check catalyst with rle {dataDict["rle"]}\n')
                exit(2)

    startAt = 1

    for i in range(1, len(entryBboxes)):

        if i == 1 and len(LifeHistoryToLife(g.getcells(entryBboxes[1]), [4,5])) >= 1:
            relevantCells = LifeHistoryToLife(g.getcells(entryBboxes[1]), [4,5])
            dataDict[f'locus'] = LifeHistoryToRLE(g.getcells(entryBboxes[i]), [4,5])
            startAt = 2
        else:
            relevantCells = LifeHistoryToLife(g.getcells(entryBboxes[i]))
            dataDict[f'forbidden {i-startAt+1}'] = LifeHistoryToRLE(g.getcells(entryBboxes[i]))

        lifeCells = LifeHistoryToLife(g.getcells(entryBboxes[i]))
        relevantX0, relevantY0 = min(relevantCells[::2]), min(relevantCells[1::2])
        lifeX0, lifeY0 = min(lifeCells[::2]), min(lifeCells[1::2])
        lifeCellsWRelevantAtOrigin = copy(lifeCells)
        for j in range(0,len(lifeCellsWRelevantAtOrigin), 2):
            lifeCellsWRelevantAtOrigin[j] -= relevantX0
            lifeCellsWRelevantAtOrigin[j+1] -= relevantY0
        # upper rigth corner of life cells is now at (lifeX0 - relevantX0, lifeY0 - relevantY0)
        translations = sorted(list(product(range(-7, 8), range(-7, 8))), key = lambda x : abs(x[0])+abs(x[1]))
        foundMatch  = False
        for transl in translations:
            shiftedLifeCoords = [(lifeCellsWRelevantAtOrigin[i]+transl[0]+dataDict['dx']+(relevantX0 - lifeX0),
                                            lifeCellsWRelevantAtOrigin[i+1]+transl[1]+dataDict['dy']+(relevantY0-lifeY0))
                                                for i in range(0,len(lifeCellsWRelevantAtOrigin), 2)]
            if all([coord in  shiftedLifeCoords for coord in patCoordsCentered]):
                # correct shift the line up the live cells is: transl[0]+dataDict['dx'], transl[1]+dataDict['dx']
                # however, upper right corner of live cells may be different than upper right corner of locus cells
                if i == 1 and len(LifeHistoryToLife(g.getcells(entryBboxes[1]), [4,5])) >= 1:
                    dataDict['locus dx'] =transl[0]+dataDict['dx'] + (relevantX0 - lifeX0)
                    dataDict['locus dy'] =transl[1]+dataDict['dy'] + (relevantY0 - lifeY0)
                else:
                    dataDict[f'forbid {i-startAt+1} dx'] = transl[0]+dataDict['dx'] + (relevantX0 - lifeX0)
                    dataDict[f'forbid {i-startAt+1} dy'] = transl[1]+dataDict['dy'] + (relevantY0 - lifeY0)
                foundMatch = True
                break
        if not foundMatch:
            msg = "Failed to be able to match up live cells in required or locus to those in the catalyst."
            if inGolly:
                g.select(entryBboxes[i])
                g.fitsel()
                g.exit(msg + " (View is centered on problematic catalyst. Check gaps between rows, entries in row.)")
            else:
                sys.stderr.write(msg+"\n")
                exit(2)
    
    maxForbidden = max(len(entryBboxes) - startAt, maxForbidden)
    
    dataDict['name'] = f'cat {catNumber}'

    catNumber += 1

    allCSVRows.append(copy(dataDict))


        
keys = ['name', 'absence', 'rle', 'dx', 'dy', 'symType']
if periodic:
    keys+=['period']
keys += ['locus', 'locus dx', 'locus dy', 'required', 'req dx', 'req dy']
if not periodic:
    keys += ['antirequired', 'antireq dx', 'antireq dy']

for i in range(1, maxForbidden+1):
    keys += [f'forbidden {i}', f'forbid {i} dx', f'forbid {i} dy']


with open(fname, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames = keys)
    writer.writeheader()
    writer.writerows(allCSVRows)

if not inGolly:
    exit(0)