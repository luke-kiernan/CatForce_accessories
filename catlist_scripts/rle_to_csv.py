from math import ceil
from copy import copy
import golly as g
import os
import csv
from itertools import product

g.note('Program for producing lists for periodic CatForce.')
g.note('1 catalyst per row. \n First entry: sparker/catalyst, with state 3/4/5 for required cells. \n Rest are forbidden.')
g.note('Include at least 5 cells between rows, and 5 cells between entries in the in the row.\n '+
        'Patterns don\'t need to exactly line up in each row, but they shouldn\'t be off by more than 7.')
g.note('Center is done automatically, but period, symmetry, and max absence time must be entered by hand.')
# states 3/4/5 = required (matching) cells

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

minRowGap = 4
minEntryGap = 4

# list of list of rectangle lists
# each list is the entries in 1 row.
rectanglesInRow = []
[xMin, yMin , w, h ] = g.getrect()

y = yMin
rowStartY = yMin
lastNonemptyRow = yMin
streakOfEmptyRows = 0

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

g.warn("analyzing entries in each row...")

for entryBboxes in rectanglesInRow:
    dataDict = {}
    
    stateRLE = LifeHistoryToRLE(g.getcells(entryBboxes[0]))
    
    patCells = LifeHistoryToLife(g.getcells(entryBboxes[0]))
    x0Pat, y0Pat = min(patCells[::2]), min(patCells[1::2])
    w, h = max(patCells[::2])-x0Pat+1,max(patCells[1::2])-y0Pat+1


    center = (ceil(x0Pat+w/2)//1, ceil(y0Pat+h/2)//1)
    # catForce universe [post-(dx,dy)] catalyst
    #          = golly universe catalyst, shifted by -center.

    dataDict['rle'] = stateRLE
    dataDict['dx'] = -(center[0]-x0Pat)
    dataDict['dy'] = -(center[1]-y0Pat)

    
    reqCells = LifeHistoryToLife(g.getcells(entryBboxes[0]), [3,4,5])
    if len(reqCells) >= 1:
        dataDict['required'] = LifeHistoryToRLE(g.getcells(entryBboxes[0]), [3,4,5])
        x0Req, y0Req = min(reqCells[::2]), min(reqCells[1::2])
        dataDict['req dx'] = x0Req - center[0]
        dataDict['req dy'] = y0Req - center[1]
    
    # want (dxReq, dyReq) such that
    # default position required shifted by (dxReq, dyReq)
    #  equals golly universe required, shifted by -center.
    # that is, (x0Req, y0Req) - center = (dxReq, dyReq)

    patCellsCentered = g.transform(g.getcells(entryBboxes[0]), -center[0], -center[1])
    patLifeCellsCentered = LifeHistoryToLife(patCellsCentered)
    patCoordsCentered = [tuple(patLifeCellsCentered[i:i+2]) for i in range(0, len(patLifeCellsCentered), 2)]

    dataDict['period'] = ''

    maxForbidden = max(len(entryBboxes) - 1, maxForbidden)

    for i in range(1, len(entryBboxes)):
        
        dataDict[f'forbidden {i}'] = LifeHistoryToRLE(g.getcells(entryBboxes[i]))

        forbiddenCells = g.getcells(entryBboxes[i])
        forbiddenLifeCells = LifeHistoryToLife(forbiddenCells)
        forbidX0, forbidY0 = min(forbiddenLifeCells[::2]), min(forbiddenLifeCells[1::2])
        forbiddenLifeCellsAtOrigin = copy(forbiddenLifeCells)
        for j in range(0,len(forbiddenLifeCellsAtOrigin), 2):
            forbiddenLifeCellsAtOrigin[j] -= forbidX0
            forbiddenLifeCellsAtOrigin[j+1] -= forbidY0
        

        translations = sorted(list(product(range(-7, 8), range(-7, 8))), key = lambda x : abs(x[0])+abs(x[1]))
        for transl in translations:
            shiftedForbiddenLifeCoords = [(forbiddenLifeCellsAtOrigin[i]+transl[0]+dataDict['dx'],
                                            forbiddenLifeCellsAtOrigin[i+1]+transl[1]+dataDict['dy'])
                                                for i in range(0,len(forbiddenLifeCellsAtOrigin), 2)]
            if all([coord in shiftedForbiddenLifeCoords for coord in patCoordsCentered]):
                dataDict[f'forbid {i} dx'] = transl[0]+dataDict['dx']
                dataDict[f'forbid {i} dy'] = transl[1]+dataDict['dy']
                break
    
    dataDict['name'] = f'cat {catNumber}'
    dataDict['symType'] = ''
    dataDict['absence'] = ''
    catNumber += 1

    allCSVRows.append(copy(dataDict))


        
keys = ['name', 'absence', 'rle', 'dx', 'dy', 'symType', 'period', 'required', 'req dx', 'req dy']
for i in range(1, maxForbidden+1):
    keys += [f'forbidden {i}', f'forbid {i} dx', f'forbid {i} dy']

fname = g.savedialog("Save csv catalyst list", "*.csv")
if fname == '':
    g.exit("File save canceled.")
with open(fname, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames = keys)
    writer.writeheader()
    writer.writerows(allCSVRows)