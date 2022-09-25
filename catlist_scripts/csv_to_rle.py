import golly as g
import itertools
import csv
def AddPadding(cellList):
    if len(cellList) % 2 == 0:
        cellList.append(0)

def ConvertToLifeHistory(cellList, state):
    newCellList = list(itertools.chain(*[[cellList[i], cellList[i+1], state] for i in range(0,len(cellList),2)]))
    AddPadding(newCellList)
    return newCellList

def GetBoundingBox(cellList):
    if g.getrect() == []:
        x0 = 0
        y0 = 0
    else:
        [x0, y0 , _, _] = g.getrect()
    g.putcells(cellList, x0 - 100, y0 -100)
    g.select([x0-150, y0-150, 100, 100])
    g.shrink()
    g.clear(0)
    if g.getselrect() == []:
        return [0,0,0,0]
    return g.getselrect()

# assumed formatting:
# name,absence,rle,x0,y0,sym,requiredRLE,reqX0,reqY0,forbiddenRLE,forbidX0,forbidY0,[more forbid]
# so: 0 for name
#     1 for absence
#     2 for rle
#     3-4 for coord
#     5 for sym
#     6-8 for required (RLE, coord)
#     9+ forbidden (RLE, coord)
# I assume catalysts aren't bigger than 20x20.
file = g.opendialog("Catalyst file as CSV", "csv")
g.new("catlist.rle")
g.setrule("LifeHistory")
with open(file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    firstRow = True
    for data in reader:

        #name = data[0] unused
        #absenceInterval = data[1] unused
        #g.warn(data[2])
        catalyst = ConvertToLifeHistory(g.parse(data['rle']),1)
        catPos = (int(data['dx']), int(data['dy']))
        #sym = data[5] unused
        required = ConvertToLifeHistory( [] if data['required'] == '' else g.parse(data['required']), 4)
        reqPos = [0,0] if data['required'] == '' else [int(data['req dx']), int(data['req dy'])]
        if 'locus' not in data or 'o' not in data['locus']:
            locus = ConvertToLifeHistory([], 1)
            locusPos = [0,0]
        else:
            # g.warn(data['locus'])
            locus = ConvertToLifeHistory(g.parse(data['locus']), 4)
            locusPos = [int(data['locus dx']), int(data['locus dy'])]
        if 'antirequired' not in data or 'o' not in data['antirequired']:
            antirequired = ConvertToLifeHistory([], 4)
            antirequiredPos = [0,0]
        else:
            # g.warn(data['locus'])
            antirequired = ConvertToLifeHistory(g.parse(data['antirequired']), 4)
            antirequiredPos = [int(data['antireq dx']), int(data['antireq dy'])]
        # if len(data) >= 12 and data[10] != '' and data[11] != '':
         #   locusPos = (int(data[10]), int(data[11]))

        # curX and curY are the multiples of 10 that act as the current origin.
        # endOfLastY keeps track of how low last row went
        # endOfLastX keeps track of how far right the last entry in the current row went.
        
        # set default values
        # due to shifting catalysts up, calculate how much extra space we need between this and the last line
        if firstRow:
            curX = 0
            curY = 0
            firstRow = False
        else:
            curX = 0
            forbidY0s = [0 if f'forbid {i} dy' not in data or data[f'forbid {i} dy'] == '' else int(data[f'forbid {i} dy']) for i in range(20)]
            locusY0 = 0 if 'locus dy' not in data or data['locus dy'] == '' else int(data['locus dy'])
            antireqY0 = 0 if 'antireq dy' not in data or data['antireq dy'] == '' else int(data['antireq dy'])
            uppermost = min(min(forbidY0s), catPos[1], locusY0, antireqY0)
            curY = 10*((endOfLastY - uppermost + 12)//10)

        g.putcells(required, curX+reqPos[0], curY+reqPos[1])
        if len(antirequired) > 1:
            g.putcells(antirequired, curX+antirequiredPos[0], curY+antirequiredPos[1])
            [_,_, widthAnti, heightAnti] = GetBoundingBox(antirequired)
        else:
            widthAnti, heightAnti = 0,0

        g.putcells(catalyst, curX+catPos[0], curY+catPos[1], 1,0,0,1,"xor")
        [_,_, widthCat, heightCat] = GetBoundingBox(catalyst)
        endOfLastX = curX+catPos[0]+max(widthCat, widthAnti)
        endOfLastY = curY+catPos[1]+max(heightCat, heightAnti)
        
        
        # catalyst with locus as state 4/5 [if there is locus]
        if len(locus) > 1:
            # need extra room, due to shifting catPos[0] to the left
            curX = 10*((endOfLastX - catPos[0] + 12) // 10)
            # g.warn(str(locus))

            g.putcells(locus, curX+locusPos[0], curY+locusPos[1])
            g.putcells(catalyst, curX+catPos[0], curY+catPos[1], 1,0,0,1,"xor")
            endOfLastX = curX+catPos[0]+width
            endOfLastY = max(endOfLastY, curY+catPos[1]+height)

        # now display forbidden\
        i = 1
        while f'forbid {i} dx' in data and data[f'forbid {i} dx'] != '':
            forbidX, forbidY = int(data[f'forbid {i} dx']), int(data[f'forbid {i} dy'])
            curX = 10*((endOfLastX - forbidX + 12) // 10)
            forbidState = ConvertToLifeHistory(g.parse(data[f'forbidden {i}']), 1)
            [_,_, width, height] = GetBoundingBox(forbidState)
            g.putcells(forbidState, curX+forbidX, curY+forbidY)
            endOfLastX = curX + forbidX + width
            endOfLastY = max(endOfLastY, curY+forbidY+height)
            i += 1
