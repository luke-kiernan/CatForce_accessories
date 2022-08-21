import golly as g
import itertools
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
with open(file, 'r') as f:
    firstRow = True
    for line in f:
    
        # parse data
        data = line.split(',')

        # check for header row.
        if len(data) >= 3 and not ('o' in data[2]):
            continue

        while(data[-1] == '' or data[-1] == '\n'):
            data.pop()
        assert(len(data) >= 5)

        #name = data[0] unused
        #absenceInterval = data[1] unused
        #g.warn(data[2])
        catalyst = ConvertToLifeHistory(g.parse(data[2]),1)
        catPos = (int(data[3]), int(data[4]))
        #sym = data[5] unused
        required = ConvertToLifeHistory( [] if len(data) < 7 else g.parse(data[6]), 5)
        reqPos = [0,0]
        if len(data) >= 9 and data[7] != '' and data[8] != '':
            reqPos = (int(data[7]), int(data[8]))
        # locus = ConvertToLifeHistory([] if len(data) < 10 else g.parse(data[9]), 3)
        # locusPos = [0,0]
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
            forbidY0s = list(map(lambda x:int(x), data[9::3])) if len(data) >= 13 else [0]
            uppermost = min(min(forbidY0s), catPos[1])
            curY = 10*((endOfLastY - uppermost + 12)//10)

        # catalyst with required as state 5
        g.putcells(catalyst, curX+catPos[0], curY+catPos[1])
        g.putcells(required, curX+reqPos[0], curY+reqPos[1])
        [_,_, width, height] = GetBoundingBox(catalyst)
        endOfLastX = curX+catPos[0]+width
        endOfLastY = curY+catPos[1]+height
        
        
        # catalyst with locus as state 3 [if there is locus]
        #if len(locus) > 1:
            # need extra room, due to shifting catPos[0] to the left
        #    curX = 10*((endOfLastX - catPos[0] + 12) // 10) 
        #    g.putcells(catalyst, curX+catPos[0], curY+catPos[1])
        #    g.putcells(locus, curX+locusPos[0], curY+locusPos[1])
        #    endOfLastX = curX+catPos[0]+width
        #    endOfLastY = max(endOfLastY, curY+catPos[1]+height)

        # now display forbidden
        for n in range(12, len(data),3):
            forbidX, forbidY = int(data[n+1]), int(data[n+2])
            curX = 10*((endOfLastX - forbidX + 12) // 10)
            forbidState = ConvertToLifeHistory(g.parse(data[n]), 1)
            [_,_, width, height] = GetBoundingBox(forbidState)
            g.putcells(forbidState, curX+forbidX, curY+forbidY)
            endOfLastX = curX + forbidX + width
            endOfLastY = max(endOfLastY, curY+forbidY+height)

