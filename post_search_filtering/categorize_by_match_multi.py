import sys
from math import ceil
import lifelib
sess = lifelib.load_rules("b3s23")
lt = sess.lifetree(n_layers=1)

zero = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$2o4b2o$o5bo$bo5bo$2o4b2o$o5bo$bo5bo$2o4b2o$o5bo$bo5bo$2o
4b2o$2b2obo$2bob2o!''')
one = lt.pattern('''x = 2, y = 14, rule = B3/S23\n
2o$bo$o$2o2$2o$bo$o$2o2$2o$bo$o$2o!''')
two = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$6b2o$6bo$7bo$6b2o$2b2obo$2bob2o$2o$o$bo$2o$2b2obo$2bob2o!''')
three = lt.pattern('''x = 6, y = 14, rule = B3/S23\n
2obo$ob2o$4b2o$4bo$5bo$4b2o$2obo$ob2o$4b2o$4bo$5bo$4b2o$2obo$ob2o!''')
four = lt.pattern('''x = 7, y = 14, rule = B3/S23\n
2o3b2o$2o3b2o2$2o3b2o$obobobo$2bobo$b2obo$5b2o$6bo$5bo$5b2o$6bo$5bo$5b2o!''')
five = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$2o$o$bo$2o$2b2obo$2bob2o$6b2o$6bo$7bo$6b2o$2b2obo$2bob2o!''')
six = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$2o$o$bo$2o$2b2obo$2bob2o$2o4b2o$o5bo$bo5bo$2o4b2o$2b2obo$2bob2o!''')
seven = lt.pattern('''x = 6, y = 14, rule = B3/S23\n
ob2o$2obo$4b2o$5bo$4bo$4b2o$2b2o$3bo$2bo$2b2o$2o$bo$o$2o!''')
eight = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$2o4b2o$o5bo$bo5bo$2o4b2o$2b2obo$2bob2o$2o4b2o$o5bo$bo5bo$2o4b2o$2b2obo$2bob2o!''')
nine = lt.pattern('''x = 8, y = 14, rule = B3/S23\n
2b2obo$2bob2o$2o4b2o$o5bo$bo5bo$2o4b2o$2b2obo$2bob2o$6b2o$6bo$7bo$6b2o$2b2obo$2bob2o!''')
widths = [8,2,8,6,7,8,8,6,8,8]
digits = {0:zero, 1:one, 2:two, 3:three, 4:four, 5:five, 6:six, 7:seven, 8:eight, 9:nine}

def NumberToLifeState(n):
    numAsState = lt.pattern('')
    xShift = 0
    kerningGap = 3
    for digit in str(n):
        numAsState += digits[int(digit)].shift(xShift,0)
        xShift += widths[int(digit)]+kerningGap
    return numAsState

# this only fixes things split between the two sides. does not fix things entirely on the wrong side.
def FixWraparound(catalysts, zoi):
    
    rightOld = catalysts[30:32, :]
    rightNew = rightOld + (rightOld.shift(-64,0).convolve(zoi) & catalysts).shift(64, 0)
    while (rightNew > rightOld):
        rightOld = rightNew.__copy__()
        rightNew = rightOld + (rightOld.shift(-64,0).convolve(zoi) & catalysts).shift(64, 0)
    leftOld = catalysts[-32:-30, :]
    leftNew = leftOld + (leftOld.shift(64,0).convolve(zoi) & catalysts).shift(-64, 0)
    while (leftNew > leftOld):
        leftOld = leftNew.__copy__()
        leftNew = leftNew + (leftOld.shift(64,0).convolve(zoi) & catalysts).shift(-64, 0)

    bottomOld = catalysts[:,30:32]
    bottomNew = bottomOld + (bottomOld.shift(0,-64).convolve(zoi) & catalysts).shift(0, 64)
    while (bottomNew > bottomOld):
        bottomOld = bottomNew.__copy__()
        bottomNew = bottomNew + (bottomOld.shift(0,-64).convolve(zoi) & catalysts).shift(0, 64)
    topOld = catalysts[:,-32:-30]
    topNew = topOld + (topOld.shift(0,64).convolve(zoi) & catalysts).shift(0, -64)
    while (topNew > topOld):
        topOld = topNew.__copy__()
        topNew = topOld + (topOld.shift(0,64).convolve(zoi) & catalysts).shift(0, -64)

    catalysts += topNew+bottomNew+rightNew+leftNew


def AddMatchesToDict(matchGenToYRowAndArgIndex, argIndex,  genZeroResults, patToMatchRLE, genRange, location, filename):


    # do a bit of prep for pattern matching.
    patToMatch = lt.pattern(patToMatchRLE)
    smallZOI = lt.pattern('')
    smallZOI[-1:2, -1:2] = 1
    deadCells = patToMatch.convolve(smallZOI)
    deadCells -= patToMatch

    allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]
    orientedPats = [] # chosen such the first on cell is at (0,0)
    orientedDeadCells = []
    for op in allOrientations:
        transfPat = patToMatch.transform(op)
        orientedPats.append(transfPat.shift(-1*transfPat.firstcell[0],-1*transfPat.firstcell[1]))
        transfDead = deadCells.transform(op).shift(-1*transfPat.firstcell[0],-1*transfPat.firstcell[1])
        orientedDeadCells.append(transfDead)
    
    if location is None:
        startGenResults = genZeroResults[genRange[0]]
        for y in range(0, genZeroResults.bounding_box[3], 100):
            candidateState = startGenResults[0:64, y:(y+64)]
            done = False
            for g in range(genRange[0], genRange[1]+1):
                for i in range(8):
                    matches = candidateState.match(orientedPats[i], dead = orientedDeadCells[i])
                    if(matches.nonempty()):
                        if g in matchGenToYRowAndArgIndex:
                            matchGenToYRowAndArgIndex[g].append((argIndex, y))
                        elif g not in matchGenToYRowAndArgIndex:
                            matchGenToYRowAndArgIndex[g] = [(argIndex, y)]
                        done = True
                        break
                if done: break
                candidateState = candidateState[1]
    else:
        if location == 'same':
            # figure out what to match where.
            firstResult = genZeroResults[0:64, 0:64]
            for i in range(8):
                matches = firstResult.match(orientedPats[i], dead = orientedDeadCells[i])
                if(matches.nonempty()):
                    aliveToMatch = orientedPats[i]
                    deadToMatch = orientedDeadCells[i]
                    relLocation = matches.firstcell
                    break
                assert(i != 7), f'Could not find match for {filename}'
        elif type(location) is tuple:
            relLocation = (32+location[0], 32+location[1])
            aliveToMatch = orientedPats[0]
            deadToMatch = transfDead[0]

        startGenResults = genZeroResults[genRange[0]]
        for y in range(0, genZeroResults.bounding_box[3], 100):
            candidateState = startGenResults[0:64, y:(y+64)]
            shiftedAliveToMatch = aliveToMatch.shift(relLocation[0],relLocation[1]+y)
            shiftedDeadToMatch = deadToMatch.shift(relLocation[0],relLocation[1]+y)
            for g in range(genRange[0], genRange[1]+1):
                if shiftedAliveToMatch <= candidateState and (candidateState & shiftedDeadToMatch).empty:
                    if g in matchGenToYRowAndArgIndex:
                        matchGenToYRowAndArgIndex[g].append((argIndex, y))
                    elif g not in matchGenToYRowAndArgIndex:
                        matchGenToYRowAndArgIndex[g] = [(argIndex, y)]
                    break
                candidateState = candidateState[1]
            


if __name__ == '__main__':

    assert(len(sys.argv) >= 4)
    # command line input: rleToMatch start-end [matchType] rleFilePath

    # read in everything
    patToMatchRLE = ''
    if sys.argv[1].endswith('.rle'):
        with open(sys.argv[1], 'r') as f:
            for line in f:
                patToMatchRLE += line
    else:
        patToMatchRLE = sys.argv[1]
    patToMatch = lt.pattern(patToMatchRLE)

    timingData = sys.argv[2].split('-')
    startGen, endGen = int(timingData[0]), int(timingData[1])

    location = None
    if sys.argv[3][0] == '(' and sys.argv[3][-1] == ')':
        data = sys.argv[1:-1].split(',')
        location = [int(data[0]), int(data[1])]
    elif sys.argv[3] == 'same':
        location = 'same'

    matchGenToYRowAndArgIndex = {}
    genZeroResultsAll = []
    for argIndex in range(4, len(sys.argv)):
        resultsRLE = ''
        with open(sys.argv[argIndex], 'r') as f:
            for line in f:
                resultsRLE += line
        genZeroResults = lt.pattern(resultsRLE)
        genZeroResults[-4:-2, -4:-2] = 1 # i seem to run into spacing issues when
        # the top left most result doesn't have the top leftmost cell. so put a block
        # in a spot that doesn't matter, so that this is the case.

        # not sure if this is needed.
        #for line in resultsRLE.split('\n'):
        #    if line.startswith('#CXRLE Pos='):
        #        genZeroResults = genZeroResults.shift(*shift)

        genZeroResultsAll.append(genZeroResults)

        AddMatchesToDict(matchGenToYRowAndArgIndex, argIndex, genZeroResults, patToMatchRLE, (startGen, endGen), location,sys.argv[argIndex])

    sortedResults = lt.pattern('')
    sortedKeys = sorted(list(matchGenToYRowAndArgIndex.keys()))
    yNew = 0
    for gen in sortedKeys:
        sortedResults += NumberToLifeState(gen).shift(-50, yNew+32)
        xEndOfLast = -25
        for (argIndex, y) in matchGenToYRowAndArgIndex[gen]:
            rightPat = genZeroResultsAll[argIndex -4]
            categoryResults = rightPat[0:rightPat.bounding_box[2], y:(y+64)]
            xStartNextResult = 100*int(ceil((xEndOfLast+20)//100))
            sortedResults += categoryResults.shift(xStartNextResult, yNew-y)
            xEndOfLast = xStartNextResult + rightPat.bounding_box[2]
        yNew += 100
    
    print(sortedResults.rle_string())






    
    


