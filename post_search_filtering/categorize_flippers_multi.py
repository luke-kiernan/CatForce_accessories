import sys
import lifelib
import os
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
allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]

def ComputeOrientationsAndDeadCells(pat):
    orientedPats = []
    orientedDeadCells = []
    smallZOI = lt.pattern('')
    smallZOI[-1:2, -1:2] = 1
    deadCells = pat.convolve(smallZOI)
    deadCells -= pat
    for op in allOrientations:
        transfPat = pat.transform(op)
        orientedPats.append(transfPat)
        transfDeadCells = deadCells.transform(op)
        orientedDeadCells.append(transfDeadCells)
    return orientedPats, orientedDeadCells


def NumberToLifeState(n):
    numAsState = lt.pattern('')
    xShift = 0
    kerningGap = 3
    for digit in str(n):
        numAsState += digits[int(digit)].shift(xShift,0)
        xShift += widths[int(digit)]+kerningGap
    return numAsState

def ChangeFirstCellToOn(rleString):
    startCellData = rleString.index('rule = B3/S23\n') + len('rule = B3/S23\n')
    ind =  startCellData
    multiplicity = ''
    while(rleString[ind].isnumeric()):
        multiplicity += rleString[ind]
        ind += 1
    newMultiplicity = 0 if multiplicity == '' else str(int(multiplicity)-1)
    if rleString[ind] == 'o':
        return rleString, False
    else:
        firstCellOnRLE = rleString[:startCellData]+f'o{newMultiplicity}{rleString[ind]}'+rleString[ind+1:]
        return firstCellOnRLE, True

def FindActivePattern(firstResult, rleToMatch):
    firstResultCopy = firstResult + lt.pattern('')
    patToMatch = lt.pattern(rleToMatch)
    orientedPats, orientedDeadCells = ComputeOrientationsAndDeadCells(patToMatch)
    activePat = lt.pattern('')
    for i in range(8):
        matches = firstResultCopy.match(orientedPats[i], dead = orientedDeadCells[i])
        if(matches.nonempty()):
            firstCell = matches.firstcell
            shouldBeContained = orientedPats[i].shift(*firstCell)
            assert(shouldBeContained <= firstResultCopy), i
            firstResultCopy -= orientedPats[i].shift(*firstCell)
            matches[firstCell[0], firstCell[1]] = 0
            activePat += orientedPats[i].shift(*firstCell)
    assert(activePat <= firstResult)
    return activePat

def CheckHowFlips(result, genRange, orientedActive, orientedActiveDead):
    workspace = result + lt.pattern('')
    for g in range(0, genRange[1]+1):
        if g >= genRange[0]:
            for i in range(8):
                for toShift in [(0,0), (-1,-1), (-1, 0), (0, -1)]:
                    if orientedActive[i].shift(*toShift) <= workspace \
                                and (orientedActiveDead[i].shift(*toShift) & workspace).empty():
                    # in rare cases, there could be multiple different transformations
                    #  that line up the active regions correctly...
                    # also, there's some redundancy here: some of these shifts
                    # don't work with certain symmetries, and certain op's end up being the same...
                    
                        return g, allOrientations[i], toShift
        workspace = workspace[1]
    return -1, None, None

if __name__ == '__main__':

    # read in everything
    patToMatchRLE = ''
    if sys.argv[1].endswith('.rle'):
        with open(sys.argv[1], 'r') as f:
            for line in f:
                patToMatchRLE += line
    else:
        patToMatchRLE = sys.argv[1]

    timingData = sys.argv[2].split('-')
    startGen, endGen = int(timingData[0]), int(timingData[1])

    complete = {}
    maybeWeldable = {}
    genZeroResultsAll = []
    activeRegionsAll = []
    counts = {'same': 0, 'parse error': 0, 'weldable': 0, 'flips ok': 0, 'barely failed':0, 'hopeless': 0}
    for argIndex in range(3, len(sys.argv)):
        resultsRLE = ''
        with open(sys.argv[argIndex], 'r') as f:
            for line in f:
                resultsRLE += line
        # manual edit to RLE it so lifelib doesn't automatically shift pattern.
        rleWithFirstCellOn, changeBack = ChangeFirstCellToOn(resultsRLE)
        # correct for extended rle
        genZeroResults = lt.pattern(rleWithFirstCellOn)
        shift = (0,0)
        for line in resultsRLE.split('\n'):
            if line.startswith('#CXRLE Pos='):
                shiftData  = line[line.index('=')+1:].split(',')
                shift = (int(shiftData[0]), int(shiftData[1]) )
                genZeroResults = genZeroResults.shift(*shift)
                break
            if line.endswith('rule = B3/S23'):
                break

        # undo change from manual edit and add block as a better solution
        genZeroResults[-4:-2, -4:-2] = 1
        if changeBack: 
            genZeroResults[shift[0], shift[1]] = 0
            assert(genZeroResults[shift[0], shift[1]] == 0)

        genZeroResultsAll.append(genZeroResults)

        # figure out active pattern
        firstResult = genZeroResults[0:64,0:64].shift(-32,-32)
        activePat = FindActivePattern(firstResult, patToMatchRLE)
        assert(activePat <= firstResult)
        activeRegionsAll.append(activePat)
        orientedActive, orientedActiveDead = ComputeOrientationsAndDeadCells(activePat)
        activeDeadCells = orientedActiveDead[0]
        for y in range(0, genZeroResults.bounding_box[3], 100):
            category = genZeroResults[0:genZeroResults.bounding_box[2], y:(y+64)].shift(0, -y-32)
            firstResult = category[0:64, -100:100].shift(-32, 0)
            assert(activePat <= firstResult)
            g, op, translation = CheckHowFlips(firstResult, (startGen, endGen), orientedActive, orientedActiveDead)
            if op != "identity" and g > 0:
                flippedActive = activePat.transform(op).shift(*translation)
                assert(flippedActive <= firstResult[g])
                for x in range(0, category.bounding_box[2], 100):
                    result = category[x:(x+64), -32:32].shift(-x-32, 0)
                    assert(activePat <= result)
                    #assert(flippedActive <= result[g]), f'{flippedActive.rle_string()}\n{result[g].rle_string()}'
                    catalysts = result - activePat
                    history = lt.pattern('')
                    for i in range(g+1):
                        history += result[i]
                    flippedCatalysts = catalysts.transform(op).shift(*translation)
                    if (flippedCatalysts & (history - catalysts)).empty(): 
                        combineThenStep = catalysts + flippedCatalysts
                        combineThenStep = combineThenStep[1]
                        stepThenCombine = catalysts[1] + flippedCatalysts[1]
                        if combineThenStep != stepThenCombine:
                             # flipped catalysts overlap with existing catalysts, but not with active region.
                            flipped = (catalysts + flippedCatalysts+activePat)[g]
                            flipped -= catalysts # this assumes all the catalysts have recovered, which might not be the case.
                            flipped += flippedCatalysts
                            if activePat <= flipped[g] and (flipped[g] & activeDeadCells).empty():
                                if 2*g in maybeWeldable:
                                    maybeWeldable[2*g].append((argIndex, x,y))
                                else:
                                    maybeWeldable[2*g] = [(argIndex, x,y)]
                                counts['weldable'] += 1
                            else: # not quite right...
                                counts['barely failed'] += 1
                        else:
                            # no interference.
                            combined = result + flippedCatalysts
                            if activePat <= combined[2*g] and (combined[2*g] & activeDeadCells).empty():
                                if 2*g in complete:
                                    complete[2*g].append((argIndex, x,y, op, translation))
                                else:
                                    complete[2*g] = [(argIndex, x,y, op, translation)]
                                counts['flips ok'] += 1
                            else: # not quite right...
                                counts['barely failed'] += 1
                    else:
                        counts['hopeless'] += 1
            # include non-flippers too.
            elif op == 'identity' and translation == (0,0):
                if g in complete:
                    complete[g].append((argIndex, 0,y, op, translation))
                else:
                    complete[g] = [(argIndex, 0,y, op, translation)]
                counts['same'] += 1
            else:
                counts['parse error'] += 1

    completePat = lt.pattern('')
    yNew = 32
    for g in sorted(complete.keys()):
        completePat += NumberToLifeState(g).shift(-50, yNew)
        xNew = 32
        for argIndex, x,y,op,transl in complete[g]:
            
            if op == 'identity':
                category = (genZeroResultsAll[argIndex-3][:, y:(y+64)]).shift(0, -y-32)
                completePat += category.shift(xNew, yNew)
                xNew += 100*( (category.bounding_box[2]+100) // 100)
            else:
                result = (genZeroResultsAll[argIndex-3][x:(x+64), y:(y+64)]).shift(-x-32, -y-32)
                catalysts = result - activeRegionsAll[argIndex-3]
                flippedCatalysts  = catalysts.transform(op).shift(*transl)
                combined = result + flippedCatalysts
                completePat += combined.shift(xNew, yNew)
                xNew += 100
        yNew += 100
    print(completePat.rle_string())
    sys.stderr.write(f'total categories processed: {sum([pair[1] for pair in counts.items()])}'+os.linesep)
    sys.stderr.write(f'\tparse error: {counts["parse error"]}'+os.linesep)
    sys.stderr.write(f'\tsame spot: {counts["same"]}'+os.linesep)
    sys.stderr.write(f'\tflipped: {counts["weldable"]+counts["hopeless"]+counts["barely failed"]+counts["flips ok"]}, with subtotals'+os.linesep)
    sys.stderr.write(f'\t\thopeless:{counts["hopeless"]}'+os.linesep)
    sys.stderr.write(f'\t\tbarely failed:{counts["barely failed"]}'+os.linesep)
    sys.stderr.write(f'\t\tmaybe weldable:{counts["weldable"]}'+os.linesep)
    sys.stderr.write(f'\t\tworks:{counts["flips ok"]}'+os.linesep)