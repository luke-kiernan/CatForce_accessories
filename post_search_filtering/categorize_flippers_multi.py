from collections import Counter
import sys
import lifelib
import os
from categorize_by_match_multi import NumberToLifeState, FixWraparound
from math import ceil
sess = lifelib.load_rules("b3s23")
lt = sess.lifetree(n_layers=1)

allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]

def HasSmallPrimeFactor(n):
    smallPrimes= [2,3,5,7,11,13]
    if any([n % prime == 0 for prime in smallPrimes]):
        return True
    return False

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

def FindActivePattern(firstResult, rleToMatch):
    firstResultCopy = firstResult.__copy__()
    patToMatch = lt.pattern(rleToMatch)
    orientedPats, orientedDeadCells = ComputeOrientationsAndDeadCells(patToMatch)
    activePat = lt.pattern('')
    for i in range(8):
        matches = firstResultCopy.match(orientedPats[i], dead = orientedDeadCells[i])
        if(matches.nonempty()):
            firstCell = matches.firstcell
            shouldBeContained = orientedPats[i].shift(*firstCell)
            assert(shouldBeContained <= firstResultCopy)
            firstResultCopy -= orientedPats[i].shift(*firstCell)
            matches[firstCell[0], firstCell[1]] = 0
            activePat += orientedPats[i].shift(*firstCell)
    assert(activePat <= firstResult)
    return activePat

def Nontrivial(result, activeRegion, lastGen):
    '''Returns True if the different images of the active region interact, False if they do not.'''
    comps = result.components()
    workspace = result[3]
    gen = 0
    activeRegionComps = [i for i in range(len(comps)) if (comps[i]&activeRegion).nonempty()]
    increment = 3
    while(len(set(activeRegionComps)) > 1 and gen < lastGen):
        toDelete = set()
        for i in range(len(comps)):
            for j in range(i+1, len(comps)):
                if workspace.component_containing(comps[i]) == workspace.component_containing(comps[j]) \
                        and workspace.component_containing(comps[i]).nonempty():
                    comps[i] += comps[j]
                    comps[j] += comps[i]
                    toDelete.add(j)
                    if j in activeRegionComps:
                        activeRegionComps[activeRegionComps.index(j)] = i
                        
        for ind in sorted(list(toDelete), reverse=True):
            del comps[ind]
        comps = [comp[increment] for comp in comps]
        workspace = workspace[increment]
        gen += increment
    if len(set(activeRegionComps)) == 1:
        return True
    return False

def WrapEdges(result):
    result += result[:, -32:-22].shift(0,64)
    result += result[-32:-22, :].shift(64,0)
    result += result[:, 22:32].shift(0,-64)
    result += result[22:32, :].shift(-64,0)

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
    sortBy = "time"
    startAt = 3
    if sys.argv[3] == "junk":
        sortBy = "junk"
        startAt = 4
    startGen, endGen = int(timingData[0]), int(timingData[1])

    results = {}
    could_not_parse = 0
    offending_files = set()
    flippers_found = 0
    thrown_out = 0
    categoryBreakdown = Counter()

    smallZOI = lt.pattern('')
    smallZOI[-1:2, -1:2] = 1
    normalZOI = lt.pattern('x = 5, y = 5, rule = B3/S23\nb3o$5o$5o$5o$b3o!').shift(-2,-2)
    bigZOI = lt.pattern('')
    bigZOI[-3:4, -3:4] = 1
    
    maxNumberToTakeFromCategory = 10

    for argIndex in range(startAt, len(sys.argv)):
        resultsRLE = ''
        with open(sys.argv[argIndex], 'r') as f:
            for line in f:
                resultsRLE += line
        genZeroResults = lt.pattern(resultsRLE)
        genZeroResults[-4:-2, -4:-2] = 1 # i seem to run into spacing issues when
        # the top left most result doesn't have the top leftmost cell. so put a block
        # in a spot that doesn't matter, so that this is the case.

        # not sure if this is needed.
        for line in resultsRLE.split('\n'):
            if line.startswith('#CXRLE Pos='):
                data = line[line.index('#CXRLE Pos=')+len('#CXRLE Pos='):].split(',')
                shift = [int(data[0]), int(data[1])]
                genZeroResults = genZeroResults.shift(*shift)

        # figure out active pattern
        firstResult = genZeroResults[0:64,0:64].shift(-32,-32)
        activePat = FindActivePattern(firstResult, patToMatchRLE)
        catalysts = firstResult - activePat
        periodic = catalysts == catalysts[1]
        orientedActive, orientedActiveDead = ComputeOrientationsAndDeadCells(activePat)
        activeDeadCells = orientedActiveDead[0]
        for y in range(0, genZeroResults.bounding_box[1]+genZeroResults.bounding_box[3], 100):
            scoredCategoryResults = []

            # figure out if/how it flips.
            category = genZeroResults[0:genZeroResults.bounding_box[2], y:(y+64)].shift(0, -y-32)
            x0 = 0
            g = -1
            testResult = category[x0:(x0+64), -32:32].shift(-x0-32, 0)
            while(g == -1 and testResult.nonempty()):
                g, op, translation = CheckHowFlips(testResult, (startGen, endGen), orientedActive, orientedActiveDead)
                if g == -1:
                    # attempt to fix wrap-around issues.
                    WrapEdges(testResult)
                    FixWraparound(testResult, normalZOI)
                    g, op, translation = CheckHowFlips(testResult, (startGen, endGen), orientedActive, orientedActiveDead)
                if g == -1:
                    # no luck. try the next result in the category, maybe it's better.
                    x0 += 100
                    testResult = category[x0:(x0+64), -32:32].shift(-x0-32, 0)
            
            if testResult.empty():
                categoryBreakdown.update(['could not parse'])
                offending_files.add(sys.argv[argIndex])
                continue # all results in category were "bad"
            if not Nontrivial(testResult,activePat,g):
                categoryBreakdown.update(['trivial'])
                continue
            if op != "identity":
                categoryBreakdown.update(['flipped'])
            else:
                categoryBreakdown.update(['same spot'])

            flippedActivePat = activePat.transform(op).shift(*translation)
            
            # the rest of the ones in the row flip the same way, which saves us some trouble.
            scores = []
            for x in range(0, category.bounding_box[0]+category.bounding_box[2], 100):
                result = category[x:(x+64), -32:32].shift(-x-32, 0)
                catalysts = result - activePat
                
                if not (flippedActivePat <= result[g]):
                    WrapEdges(result)
                    FixWraparound(result, normalZOI)
                if not (flippedActivePat <= result[g]):
                    continue

                # TODO: if periodic, check if it interacts more than once,
                # and if the period is one for which we have some sparkers.

                if op != "identity":
                    history = lt.pattern('')
                    for i in range(g+1):
                        history = history + result[i]
                    flippedCatalysts = catalysts.transform(op).shift(*translation)
                    if not periodic:
                        if (flippedCatalysts & (history-catalysts)).nonempty():
                            thrown_out += 1
                            # flipped catalysts overlap with cells that were active. probably worthless
                            continue
                    else:
                        catHistory = lt.pattern('')
                        for i in range(g+1):
                            catHistory = catHistory + catalysts[i]
                        if (flippedCatalysts & (history-catHistory)).population > 12:
                            # 12 here more or less arbitrary: oscs can dodge some overlap via 
                            # fortunate timing but not a lot.
                            thrown_out += 1
                            continue
                    # TODO: check for catalyst-flipped catalyst interference
                    # maybe it doesn't work, but it's plausibly weldable.

                    testing = testResult[g]+flippedCatalysts
                    if activePat <= testing[g]:
                        # success!!
                        flippers_found += 1
                        result += flippedCatalysts[-g]
                        g *= 2
                        op = "identity"
                        translation = (0,0)
                        flippedActivePat = activePat.__copy__()
                        flippedCatalysts = catalysts.__copy__()
                        flippers_found += 1
                
                score = 0
                phase = result[g]
                if sortBy == "junk":
                    score = (phase-flippedActivePat-catalysts).population
                    cutoff = -1
                else:
                    flippedActiveEvolved = flippedActivePat.__copy__()
                    flippedActiveDead = flippedActiveEvolved.convolve(smallZOI)-flippedActiveEvolved
                    # TODO: helper function that computes generation of catalyst destruction
                    # TODO: test that the below calculate k, k2 correctly
                    #       switch to using k2 only?
                    k = 0
                    # check for how long the active region evolves "correctly"
                    while(flippedActiveEvolved <= phase[k] and (flippedActiveDead & phase[k]).empty() and k < 2*g):
                        k += 1
                        flippedActiveEvolved = flippedActiveEvolved[1]
                        flippedActiveDead = flippedActiveEvolved.convolve(smallZOI)-flippedActiveEvolved
                    
                    # we count the number of generations for which the active region evolves
                    #  in the same way as it did g gens ago
                    k2 = 0
                    kthGenActivePart = (result[k2]-(result[k2] & catalysts[k2].convolve(smallZOI)))\
                                                                .transform(op).shift(*translation)
                    kthGenHalo = kthGenActivePart.convolve(smallZOI)-kthGenActivePart
                    matchedGen = result[g+k2]
                    while( (kthGenHalo & matchedGen).empty() and kthGenActivePart <= matchedGen and k2 < 2*g):
                        k2 += 1
                        kthGenActivePart = (result[k2]-(result[k2] & catalysts[k2].convolve(smallZOI)))\
                                                                .transform(op).shift(*translation)
                        kthGenHalo = kthGenActivePart.convolve(smallZOI)-kthGenActivePart
                        matchedGen = result[g+k2]

                    k = max(k, k2)
                    if (k < 1):
                        sys.stderr.write(f"problem in {sys.argv[argIndex]} near ({x},{y})"+os.linesep)
                        sys.stderr.write(kthGenActivePart.rle_string())
                        sys.stderr.flush()
                    assert(k >= 1)
                    score = k - 1
                scoredCategoryResults.append((score,result.__copy__()))
            if len(scoredCategoryResults) > 0:
                scoredCategoryResults.sort(key = lambda pair : -1*pair[0])
                if not g in results:
                    results[g] = []
                results[g] += scoredCategoryResults[:maxNumberToTakeFromCategory]

    for _, listOfResults in results.items():
        listOfResults.sort(key = lambda pair : -1*pair[0])
    sortedResults = lt.pattern('')
    sortedKeys = sorted(results.keys())
    yNew = 0
    for gen in sortedKeys:
        sortedResults += NumberToLifeState(gen).shift(-50, yNew+32)
        xEndOfLast = -25
        for _, result in results[gen]:
            xStartNextResult = 100*int(ceil((xEndOfLast+20)/100))
            sortedResults += result.shift(xStartNextResult+32, yNew+32)
            xEndOfLast = xStartNextResult + result.bounding_box[2]
        yNew += 100
    print(sortedResults.rle_string())
    sys.stderr.write('Category breakdown:'+os.linesep)
    for key, val in dict(categoryBreakdown).items():
        sys.stderr.write(f'\t{key}: {val}' + os.linesep)
    
    sys.stderr.write('Results:'+os.linesep)
    sys.stderr.write(f'\t{thrown_out} results where flipped catalysts overlapped with active'+os.linesep)
    if flippers_found > 0:
        sys.stderr.write(f'\t{flippers_found} results where flipped placements worked'+os.linesep)