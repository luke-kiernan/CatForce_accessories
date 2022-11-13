# arguments: sampleSoupFile outputFile
# input file formatting assumed to be
# [catalyst rle]:
# \t[soup1]
# \t[soup2]
import sys
import random
from statistics import mode
from math import floor,ceil
import lifelib
import csv
sess = lifelib.load_rules("b3s23")
lt = sess.lifetree(n_layers=1)
allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]

inverse = {}
for orient in allOrientations:
    if orient in ["rot270", "rot90"]:
        otherInd = (["rot270", "rot90"].index(orient) + 1) % 2
        inverse[orient] = ["rot270", "rot90"][otherInd]
    else:
        inverse[orient] = orient

def Similarity(pat1, pat2):
    if (pat1 + pat2).population > 0:
        return (pat1&pat2).population / (pat1 + pat2).population
    return 1

def TilePat(pat, bbox=[-64,-64,3*64,3*64]):
    tiled = lt.pattern('')
    xStart = 64*floor(bbox[0]/64) if bbox[0] < 0 else 64*ceil(bbox[0]/64)
    xStop = 64*floor((bbox[0]+bbox[2])/64) if bbox[0]+bbox[2] < 0 else 64*ceil((bbox[0]+bbox[2])/64)
    yStart = 64*floor(bbox[1]/64) if bbox[1] < 0 else 64*ceil(bbox[1]/64)
    yStop = 64*floor((bbox[1]+bbox[3])/64) if bbox[1]+bbox[3] < 0 else 64*ceil((bbox[1]+bbox[3])/64)
    intersectWith = lt.pattern('')
    intersectWith[bbox[0]:(bbox[0]+bbox[2]), bbox[1]:(bbox[1]+bbox[3])] = 1
    for i in range(xStart, xStop+64, 64):
        for j in range(yStart, yStop+64, 64):
            tiled+= pat(i,j) & intersectWith
    return tiled

def ComputeOrientationsAndDeadCells(pat):
    orientedPats = []
    orientedDeadCells = []
    smallZOI = lt.pattern('')
    smallZOI[-1:2, -1:2] = 1
    for op in allOrientations:
        transfPat = lt.pattern('')
        transfPat |= pat.transform(op)
        orientedPats.append(transfPat)
        orientedDeadCells.append(transfPat.convolve(smallZOI)-transfPat)
    return orientedPats, orientedDeadCells

def TrimRLE(rleString):
    trimmed = ''
    for line in rleString.split('\n'):
        if line.startswith('x') or line.startswith('#'):
            pass
        elif 'o' in line or 'b' in line or '$' in line:
            trimmed += line.rstrip().lstrip()
    return trimmed

def LocateCatalyst(catPat, reactionPat):
    orientedPats, orientedDeadCells = ComputeOrientationsAndDeadCells(catPat)
    
    bestOrientation = None
    closest = [None, None]
    resultPat = reactionPat[:,:]
    # shouldBeContained = lt.pattern('$'.join(64*['64o'])).shift(-32,-32)
    for tryNum in range(2):
        for i in range(8):
            matches = resultPat.match(orientedPats[i], dead = orientedDeadCells[i])
            if(matches.nonempty()):
                bestOrientation = i
                shouldBeContained = orientedPats[i].shift(*matches.firstcell)
                closest = (matches.firstcell[0], matches.firstcell[1])
                break
        if tryNum == 0 and closest[0] is None:
            resultPat = lt.pattern('')
            for x in [-64, 0,64]:
                for y in [-64, 0,64]:
                    resultPat += reactionPat.shift(x,y) & lt.pattern('$'.join(64*['64o'])).shift(-32,-32)
    
    #assert(shouldBeContained <= resultPat and not (bestOrientation is None))
    return closest, allOrientations[bestOrientation]

def GetOffset(otherPat, origin):
    if otherPat.empty():
        return None
    upperRightCorner = otherPat.bounding_box[:2]
    return (upperRightCorner[0]-origin[0], upperRightCorner[1]-origin[1])

def ComputeStabilizer(pat, patHalo, stabilizedBy, shiftToMatchUp):
    for transf in allOrientations:
        if (pat.transform(transf).match(pat, patHalo)).nonempty():
            stabilizedBy.add(transf)
            shiftToMatchUp[transf] = pat.match(pat.transform(transf),\
                                                patHalo.transform(transf)).firstcell

# read in sample soups
x = 0
y = 0
catsToReactions = {}
currentCat = None
with open(sys.argv[1], "r") as f:
    for line in f:
        if ':' in line and 'o' in line.split(':')[0]:
            catsToReactions[line.split(':')[0]] = []
            currentCat = line.split(':')[0]
        elif 'o' in line.rstrip().lstrip() and not (currentCat is None):
            catsToReactions[currentCat].append(line.rstrip().lstrip())
# econ catalysts
econCats = []
with open('/Users/lukekiernan/Downloads/eonomical_list.txt', 'r') as f:
    for line in f:
        econCats.append(lt.pattern(line.rstrip()))
keys = ['name', 'absence','rle','dx','dy', 'symType']
rows = []
doneWith  = 0

# some constants
smallZOI = lt.pattern('')
smallZOI[-1:2, -1:2] = 1
cross = lt.pattern('')
cross[0, -1:2] = 1
cross[ -1:2,0] = 1
ZOI = cross.convolve(smallZOI)
snake = lt.pattern('2o$o$bo$2o')
snakeHalo = snake.convolve(smallZOI) - snake

for cat in catsToReactions:
    absence = 5

    catPat = lt.pattern(cat)
    catbbox = catPat.bounding_box
    assert(catbbox[3] < 64 and catbbox[2] < 64)
    

    
    catHalo = catPat.convolve(smallZOI) - catPat
    # econ list
    inEcon = False
    for orient in allOrientations:
        for econCat in econCats:
            if econCat.population == catPat.population and\
                 (econCat.match(catPat.transform(orient), catHalo.transform(orient))).nonempty():
                inEcon = True
                break
        if inEcon: break
    if not inEcon:
        sys.stderr.write('  skipping cat ' + str(doneWith+1) + ' since not in econ\n')
        doneWith += 1
        continue
        

    startAntiReqPat = catPat.convolve(smallZOI)
    catbbox = catPat.bounding_box
    bigBox = lt.pattern('')
    bigBox[(catbbox[0]-64):(catbbox[0]+catbbox[2]+64),
           (catbbox[1]-64):(catbbox[1]+catbbox[3]+64) ] = 1
    comps = (bigBox - startAntiReqPat).components()
    for comp in comps:
        if comp.population < 8:
            startAntiReqPat += comp
    startAntiReqPat -= catPat

    mooreComps = catPat.components(ZOI)
    
    stabilizedBy = set()
    shiftToMatchUp = {}
    ComputeStabilizer(catPat, catHalo, stabilizedBy, shiftToMatchUp)
    
    comps =  catPat.components()
    isSnake = [any([comp.match(snake.transform(transf), snakeHalo.transform(transf))
                 for transf in allOrientations]) for comp in comps]
    noSymLocus = len(stabilizedBy) > 1 and len(comps) > 1 and \
                    all([comp.population < 7 for comp in comps]) and \
                        not all(isSnake)
    # noSymLocus here means: don't symmetricize the locus/required/antirequired
    # (guessing it's a symmetric catalyst that acts in an asymmetric way)
    # (snakes might be an exception here, since those aren't transparent.)
    
    reqTries = []
    antiReqTries = []
    locusTries = []

    skippedSoups = []
    for tryNum in range(5):
        compTransparentCount = [0 for _ in mooreComps]
        skippedSoups.append(set())    
        transpLocus = lt.pattern('')
        transpAntiReq= lt.pattern('')
        transpAntiReq |= startAntiReqPat

        antiReqPat = lt.pattern('')
        antiReqPat |= startAntiReqPat

        reqPat = lt.pattern('')
        reqPat |= catPat

        locusPat = lt.pattern('')
        if tryNum > 0:
            random.shuffle(catsToReactions[cat])
        for soup in catsToReactions[cat]:
            soupPat = lt.pattern(soup)
            soupBbox = soupPat.bounding_box
            location, orient = LocateCatalyst(catPat, soupPat)
            centeredCatSoup = TilePat(soupPat.shift(-1*location[0],
                    -1*location[1]).transform(inverse[orient]),
                    [-4*64,-4*64,5*64,5*64])
            # assert(catPat <= centeredCatSoup)
            # assert((centeredCatSoup & catHalo).empty())
            soupCenteredNoCat = centeredCatSoup - catPat

            # we update these each gen.
            reqForRun = lt.pattern('')
            reqForRun |= catPat
            antiReqForRun = lt.pattern('')
            antiReqForRun |= startAntiReqPat
            locusForRun = lt.pattern('')
            lastActivated = -1
            absentFor = 0
            recoveredFor = 0

            # we only update these once it's recovered
            savedReqForRun = lt.pattern('')
            savedReqForRun |= catPat
            savedAntiReqForRun = lt.pattern('')
            savedAntiReqForRun |= startAntiReqPat
            maxAbsentFor = 0

            for n in range(200):
                centeredCatSoup = centeredCatSoup[1]
                catPresent = (centeredCatSoup & catHalo).empty() and catPat <= centeredCatSoup
                if not catPresent:
                    if lastActivated == -1:
                        locusForRun += (centeredCatSoup & catHalo).convolve(smallZOI) & catPat
                        lastActivated = n
                    elif recoveredFor > 0:
                        lastActivated = n
                    antiReqForRun -= centeredCatSoup
                    reqForRun &= centeredCatSoup
                    recoveredFor = 0
                    absentFor += 1
                elif lastActivated != -1: # present and has been activated.
                    recoveredFor += 1
                    if recoveredFor == 5: 
                        savedReqForRun &= reqForRun
                        savedAntiReqForRun &= antiReqForRun
                        maxAbsentFor = max(maxAbsentFor, n-4-lastActivated)
                else: # in case the catalyst prevents adjacent births
                    # but doesn't have cells born next to it (so never absent)
                    soupCenteredNoCat = soupCenteredNoCat[1]
                    if (soupCenteredNoCat & catHalo).nonempty():
                        lastActivated = n
                        locusForRun += (soupCenteredNoCat & catHalo).convolve(smallZOI) & catPat
                if absentFor > 30 or recoveredFor > 15:
                    break
            if maxAbsentFor > 0:
                absence = max(absence, maxAbsentFor)
                bestTransf = 'identity'
                if noSymLocus and locusPat.nonempty():
                    # if there's nontrivial symmetry, find transformation
                    # that best matches up locus/required for run with locus/required so far.
                    bestPercent = 0
                    bestTransf = None
                    for transf in stabilizedBy:
                        assert(catPat == catPat.transform(transf).shift(*shiftToMatchUp[transf]))
                        overlapPopPercent = 1
                        for i in range(2):
                            forRun = [locusForRun, savedReqForRun][i]
                            overall = [locusPat, reqPat][i]
                            matched = forRun[i].transform(transf).shift(*shiftToMatchUp[transf])
                            if (matched & overall).nonempty():
                                overlapPopPercent *= ((matched & overall).population)/((matched+overall).population)
                            else:
                                overlapPopPercent *= 0.01
                        if overlapPopPercent > bestPercent:
                            bestTransf = transf
                            bestPercent = overlapPopPercent
                # transform stuff so they match up
                transfMooreComps = [comp.transform(bestTransf).shift(*shiftToMatchUp[bestTransf])
                                for comp in mooreComps]
                reqForRun = savedReqForRun.transform(bestTransf).shift(*shiftToMatchUp[bestTransf])
                locusForRun = locusForRun.transform(bestTransf).shift(*shiftToMatchUp[bestTransf])
                antiReqForRun =  savedAntiReqForRun.transform(bestTransf).shift(*shiftToMatchUp[bestTransf])
                # analyze if components are transparent--if they are, save antirequired/locus
                # separately (maybe it's an outlier)
                # TODO: maybe we should really have a list of transpAntiReq, transpLocus, one per component
                for i in range(len(mooreComps)):
                    if (reqForRun & transfMooreComps[i]).empty():
                        compTransparentCount[i] += 1
                        transpLocus += (locusForRun & transfMooreComps[i])
                        transpAntiReq -= (transfMooreComps[i].convolve(smallZOI) - antiReqForRun)
                    elif (reqForRun & transfMooreComps[i] & reqPat).nonempty():
                        reqPat -= (transfMooreComps[i] - reqForRun)
                        antiReqPat -= (transfMooreComps[i].convolve(smallZOI) - antiReqForRun)
                        locusPat += (locusForRun & transfMooreComps[i])
                    elif len(stabilizedBy) == 1 or catPat.population > 8:
                        skippedSoups[-1].add(soup)
            else:
                if tryNum == 0:
                    sys.stderr.write("\n   skipping this soup because catalyst failed to recover:\n")
                    sys.stderr.write("   " + soup+"\n")
                continue
        
        for i in range(len(mooreComps)):
            if compTransparentCount[i] > len(catsToReactions[cat])/2:
                #if tryNum == 0:
                #    sys.stderr.write(f"\n   component was transparent for catalyst {cat}\n")
                reqPat -= mooreComps[i]
                locusPat += transpLocus
                antiReqPat -= (mooreComps[i].convolve(smallZOI) - transpAntiReq)
    
        if not noSymLocus:
            # make the required, anti-required, locus symmetrical
            for transf in stabilizedBy:
                    if transf != "identity":
                        locusPat |= locusPat.transform(transf).shift(*shiftToMatchUp[transf])
                        reqPat &= reqPat.transform(transf).shift(*shiftToMatchUp[transf])
                        antiReqPat &= antiReqPat.transform(transf).shift(*shiftToMatchUp[transf])
        
        reqTries.append(lt.pattern(''))
        reqTries[-1] |= reqPat
        antiReqTries.append(lt.pattern(''))
        antiReqTries[-1] |= antiReqPat
        locusTries.append(lt.pattern(''))
        locusTries[-1] |= locusPat
        # only do multiple trials if the first trial "looks fishy"
        if len(skippedSoups[-1]) < 3 and tryNum == 0 \
                 and all([compTransparentCount[i] < len(catsToReactions[cat])/4 or \
                    compTransparentCount[i] > 3*len(catsToReactions[cat])/4 for i in range(len(mooreComps))]):
            break
    
    bestOnes = []
    if len(locusTries) > 1:

        # for each of locus/required/antirequired,
        # find trial such that dist to 2 closest neighbors is smallest.
        for k in range(3):
            bestAvg = 0
            best = None
            diffTries = [reqTries, antiReqTries, locusTries][k]
            for i in range(5):
                similarities = sorted([Similarity(diffTries[i], diffTries[j]) 
                                                    for j in range(5) if i != j])
                if bestAvg < (similarities[-1]+similarities[-2])/2:
                    best = i
                    bestAvg = (similarities[-1]+similarities[-2])/2
            bestOnes.append(best)
        # hopefully these agree.
        if len(set(bestOnes)) == 1:
            reqPat = reqTries[mode(bestOnes)]
            antiReqPat = antiReqTries[mode(bestOnes)]
            locusPat = locusTries[mode(bestOnes)]
        else:
            sys.stderr.write("\nWARNING: best ones did not agree for catalyst " + cat + "\n")
            
            reqPat = reqTries[bestOnes[0]]
            antiReqPat = antiReqTries[bestOnes[1]]
            locusPat = locusTries[bestOnes[2]]

        soupsSkippedOften = set()
        allSoupsSkipped = skippedSoups[0] | skippedSoups[1] | skippedSoups[2] | skippedSoups[3] | skippedSoups[4] 
        sys.stderr.write(f'\n  for catalyst {cat}, we skipped these soups repeatedly (due to discrepancies in required)\n')
        numSoupsSkippedRepeatedly = 0
        for soup in allSoupsSkipped:
            if sum([int(soup in skipped) for skipped in skippedSoups]) >= 3:
                numSoupsSkippedRepeatedly += 1
                sys.stderr.write('      '+soup+"\n")
        if numSoupsSkippedRepeatedly > 3:
            sys.stderr.write("  skipped more than 3 soups repeatedly: skipping catalyst\n")
            doneWith += 1
            continue
        
    else:
        reqPat = reqTries[0]
        antiReqPat = antiReqTries[0]
        locusPat = locusTries[0]
    

    symType = ''
    if not noSymLocus:
        if len(stabilizedBy) == 1:
            symType = '*'
        elif len(stabilizedBy) == 2:
            if 'rot180' in stabilizedBy:
                symType = 'x'
            else:
                symType = '@'
        elif len(stabilizedBy) == 4:
            if 'flip_x' not in stabilizedBy:
                symType = '-'
            else:
                assert('swap_xy_flip' not in stabilizedBy)
                symType = '/'
        else:
            symType = '.'
    else:
        symType = '*'
    # generate info for csv [besides forbidden]
    bbox = locusPat.bounding_box if locusPat.nonempty() else catPat.bounding_box
    omitLocus = (catPat <= locusPat.convolve(ZOI))
    origin = ((bbox[0]+bbox[2])//2, (bbox[1]+bbox[3])//2)
    if omitLocus:
        origin = ((catbbox[0]+catbbox[2])//2, (catbbox[1]+catbbox[3])//2)
    catOffsets = GetOffset(catPat, origin)
    catInfoDict = {}
    catInfoDict['name'] = f"cat {doneWith+1}"
    catInfoDict['rle'] = TrimRLE(catPat.rle_string())
    catInfoDict['dx'] = str(catOffsets[0]); catInfoDict['dy'] = str(catOffsets[1])
    catInfoDict['symType'] = symType
    catInfoDict['absence'] = str(absence)
    partWords = ["required", "antirequired", "locus"]
    partPats = [reqPat, antiReqPat, locusPat]
    partAbbrevs = ["req", "antireq", "locus"]
    for i in range(3):
        if partPats[i].nonempty() and (partWords[i] != "locus" or not omitLocus):
            catInfoDict[partWords[i]] = TrimRLE(partPats[i].rle_string())
            offset = GetOffset(partPats[i],origin)
            catInfoDict[partAbbrevs[i]+' dx'] = str(offset[0]); catInfoDict[partAbbrevs[i]+' dy'] = str(offset[1])
            if partWords[i] not in keys:
                keys += [partWords[i], partAbbrevs[i]+' dx', partAbbrevs[i]+' dy']

    # generate forbidden
    glider = lt.pattern('bo$2bo$3o')
    forbidden = []
    forbidOffsets = []
    for orient in ['identity', 'rot90', 'rot180', 'rot270']:
        # transform pattern, instead of transforming glider. [undo transformation before getting RLE]
        # due to glide symmetry, only do rotations not reflections.
        transfPat = catPat.transform(orient)
        transfHalo = catHalo.transform(orient)
        transfReq = reqPat.transform(orient)
        transfAntiReq = antiReqPat.transform(orient)
        bbox = transfPat.bounding_box
        # glider 3x3, going down and right, so start 5 away from left or top edge
        xRanges = [[bbox[0]-5], range(bbox[0]-5,bbox[0]+bbox[2])]
        yRanges = [range(bbox[1]-5,bbox[1]+bbox[3]),[bbox[1]-5]]
        for whichRange in range(2):
            for x in xRanges[whichRange]:
                for y in yRanges[whichRange]:
                    shiftedGlider = glider.shift(x,y)
                    gen = 0
                    gliderAt = [x,y]
                    # step until right before glider collides with catalyst. Keep in mind glider
                    # might "miss" and not collide at all.
                    while shiftedGlider[1] + transfPat == (shiftedGlider+transfPat)[1] and \
                            gliderAt[0] <= bbox[0]+bbox[2]+2 and gliderAt[1] <= bbox[1]+bbox[3]+2:
                        shiftedGlider = shiftedGlider[1]
                        gen += 1
                        if gen % 4 == 0: gliderAt[0] += 1; gliderAt[1] += 1
                    if not (gliderAt[0] <= bbox[0]+bbox[2]+2 and gliderAt[1] <= bbox[1]+bbox[3]+2):
                        continue # glider missed, exited bounding box.
                    combinedStart = transfPat + shiftedGlider
                    combined = lt.pattern('')
                    combined += combinedStart
                    absentFor = 0
                    recoveredFor = 0
                    # run and check if catalyst recovers
                    while absentFor < absence and transfReq <= combined and (combined & transfAntiReq).empty():
                        combined = combined[1]
                        if transfPat <= combined and (transfHalo & combined).empty():
                            absentFor = 0
                            recoveredFor += 1
                            if recoveredFor > 4:
                                toAdd = combinedStart.transform(inverse[orient])
                                # with symmetry, possible we get the same thing twice
                                if TrimRLE(toAdd.rle_string()) not in forbidden:
                                    forbidden.append(TrimRLE(toAdd.rle_string()))
                                    forbidOffsets.append(GetOffset(toAdd,origin))
                                break
                        else:
                            recoveredFor = 0
                            absentFor += 1
    for i in range(len(forbidden)):
        catInfoDict[f"forbidden {i+1}"] = forbidden[i]
        catInfoDict[f"forbid {i+1} dx"] = forbidOffsets[i][0]
        catInfoDict[f"forbid {i+1} dy"] = forbidOffsets[i][1]
        if f"forbidden {i+1}" not in keys:
            keys += [f"forbidden {i+1}", f"forbid {i+1} dx",f"forbid {i+1} dy"]
    doneWith += 1
    sys.stderr.write(f"\rdone with catalyst {doneWith} of {len(catsToReactions)}")
    rows.append(catInfoDict)

with open(sys.argv[2], "w") as f:
    writer = csv.DictWriter(f, fieldnames = keys)
    writer.writeheader()
    writer.writerows(rows)