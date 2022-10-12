# arguments: sampleSoupFile outputFile
# input file formatting assumed to be
# [catalyst rle]:
# \t[soup1]
# \t[soup2]
import sys
from math import floor,ceil
import lifelib
import csv
sess = lifelib.load_rules("b3s23")
lt = sess.lifetree(n_layers=1)
allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]
#inverses = ["identity", "rot90", "rot180", "rot270" , "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]
inverse = {}
for orient in allOrientations:
    if orient in ["rot270", "rot90"]:
        otherInd = (["rot270", "rot90"].index(orient) + 1) % 2
        inverse[orient] = ["rot270", "rot90"][otherInd]
    else:
        inverse[orient] = orient

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
    shouldBeContained = lt.pattern('$'.join(64*['64o'])).shift(-32,-32)
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
                
    if not (shouldBeContained <= resultPat):
        sys.stderr.write("problematic soup: \n")
        sys.stderr.write(soup+"\n")
        sys.stderr.write("problematic catalyst: \n")
        sys.stderr.write(cat+"\n")
        sys.stderr.write("problematic result pat: \n")
        sys.stderr.write(resultPat.rle_string()+"\n")
        sys.stderr.write("problematic oriented patterns: \n")
        sys.stderr.write(resultPat.rle_string()+"\n")
    assert(shouldBeContained <= resultPat)
    if not( bestOrientation is None):
        return closest, allOrientations[bestOrientation]
    sys.stderr.write("troubles with catalyst " + cat + "\n")
    sys.stderr.write("      and soup " + soup + "\n")
    return None, None

def GetOffset(otherPat, origin):
    if otherPat.empty():
        return None
    upperRightCorner = otherPat.bounding_box[:2]
    return (upperRightCorner[0]-origin[0], upperRightCorner[1]-origin[1])

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

keys = ['name', 'absence','rle','dx','dy', 'symType']
rows = []
doneWith  = 0
for cat in catsToReactions:
    catPat = lt.pattern(cat)
    assert(catPat.nonempty())
    catbbox = catPat.bounding_box
    assert(catbbox[3] < 64 and catbbox[2] < 64)

    reqPat = lt.pattern('')
    reqPat |= catPat
    
    smallZOI = lt.pattern('')
    smallZOI[-1:2, -1:2] = 1
    cross = lt.pattern('')
    cross[0, -1:2] = 1
    cross[ -1:2,0] = 1
    ZOI = cross.convolve(smallZOI)
    absence = 5
    
    catHalo = catPat.convolve(smallZOI) - catPat
    antiReqPat = catPat.convolve(smallZOI)
    locusPat = lt.pattern('')
    catbbox = catPat.bounding_box
    bigBox = lt.pattern('')
    bigBox[(catbbox[0]-64):(catbbox[0]+catbbox[2]+64),
           (catbbox[1]-64):(catbbox[1]+catbbox[3]+64) ] = 1
    comps = (bigBox - antiReqPat).components()
    for comp in comps:
        if comp.population < 8:
            antiReqPat += comp
    antiReqPat -= catPat
    startAntiReqPat = lt.pattern('')
    startAntiReqPat += antiReqPat
    transparentCount = 0

    symType = ''
    stabilizedBy = set()
    for transf in allOrientations:
        if (catPat.transform(transf).match(catPat)).nonempty():
            stabilizedBy.add(transf)
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
    
    for soup in catsToReactions[cat]:
        soupPat = lt.pattern(soup)
        soupBbox = soupPat.bounding_box
        location, orient = LocateCatalyst(catPat, soupPat)
        centeredCatSoup = TilePat(soupPat.shift(-1*location[0],
                -1*location[1]).transform(inverse[orient]),
                [-4*64,-4*64,5*64,5*64])
        if (centeredCatSoup & catHalo).nonempty():
            sys.stderr.write(catHalo.rle_string())
            sys.stderr.write(centeredCatSoup.rle_string())
        assert(catPat <= centeredCatSoup)
        assert((centeredCatSoup & catHalo).empty())
        soupCenteredNoCat = centeredCatSoup - catPat
        activated = False
        activatedGen = -1
        recoveredFor = -1
        reqForRun = lt.pattern('')
        reqForRun |= catPat
        antiReqForRun = lt.pattern('')
        antiReqForRun |= startAntiReqPat
        for n in range(200):
            centeredCatSoup = centeredCatSoup[1]
            if ((centeredCatSoup & catHalo).nonempty()) or not (catPat <=centeredCatSoup):
                if not activated:
                    locusPat += (centeredCatSoup & catHalo).convolve(smallZOI) & catPat
                    activated = True
                    activatedGen = n
                antiReqForRun -= centeredCatSoup
                reqForRun &= centeredCatSoup
                recoveredFor = 0
            elif activated:
                recoveredFor += 1
                if recoveredFor > 4:
                    absence = max(absence, n-activatedGen-3)
                    break
            elif not activated: # the case where the catalyst prevents adjacent births, but doesn't
                # have cells born next to it.
                soupCenteredNoCat = soupCenteredNoCat[1]
                if (soupCenteredNoCat & catHalo).nonempty():
                    activated = True
                    activatedGen = n
                    locusPat += (soupCenteredNoCat & catHalo).convolve(smallZOI) & catPat
        if recoveredFor > 0:
            if reqForRun.nonempty() and (reqForRun & reqPat).nonempty():
                reqPat &= reqForRun
                antiReqPat &= antiReqForRun
            elif reqForRun.empty():
                transparentCount += 1
            elif len(stabilizedBy) == 1 or catPat.population > 8:
                sys.stderr.write("\n   skipping this soup because it'd make the required empty:\n")
                sys.stderr.write("   " + soup+"\n")
        else:
            sys.stderr.write("\n   skipping this soup because catalyst failed to recover:\n")
            sys.stderr.write("   " + soup+"\n")
            continue
    if transparentCount > len(catsToReactions[cat])/2:
        reqPat = lt.pattern('')
        antiReqPat = lt.pattern('')
    
    # make the required, anti-required, locus symmetrical if the catalyst has symmetry
    for transf in stabilizedBy:
        if transf != "identity":
            shiftBy = (catPat.transform(transf).match(catPat)).firstcell
            locusPat |= locusPat.transform(transf).shift(-shiftBy[0], -shiftBy[1])
            reqPat &= reqPat.transform(transf).shift(-shiftBy[0], -shiftBy[1])
            antiReqPat &= antiReqPat.transform(transf).shift(-shiftBy[0], -shiftBy[1])

    bbox = locusPat.bounding_box
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
    #catData.append(f"cat {TrimRLE(catPat.rle_string())} {absence} {catOffsets[0]} {catOffsets[1]} {symType}")
    partWords = ["required", "antirequired", "locus"]
    partPats = [reqPat, antiReqPat, locusPat]
    partAbbrevs = ["req", "antireq", "locus"]
    #partOffsets = [GetOffset(aPat,origin) for aPat in partPats]
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
                    while absentFor < absence and transfReq <= combined:
                                # anti-required tends to be aggressive, not allow for boat bits,
                                # so leave out for now.
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