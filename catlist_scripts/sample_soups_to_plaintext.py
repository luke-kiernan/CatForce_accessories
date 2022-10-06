# arguments: sampleSoupFile outputFile
# input file formatting assumed to be
# [catalyst rle]:
# \t[soup1]
# \t[soup2]
import sys
from math import floor,ceil
import lifelib
sess = lifelib.load_rules("b3s23")
lt = sess.lifetree(n_layers=1)
allOrientations = ["identity", "rot270", "rot180", "rot90", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]
inverses = ["identity", "rot90", "rot180", "rot270" , "flip_x", "flip_y", "swap_xy", "swap_xy_flip"]

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
        transfPat = pat.transform(op)
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

    soupBbox = reactionPat.bounding_box
    middleOfSoup = ((soupBbox[0]+soupBbox[2])//2,(soupBbox[1]+soupBbox[3])//2 )
    resultPat = TilePat(reactionPat)
    bestOrientation = None
    dist = 1000
    # for large catalysts, might be multiple matches due to wrap-around
    # and TilePat, so we prefer those close to the middle of the tile.
    closest = [None, None]
    for i in range(8):
        matches = resultPat.match(orientedPats[i], dead = orientedDeadCells[i])
        if(matches.nonempty()):
            for cell in matches.coords():
                if max(abs(cell[0]-middleOfSoup[0]),abs(cell[1]-middleOfSoup[1])) < dist:
                    closest = cell[:]
                    dist = max(abs(cell[0]-middleOfSoup[0]),abs(cell[1]-middleOfSoup[1]))
                    bestOrientation = i
                    shouldBeContained = orientedPats[i].shift(*closest)
                    if dist <  max(abs(middleOfSoup[0]),abs(middleOfSoup[1]))\
                             and shouldBeContained <= resultPat:
                        break
        if dist < max(abs(middleOfSoup[0]),abs(middleOfSoup[1])):
            break
    assert(shouldBeContained <= resultPat)
    if not( bestOrientation is None):
        return closest, bestOrientation
    sys.stderr.write("troubles with catalyst " + cat + "\n")
    sys.stderr.write("      and soup " + soup + "\n")
    return None, None

def GetOffset(otherPat, origin):
    if otherPat.empty():
        return None
    upperRightCorner = otherPat.bounding_box[:2]
    return (upperRightCorner[0]-origin[0], upperRightCorner[1]-origin[1])
x = 0
y = 0
catsToReactions = {}
catData = []
currentCat = None
with open(sys.argv[1], "r") as f:
    for line in f:
        if ':' in line and 'o' in line.split(':')[0]:
            catsToReactions[line.split(':')[0]] = []
            currentCat = line.split(':')[0]
        elif 'o' in line.rstrip().lstrip() and not (currentCat is None):
            catsToReactions[currentCat].append(line.rstrip().lstrip())
catEntries = []

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
    bigBox[(catbbox[0]-12):(catbbox[0]+catbbox[2]+24),
           (catbbox[1]-12):(catbbox[1]+catbbox[3]+24) ] = 1
    comps = (bigBox - antiReqPat).components()
    for comp in comps:
        if comp.population < 8:
            antiReqPat += comp
    antiReqPat -= catPat
    for soup in catsToReactions[cat]:
        soupPat = lt.pattern(soup)
        soupBbox = soupPat.bounding_box
        middleOfSoup = ((soupBbox[0]+soupBbox[2])//2,(soupBbox[1]+soupBbox[3])//2 )
        location, orientNum = LocateCatalyst(catPat, soupPat)
        centeredCatSoup = TilePat(soupPat.shift(-1*location[0],
                -1*location[1]).transform(inverses[orientNum]),
                [-4*64,-4*64,5*64,5*64])
        if (centeredCatSoup & catHalo).nonempty():
            sys.stderr.write(catHalo.rle_string())
            sys.stderr.write(centeredCatSoup.rle_string())
        assert(catPat <= centeredCatSoup)
        assert((centeredCatSoup & catHalo).empty())
        activated = False
        activatedGen = -1
        recoveredFor = -1
        for n in range(100):
            centeredCatSoup = centeredCatSoup[1]
            if ((centeredCatSoup & catHalo).nonempty()) or not (catPat <=centeredCatSoup):
                if not activated:
                    locusPat += (centeredCatSoup & catHalo).convolve(smallZOI) & catPat
                    activated = True
                    activatedGen = n
                antiReqPat -= centeredCatSoup
                reqPat &= centeredCatSoup
                recoveredFor = 0
            elif activated:
                recoveredFor += 1
                if recoveredFor > 4:
                    absence = max(absence, n-activatedGen)
                    break
    # work out symmetry
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
    # TODO:
    # make the required, anti-required, locus symmetrical if the catalyst has symmetry
    #for transf in stabilizedBy:
    #    if transf != "identity":
    #        shiftBy = (catPat.transform(transf).match(catPat)).firstcell
    bbox = locusPat.bounding_box
    omitLocus = (catPat <= locusPat.convolve(ZOI))
    origin = ((bbox[0]+bbox[2])//2, (bbox[1]+bbox[3])//2)
    if omitLocus:
        origin = ((catbbox[0]+catbbox[2])//2, (catbbox[1]+catbbox[3])//2)
    catOffsets = GetOffset(catPat, origin)
    catData.append(f"cat {TrimRLE(catPat.rle_string())} {absence} {catOffsets[0]} {catOffsets[1]} {symType}")
    partWords = ["required", "antirequired", "locus"]
    partPats = [reqPat, antiReqPat, locusPat]
    partOffsets = [GetOffset(aPat,origin) for aPat in partPats]
    for i in range(3):
        if partPats[i].nonempty() and (partWords[i] != "locus" or not omitLocus):
            catData[-1] += f" {partWords[i]} {TrimRLE(partPats[i].rle_string())} {partOffsets[i][0]} {partOffsets[i][1]}"
    sys.stderr.write("done with catalyst " + TrimRLE(catPat.rle_string()) + "\n")
with open(sys.argv[2], "w") as f:
    for line in catData:
        f.write(line+"\n")   
