import sys
import re
def RLE_to_array(rleString):
    rleString = re.sub(r"\s+", "", ''.join(rleString))
    if rleString[-1] == '!':
        rleString = rleString.rstrip('!')
    cellString = ''
    def decodeRLEchar(char):
        if char == 'b': return '0'
        if char == 'o': return '1'
        if char == '$': return '\n'
        print(f"unrecognized character {repr(char)}")
        assert(False)
    
    i=0
    # deal with numbers: need to turn 13b into [13 zeros]
    while(i < len(rleString)):
        while(i < len(rleString) and not rleString[i].isnumeric()):
            cellString += decodeRLEchar(rleString[i])
            i += 1
        if i < len(rleString):
            numString = '' # number could have multiple digits.
            while(rleString[i].isnumeric()):
                numString += rleString[i]
                i += 1
            cellString += int(numString)*decodeRLEchar(rleString[i])
            i += 1
    genZero = [list(row) for row in cellString.split('\n')]

    width = max([len(row) for row in genZero]) # pad with zeros.
    for row in genZero:
        while len(row) < width:
            row.append(0)
    
    return genZero

def matrixToRLE(matrix):
    def encodeRLEchar(num):
        if num == 1: return 'o'
        if num == 0: return  'b'
    rleString = ''
    for row in matrix:
        rleString += ''.join(row).rstrip('0')
        rleString += '\n'
    rleString = rleString.replace('0', 'b')
    rleString = rleString.replace('1', 'o')
    rleString = rleString.replace('\n', '$')
    rleString = rleString.rstrip('$').lstrip('$')
    condensedRLEstring = ''
    i = 0
    while i < len(rleString):
        runOf = rleString[i]
        runLength = 1
        i += 1
        while ( i < len(rleString) and rleString[i] == runOf):
            i += 1
            runLength += 1
        if runLength > 1:
            condensedRLEstring += f'{runLength}{runOf}'
        else:
            condensedRLEstring += runOf
    condensedRLEstring += '!'
    return condensedRLEstring



if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Wrong length of argv??? Here\'s what we got, 1 item per line:')
        for item in sys.argv: print(repr(item))
    rleMatrix = RLE_to_array(sys.argv[1])
    print(f"{len(rleMatrix[0])},{len(rleMatrix)}")

