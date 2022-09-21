#!/bin/bash

# format: [name]&[rle]
# listed in order of https://conwaylife.com/wiki/List_of_common_evolutionary_sequences
# to try to catch as many thing as possible with fitlers, some of these are descendants
# uses forms recommended by HotDogPi [thanks hdp!]
HF=$"Honeyfarm&2b3o\$bo3bo\$o5bo\$o5bo\$o5bo\$bo3bo\$2b3o!"
R=$"R&b2o\$bobo\$o2bo\$b2o\$bo!"
Pi=$"Pi&b3o\$o3bo\$2ob2o!"
S=$"Stairstep&3o\$o2b2o\$2o2bo\$2b3o!"
Century=$"Century&b3o\$3bo\$ob2o\$2o!"
B=$"B&2bo\$b3o\$2o2bo\$2b3o!"
W=$"Wing&b2o\$o2bo\$bobo\$2b2o!"
TL=$"TrafficLight&3o\$obo\$3o!"
TD=$"Teardrop&b3o\$o2bo\$o2bo\$b2o!"
I=$"I&3b2o\$bo2bo\$o3bo\$bobo\$2bo!"
O=$"Octomino&b3o\$o\$obo\$3o!"
J33=$"Jason33&2bo\$2obo\$o3bo\$b3o!"
Pro=$"Procrastinator&o\$3o\$3bo\$b3o!"
BT=$"BlonkTie&bobo\$o3bo\$o4bo\$bo3bo\$2b3o!"
IW=$"Iowona&2b2o\$o3bo\$o2bo\$b2o!"
U=$"UTurner&2bo\$bobo\$o\$o2bo\$2o!"
HDesc=$"HDescendant&2o\$o2b2o\$bo2bo\$b3o!"
RT=$"RTurner&2b2o\$4bo\$o2bo\$3o!"

# block and fishhook. 4 catalyst searching with this is  feasible.
read -r -d '' topTwo << "EOM"
cat 2o$2o 10 -1 -1 + forbidden obo$b2o$bo2$2o$2o! 0 -4 forbidden 2o3bo$2ob2o$4b2o! 0 0 forbidden b2o$b2o2$bo$2o$obo! -1 0 forbidden 2o$b2ob2o$o3b2o! -4 -1 forbidden o3b2o$b2ob2o$2o! -4 0 forbidden obo$2o$bo2$b2o$b2o! -1 -4 forbidden 4b2o$2ob2o$2o3bo! 0 -1 forbidden 2o$2o2$bo$b2o$obo! 0 0 required o! -1 -1
cat 2o$o$b3o$3bo! 10 -2 -2 * forbidden 2o$o$b3ob2o$3bobobo$6bo! -2 -2 forbidden bo$obo$2o2$2o$o$b3o$3bo! -2 -6 forbidden bo$2bo$3o2$3b2o$3bo$4b3o$6bo! -5 -6 forbidden 2bo$obo$b2o2$3b2o$3bo$4b3o$6bo! -5 -6 required bo$o$b3o$3bo! -2 -2 locus 2o$o! -2 -2
EOM
topTwo+=$'\n'

# block, fishhook, tub, boat, and eater2. 
# required cells and symmetries optimized for high-volume searching. Long tub instead of boat due to Kazyan, saves 2 catalysts.
read -r -d '' Kazyan << "EOM"
cat 2o$2o 10 -1 -1 / forbidden obo$b2o$bo2$2o$2o! 0 -4 forbidden 2o3bo$2ob2o$4b2o! 0 0 forbidden b2o$b2o2$bo$2o$obo! -1 0 forbidden 2o$b2ob2o$o3b2o! -4 -1 forbidden o3b2o$b2ob2o$2o! -4 0 forbidden obo$2o$bo2$b2o$b2o! -1 -4 forbidden 4b2o$2ob2o$2o3bo! 0 -1 forbidden 2o$2o2$bo$b2o$obo! 0 0 required o! -1 -1
cat 2o$o$b3o$3bo! 10 -2 -2 * forbidden 2o$o$b3ob2o$3bobobo$6bo! -2 -2 forbidden bo$obo$2o2$2o$o$b3o$3bo! -2 -6 forbidden bo$2bo$3o2$3b2o$3bo$4b3o$6bo! -5 -6 forbidden 2bo$obo$b2o2$3b2o$3bo$4b3o$6bo! -5 -6 required bo$o$b3o$3bo! -2 -2 locus 2o$o! -2 -2
cat bo$obo$bobo$2bo! 6 -1 -1 | required obo$bobo! -1 0
cat bo$obo$bo! 10 -1 -1 | required o! -1 0
cat 2obo$2ob3o$6bo$2ob2obo$bob2ob2o$bo$2b3ob2o$4bob2o! 10 -3 -3 | forbidden 2bo$obo$b2o2$3b2obo$3b2ob3o$9bo$3b2ob2obo$4bob2ob2o$4bo$5b3ob2o$7bob2o! -6 -7 forbidden bo$2bo$3o2$3b2obo$3b2ob3o$9bo$3b2ob2obo$4bob2ob2o$4bo$5b3ob2o$7bob2o! -6 -7 forbidden 2bo$obo$b2o$4b2obo$4b2ob3o$10bo$4b2ob2obo$5bob2ob2o$5bo$6b3ob2o$8bob2o! -7 -6 forbidden bo$2bo$3o$4b2obo$4b2ob3o$10bo$4b2ob2obo$5bob2ob2o$5bo$6b3ob2o$8bob2o! -7 -6 forbidden 2obo$2ob3o$6bo$2ob2obo$bob2ob2o$bo$2b3ob2o$4bob2o$9b3o$9bo$10bo! -3 -3 forbidden 2obo$2ob3o$6bo$2ob2obo$bob2ob2o$bo$2b3ob2o$4bob2o$9b2o$9bobo$9bo! -3 -3 forbidden 2obo$2ob3o$6bo$2ob2obo$bob2ob2o$bo$2b3ob2o$4bob2o2$8b3o$8bo$9bo! -3 -3 forbidden 2obo$2ob3o$6bo$2ob2obo$bob2ob2o$bo$2b3ob2o$4bob2o2$8b2o$8bobo$8bo! -3 -3 required 3bo$3b3o$6bo$2ob2obo$bob2ob2o$bo$2b3o$4bo! -3 -3 locus 2o$2o5$6b2o$6b2o! -3 -3
EOM
Kazyan+="\n"


# begin main parameters_to_edit

# only used if CatForce.cpp and LifeAPI.h aren't present, and to set tinyList.
CATFORCE_DIRECTORY="../../CatForceClean"

# read in from text file
tinyList=$(<$CATFORCE_DIRECTORY/catlists/The_Tiny_List_of_Catalysts)

# CatForce input file parameters
whichCatalysts="$Kazyan"
maxGen=200
lastGen=50
numCats=2
earlyParams="start-gen 1\nmax-gen $maxGen\nlast-gen $lastGen\nstable-interval 18\nnum-catalyst $numCats\nnum-transparent 1\n"
oscParams=$'stop-after-cats-destroyed 10\n'

# script parameters: which offests, regions, symmetries to iterate through.
searchName="2Kazyan" # file naming purposes. the temp file CatForce reads from is 
# named with this suffix, so that you can have multiple copies of this script running at once.
stopAt=15 # search covers: active region has center 0 <= x <= stopAt, 0<= y <= stopAt
activeRegions=("$HF" "$R" "$Pi" "$S" "$Century" "$B" "$W" "$TL" "$TD"  "$I" "$O" "$J33" "$Pro" "$BT" "$IW" "$U" "$HDesc" "$RT" )
symOptions='D2| D2|even D2/ C2 C2|even C2-even C2even C4 C4even  D4+ D4+|even D4+-even D4+even D4x D4xeven D8 D8even'
# try to catch flippers. does not take into account whether flipped catalyst placements are viable
# (but that's why we have post-search filtering scripts)
flippers="true"

debug="false" # print CatForce progress messages and save full-report files. deletes them if they're empty, though.

# end main parameters to edit

# notes/clarifications:

# if symmetry isn't a D2 variant, via composing with reflections across the axes, we can assume the center of the active region is below and
# to the right of the center of symmetry. Via further composing with a reflection across y=x, we can assume x >= y.
# however, that last part swaps | with -. so if you include D4+|even, for example, you should probably also include D4+-even.
# for D2, it suffices to have D2- and D2/.

# there will be some overlap: cases where the center of the active region lies on the diagonal x = y 
# or along the horizontal line thru the center of symmetry are covered twice, but oh well.
# Working in terms of the vector between the 2+ active regions, instead of the center of the active 
# region, would probably be cleaner, but rewriting this isn't high on the priority list.

# first 8 active regions and C2 variants seem to be the most productive.
# As far as I know, something has been found via symmetric CatForce search in all symmetries all besides the D8 variants.

if ! [[ -f 'CatForce.cpp' ]]; then cp $CATFORCE_DIRECTORY/CatForce.cpp .; fi
if ! [[ -f 'LifeAPI.h' ]]; then cp $CATFORCE_DIRECTORY/LifeAPI.h .; fi

if ! [[ -f 'active_region_occurs' ]] || ! [[ -f 'CatForce' ]] || [[ ! -f 'compute_flipper_filters' && $flippers == "true" ]]
then
    echo "executables not found...trying running make"
    if ! make
    then
        echo "make did not run successfully. exiting script"
        exit 0;
    fi
fi

if ! [[ -d masks ]]
then
    mkdir masks
fi

for obj in "${activeRegions[@]}"
do
    rle=${obj#*&} # regex for drop everything before the first &
    name=${obj%&*} # regex for drop everything after the last &

    echo "Starting active object ${name}"

    if [ ! -d "$name" ]
    then
        mkdir "$name"
    fi
    rotatedObj=$( python3 generateOrientedRLEs.py "$rle" )
    IFS=$'\n' read -d '' -ra rotatedRLEs < <(printf '%s\0' "$rotatedObj") # found via stack overflow, split up by \n.
    # rotatedRLEs = all orientations of the object [as a set: ie it detects redundancy, only 1 for honeyfarm, 4 for pi, etc]

    for sym in $symOptions
    do
        if [ $sym = 'D2backslash' ]
        then
            echo 'WARNING: D2\ no longer supported. It causes placement headaches, due to wanting to put active region in 2nd quadrant not 1st.'
            echo '         Searching D2/ intead.'
            sym="D2/"
        fi

        echo "    Starting symmetry ${sym}"
        soFar=0

        symParams="symmetry ${sym}\n"

        if  [ "${sym:2:4}" = "even" ] || [ $sym = 'D4+even' ] || [ $sym = 'D4xeven' ]
        then
            centerOfSymX="-0.5"
            centerOfSymY="-0.5"
        elif [ $sym = 'C2|even' ] || [ $sym = 'D4+|even' ]
        then
            centerOfSymX="-0.5"
            centerOfSymY="0"
        elif [ $sym = 'C2-even' ] || [ $sym = 'D4+-even' ]
        then
            centerOfSymX="0"
            centerOfSymY="-0.5"
        else
            centerOfSymX="0"
            centerOfSymY="0"
        fi

        startTimeForSymAndObj=$SECONDS

        for orientNum in ${!rotatedRLEs[@]}
        do
            orientedRLE="${rotatedRLEs[$orientNum]}"

            dims=$( python3 getWidthHeight.py $orientedRLE )

            height=${dims#*,}
            width=${dims%,*}
            
            if [ "${sym:0:2}" != "D2" ]
            then
                # naive total,  (area of 0 <= y <= x <= stopAt) * (num orientations)
                approxTotal=$(( ($stopAt+1) * ($stopAt+1) * ${#rotatedRLEs[@]} / 2 ))
                # rough correction for overlap
                if [[ "${sym:0:1}" == "C" ]]
                then
                    approxTotal=$(($approxTotal - $width / 2 * ($height / 2) ))
                elif [ "${sym:0:3}" == "D4+" ] ||  [ "${sym:0:2}" == "D8" ]
                then
                    approxTotal=$(($approxTotal - $height / 2 * $stopAt))
                fi
                if [ "${sym:0:3}" == "D4x" ] || [ "${sym:0:2}" == "D8" ]
                then
                    greaterDim=$( [[ $width > $height ]] && echo $width || echo $height)
                    approxTotal=$(($approxTotal - $greaterDim / 2 * $stopAt))
                fi
            else
                approxTotal=$(( ($stopAt+1) * ${#rotatedRLEs[@]} ))
            fi

            # because we count points not edges, something with width 2, height 3, corner (0,0)
            # has center (0.5, 1), not (1, 1.5). I want x0 + width/2 to be the center.
            # So I'll be non-standard and count edges, not points.
            height=$(( $height - 1 ))
            width=$(( $width - 1 ))

            # due to the above, odd width things have center at half-integer.
            if [ $(( $width % 2 )) == "1" ]
            then
                startCenterX="-0.5"
            else
                startCenterX="0"
            fi

            if [ $(( $height % 2 )) == "1" ]
            then
                startCenterY="-0.5"
            else
                startCenterY="0"
            fi


            for centerX in $(seq $startCenterX 1 $stopAt)
            do
                for centerY in $(seq $startCenterY 1 $stopAt)
                do
                    # upper left corner.
                    x0=$( echo "$centerX - $width / 2" | bc -l )
                    y0=$( echo "$centerY - $height / 2" | bc -l )
                    x0=$( printf "%.0f" $x0 )
                    y0=$( printf "%.0f" $y0 )

                    activeRegionOccurs=$(./active_region_occurs ${orientedRLE} $x0 $y0 $sym)

                    if [[ "${sym:0:2}" != "D2" ]] && [[ $(echo "$centerX >= $centerOfSymX" |bc -l) = "1" ]] && [[ $(echo "$centerY >= $centerOfSymY" |bc -l) = "1" ]] && [[ $(echo "$centerX >= $centerY" |bc -l) = "1" ]]
                    then
                        inCorrectRegion="true"
                    elif [[ $sym = "D2|" || $sym = "D2|even" ]] && [[ $(echo "$centerY == $startCenterY" |bc -l) = "1" ]]
                    then
                        inCorrectRegion="true"
                    elif [[ $sym = "D2-" || $sym = "D2-even" ]] && [[ $(echo "$centerX == $startCenterX" |bc -l) = "1" ]]
                    then
                        inCorrectRegion="true"
                    elif [[ $sym = 'D2/' ]] && [[ $(echo "$centerY == $centerX" |bc -l) = "1" || $(echo "$centerX - $centerY == 0.5" |bc -l) = "1"  ]]
                    then
                        inCorrectRegion="true"
                    else
                        inCorrectRegion="false"
                    fi

                    if [[ "$inCorrectRegion" == "true" && "$activeRegionOccurs" == "true" ]]
                    then

                        # mvr's CatForce automatically intersects with a fundamental domain.
                        searchArea="search-area -30 -30 60 60\n"
                        
                        #pattern
                        patParams="pat ${orientedRLE} $x0 ${y0}\n"

                        # output file naming (avoid |, / in file names)
                        symLen=${#sym}
                        if [[ $sym = 'D2|' || $sym = 'D2|even' ]]
                        then
                            symAsString="${sym/\|/ReflectAcrossYAxis}"
                        elif [[ symLen -ge 5  && "${sym:symLen-5:5}" = '|even' ]] 
                        then
                            symAsString="${sym/\|/horizontal}"
                        elif [[ $sym = 'D2/' ]]
                        then
                            symAsString='D2slash'
                        else
                            symAsString="$sym"
                        fi
                    
                        # filters
                        if [ "$flippers" = "true" ]
                        then
                            reoccurFilter=$(./compute_flipper_filters $orientedRLE $x0 $y0 $sym 5 $maxGen)
                            reoccurFilter+="\n"
                        else
                            reoccurFilter="filter 5-$maxGen $orientedRLE $x0 ${y0}\n"
                        fi
                    
                        outfileName="${name}_${symAsString}_orient${orientNum}_${centerX/./,}_${centerY/./,}_${searchName}.rle"
                        outfileFullName="${name}_${symAsString}_orient${orientNum}_${centerX/./,}_${centerY/./,}_${searchName}_full.rle"

                        if [ "$debug" == "true" ]
                        then
                            outfileParams="output ${outfileName}\nfull-report ${outfileFullName}\n" # debugging purposes, save full report
                        else
                            outfileParams="output ${outfileName}\n"
                        fi

                        # create search file
                        echo -e "${earlyParams}${symParams}${searchArea}${whichCatalysts}${patParams}${reoccurFilter}${outfileParams}${oscParams}" > temp_${searchName}.txt

                        # run the search (don't display output if debug is false)
                        if [ "$debug" == "false" ]
                        then
                            ./CatForce temp_${searchName}.txt > /dev/null
                        else
                            ./CatForce temp_${searchName}.txt
                        fi

                        # delete RLEs with no results
                        if [ "$(wc -c ${outfileName} | awk '{print $1}')" = "28" ]
                        then
                            rm $outfileName
                        else
                            mv $outfileName "${name}/${newOutfileName}"
                            printf "\n        Possible results: see ${outfileName}\n"
                        fi

                        # debug mode and no full results => delete full results.
                        if [ -f $outfileFullName ] && [ "$(wc -c ${outfileFullName} | awk '{print $1}')" = "28" ]
                        then
                            rm $outfileFullName
                        fi

                        # update count.
                        soFar=$(($soFar + 1))
                        timeForSymAndObj=$(echo "$SECONDS - $startTimeForSymAndObj" |bc -l)

                        # display progress message.
                        toPrint="\r        ${soFar} done of roughly $approxTotal for $name with symmetry $sym. "
                        if [[ $(echo "$timeForSymAndObj > 3600 * $soFar" |bc -l) = "1" ]]
                        then
                            toPrint+="Average $(($timeForSymAndObj  / ($soFar * 3600) )) h $((($timeForSymAndObj / (60 * $soFar) ) % 60)) min $(( ($timeForSymAndObj / $soFar ) % 60)) sec per search"
                        elif [[ $(echo "$timeForSymAndObj > 60 * $soFar" |bc -l) = "1" ]]
                        then
                            toPrint+="Average $((($timeForSymAndObj / (60 * $soFar) ) % 60)) min $(( ($timeForSymAndObj / $soFar ) % 60)) sec per search"
                        else
                            toPrint+="Average $(( ($timeForSymAndObj / $soFar ) % 60)) sec per search"
                        fi
                        printf "$toPrint"
                    fi
                done
            done # locations
        done # orientations
        printf " DONE, after $(($timeForSymAndObj  /  3600)) h $((($timeForSymAndObj / 60 ) % 60)) min $(($timeForSymAndObj % 60)) sec.\n"
    done # symmetries
done # objects

