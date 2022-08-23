# CatForce_accessories
Helper scripts for CatForce searching and some additional catalyst lists.

**Catlist Scripts**
These are scripts for generating CatForce input files. There are three main file formats: rle, csv, and plaintext. 
- csv format: the data for each catalyst is stored as entries in a table, 1 row per catalyst. Good for editing or splicing together lists.
- rle format: one catalyst per line. The leftmost entry should be the catalyst, with required cells as states 3, 4, or 5. The rest then are the forbidden patterns. Good for adding things to lists.
- plaintext format: the format accepted by CatForce.
Four scripts are included:
- `csv_to_plaintext.py` example usage `python3 csv_to_plaintext.py input.csv [name] > output.txt`. If the optional `name` parameter is provided, then the catalysts are labeled with their name (from the csv) and period. Works for peridoic and stable. 
- `plaintext_to_csv.py` example usage `python3 plaintext_to_csv.py input.txt output.csv`. Creates a column for name, but that must be entered by hand. Works for peridoic and stable. 
- `csv_to_rle.py` Golly script. Columns are assumed to be: name, rle, absence, dx, dy, symType, required, req dx, req dy, forbidden 1, forbid 1 dx, forbid 1 dy, forbidden 2, etc. in that order. Works for peridoic and stable. 
- `rle_to_csv.py` Golly script. Columns for name, absence, period, and symmetry type are created, but those fields must be entered by hand. Works for peridoic and stable: just delete the period column for stable. Make sure to include at least 5 cells of empty space between rows, and between patterns within each row.

**Additional Catalyst Lists**
These are some additional lists of catalysts or sparkers that I've made. For each list, all three file formats are included: add catalysts to the rle or alter parameters in the csv, then re-run the scripts to generate new lists. Included are:
- `high_period_sparkers.txt`: sparkers of period 6-10, with a few of period 16 too, for use with periodic CatForce. Not intended to be comprehensive: contains one of each type, for the most part. Sourced from Lifewiki.
- in the `by-period` folder, there's also a compilation of sparkers of period 3 and period 4. Periods 6-10 and 16 are included in that folder well, but `high_period_sparkers.txt` is just the combined contents of those lists. The period 4 ones are sourced from [this forum post](https://conwaylife.com/forums/viewtopic.php?f=2&t=1437&p=85170&hilit=sparker+collection+p4#p85170) and [Lifewiki](https://conwaylife.com/wiki/Category:Sparkers_with_period_4); the period 3 ones are from Lifewiki.
- `mononers_and_variants`: a list of stable catalysts. Intended for completing partials, where space is tight and the usual dimers may not fit. Again, not comprehensive: this list omits several large catalysts on `The_Big_List_of_Catalysts`, such as the hivepush catalysts. Most of these catalysts come from [Tutorials/Catalyses](https://conwaylife.com/wiki/Tutorials/Catalyses) and [Catalyst Discussion Thread](https://conwaylife.com/forums/viewtopic.php?f=2&t=1878&p=24656&hilit=catalyst+testing#p24440).

**Bash Search Script**
Will add soon.

**Post Search Filtering**
When searching for oscillators with the bash script, often one mechanism will show up at different offsets, creating a boatload of uninteresting result files that are tedious to go through. Both of these scripts require [Lifelib](https://gitlab.com/apgoucher/lifelib/-/tree/master), and only work with stable catalysts.

`categorize_by_match_multi.py`: combine results between multiple CatForce results files, based off in which generation
a specific subpattern occurs.

Example usage `python3 categorize_by_match_multi.py [rleToMatch] [start]-[end] [matchType] results1.rle...resultsN.rle`.

- rleToMatch: an rle string, like `3o$2bo$3o!`. (Rmemeber to escape the dollar signs or enclose in single quotes.)
- start-end: which generations to look for the pattern.
- matchType: options are `same` (look for it in the same orientation and location as generation 0), `([x],[y])` (look for it in a fixed spot, same orientation as provided). If omitted, the first match in any location or orientation is used.
- RLE files: all the categories from these results files will be combined into one RLE, on the basis
of in which generation a match is found. (If no match is found, the pattern is omitted.)

`categorize_flippers_multi.py`: a version of `categorize_by_match_multi.py` that also checks for flippers.
If the active region "flips" in g generations, it tries to place the transformed catalysts, and then sees if the active region
reoccurs in its original location after g more generations. There is no matchType argument: it essentially works as if matchType is `same`.

Usage `python3 categorize_flippers_multi.py [rleToMatch] [start]-[end] results1.rle...resultsN.rle`. 