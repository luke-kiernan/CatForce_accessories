# CatForce_accessories
Helper scripts for searching for oscillators in Conways Game of Life via mvr's depth-first and periodic modifications
 of simsim's [CatForce](https://github.com/mvr/CatForce). For more information on these things, see those repos. 
 For information on how to use those programs or similar topics, see [LifeWiki](https://conwaylife.com/wiki/Main_Page) and the [conwaylife forum](https://conwaylife.com/forums/).

**Catlist Scripts**
These are scripts for generating CatForce input files. There are three main file formats: rle, csv, and plaintext. 
- csv format: the data for each catalyst is stored as entries in a table, 1 row per catalyst. Good for editing or splicing together lists.
- rle format: one catalyst per line. The leftmost entry should be the catalyst, with required cells as states 3, 4, or 5. The rest then are the forbidden patterns. Good for adding things to lists--almost all the lists in this repo were created using this script.
- plaintext format: the format accepted by CatForce.
Four python scripts that convert between these formats are included:
- `csv_to_plaintext.py` example usage `python3 csv_to_plaintext.py input.csv [name] > output.txt`. If the optional `name` parameter is provided, then the catalysts are labeled with their name (from the csv) and period. Works for peridoic and stable. 
- `plaintext_to_csv.py` example usage `python3 plaintext_to_csv.py input.txt output.csv`. Creates a column for name, but that must be entered by hand. Works for peridoic and stable. 
- `csv_to_rle.py` [Golly](http://golly.sourceforge.net) script. Columns are assumed to be: name, absence, rle, dx, dy, symType, required, req dx, req dy, forbidden 1, forbid 1 dx, forbid 1 dy, forbidden 2, etc. in that order. Works for peridoic and stable. 
- `rle_to_csv.py` [Golly](http://golly.sourceforge.net) script. Columns for name, absence, period, and symmetry type are created, but those fields must be entered by hand. Works for peridoic and stable: just delete the period column for stable. Make sure to include at least 5 cells of empty space between rows, and between patterns within each row.

**Additional Catalyst Lists**
These are some additional lists of catalysts or sparkers that I've made. Add catalysts to the rle or alter parameters in the csv, then re-run the scripts to generate new lists. Included are:
- in the `by-period` folder, there's csv's of sparkers of each period.
- `all_sparkers.txt`: all the sparkers in the by-period folder, in one big list.
- `monomers_and_variants`: a list of stable catalysts. Intended for completing partials, where space is tight and the usual dimers may not fit. Again, not comprehensive: this list omits several large catalysts on `The_Big_List_of_Catalysts`, such as the hivepush catalysts. Most of these catalysts come from [Tutorials/Catalyses](https://conwaylife.com/wiki/Tutorials/Catalyses) and [Catalyst Discussion Thread](https://conwaylife.com/forums/viewtopic.php?f=2&t=1878&p=24656&hilit=catalyst+testing#p24440).

**Bash Search Script**
A driver script that calls CatForce with different inputs and prints the results. There is no good interface currently;
 using it requires setting some variables in the script. Process:
0. Put a copy of `CatForce.cpp` and `LifeAPI.h` in the directory, either the `depth-first` or `depth-first-periodic` branch.
1. Set up catalyst lists. Follow the pattern of the `topTwo` and `Kazyan` variables to set up catalyst or sparker lists
 from within bash, or read in from a text file, like the `tinyList` variable.
2. Configure the search parameters for the input file that's created: which catalysts, start generation, etc.
3. Set meta-search parameters: which offests, regions, symmetries to iterate through.
4. Run the executable. It'll delete empty result files if they're empty, and move non-empty ones 
to a folder named after the active regions. It also prints the average time per search.

The bash script also relies on a handful of other helper scripts--these are rough things quickly thrown together to 
compute pieces of information that I needed when writing the bash script.  These are:
- `active_region_occurs.cpp`
- `compute_flipper_filters.cpp`
- `getWidthHeight.py`
- `generateOrientedRLEs.py`
The first two rely on `LifeAPI.h`, so a Makefile is included as well. (The rules for these are the same 
as for building CatForce, just different .cpp file. The bash script will run `make` if it doesn't see the executables it needs in its directory.)

I typically aim to have each search finish on the order of magnitude of a few minutes. Currently, this means that 4-catalyst
searching is limited to `topTwo` and 3-catalyst searching to `Kazyan` (or lists of similar lengths, where
we count with multiplicity, ie taking symmetry and period into account). For 2-catalyst searching, 50
catalysts is about right.

If you're having trouble with the script, set `debug=true` so that it'll print the CatForce output. 
Also take a look at the temp file created and make sure it's properly formatted. The script 
includes some more comments about how it works and such.

The majority of oscillators found so far (as of late August 2022) via CatForce oscillator search 
have been found wtih the first four active regions and C2 symmetries. For stable-supported, 3 catalysts 
with the Kazyan list has been fruitful, though I'm starting to recognize partials as ones that I've seen before.
For sparker-supported, there's loads of results with just one sparker, or one sparker and common one catalyst.

**Post Search Filtering**
The bash script simply calls CatForce on each each offset, symmetry, and active region. Results end up in separate files. 
When searching at volume, going through these raw result files by hand isn't feasible: often one mechanism will show up 
at several of offests, creating a boatload of uninteresting result files. Both of these scripts require [Lifelib](https://gitlab.com/apgoucher/lifelib/-/tree/master), and only work with stable catalysts.

`categorize_by_match_multi.py`: combine results between multiple CatForce results files, based off in which generation
a specific subpattern occurs. Decently useful for processing CatForce search results in general, even not from the bash script.

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

Another option for post-search processing is to use apgsearch's pattern recognition. mvr's [golchemy repo](https://github.com/mvr/golchemy) includes `catpipe.py`, which separates out a CatForce `.rle` file into different results, 
which then can be piped into apgsearch. However, this will only find complete oscillators: near-misses or 
flippers won't be recognized.