# Viewing repeats in consensus maps

Instant *de novo* detection and visualization of 
repeat structures from optical mapping data

## Opening a repeat map web view via command line
- `python3 -m omacase -m view input.cmap -r`
- The `-r` option filters viewing to only repeat-containing maps 
and highlights the repeats. It can be unset and reverted interactively 

## Opening in the OMACASE web UI
1. Open the OMACASE Web app
2. Go to the **"Upload Dataset"** tab. 
   Upload the file into the storage folder if it is not already present.
3. Go to the **"Select Uploaded Dataset"** tab. Select the file to analyze.
4. Go to **"Map Viewer"**. Refresh once (as needed due to Dash behaviour).
5. Select **"Repeat Only"** under "Other Filtering" method

Maps can be sorted by the IDs, lengths, and repeat scores of maps.

The repeat score is calculated based on the number of distinct segments 
in a repeat unit, copy number, and length variance of segments in a repeat unit.
The score is used to rank repeat candidates according to their potential 
level of interest, as the more distinct segment a repeat unit possesses, 
the more likely it can be aligned to an *in silico* reference or other maps 
for location and identification. The higher copy number also presents more 
potential functional impact and hence a higher interest for study.

The OMACASE view module also provides page and view region control.
Map with zooming in/out capability.