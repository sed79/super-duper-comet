# Data analysis of a super duper comet - 2I/Borisov

This repo contains scripts and associated files for the data analysis of VLT/MUSE observations of the interstellar comet 2I/Borisov.
Scripts are excecuted individually to perform separate steps of the analysis pipeline. They should be run in order as outlines below as scripts use files created by the previous step's script. Figures of data outputs are also created using scripts beginning with 'Display'. The pipeline begins with reduced IFU datacubes and create the following data products and associated figures.

## Dust Colour (S')
Input: reduced data cube
1. MaskingExtremes.py
2. Shift.py
3. SolarReflectanceSlopes.py
4. S'aperture.py

<br>Output: S' values for each observing night

5. DisplayDustColour.py
6. DisplayS'plotsMUSEvsLiterature.py

## Dust maps
Input: reduced data cube
1. MaskingExtremes.py
2. Shift.py
3. MapExtraction.py 
4. MedianCoaddMaps.py
5. CropMaps.py

<br>Output: .fits files of dust emission maps

6. DisplayDustGrid.py

## Gas maps
Input: reduced data cube
1. MaskingExtremes.py
2. Shift.py
3. DustSubtraction.py
4. CNSubtraction.py
5. MapExtraction.py
6. MedianCoaddMaps.py
7. CropMaps.py

<br>Output: .fits files of gas emission maps

8. DisplayGasC2Grid.py
9. DisplayGasNH2Grid.py
10. DisplayGasCNGrid.py

## Extra Figures/Tables/Data
Figure of 2I's solar system trajectory: BigHero.py

Figure comparing gases radial distribution: DisplaySpatialDistribution.py

Grid of enhanced dust maps: DisplayDustEnhancedGrid.py

Figure of production rates: DisplayProductionRates&BothRatios.py
Comparison of gas map with stars in frame: 1. LambdaSum.py 2. DisplayCNvsStars.py
Grid of enhanced gas maps: DisplayEnhancedGasMaps.py
Table of key observating parameters: ObservationTable.py
Calculating optocentre: 1. LambdaSum.py 2. CentroidCalc.py
