BEGIN BULK
GRID,5,,2.0,0.0
GRID,6,,2.0,1.0
GRID,7,,3.0,0.5
CBEAM,3,1,5,7
CBEAM,4,1,6,7
PBEAML,1,1,MSCBML0,L,,,,,0.1,0.4,0.018,0.013,,YES,1.0,0.1,0.4,0.018,0.013
MAT1,1,2.0e11,8.0e9,0.3,5000
FORCE,1,7,,500000.0,,,1.0
ASET1,1,5,6
ENDDATA
