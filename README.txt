To run Hit Removal you will need to execute two scripts in succession.

The first applies the CosmicRemoval algorithms:

    python $LARLITE_USERDEVDIR/HitRemoval/CosmicRemoval/mac/cosmicremoval.py $INPUT

The output will be created in your current directory: cosmicremoval.root, referred to as $OUT1

The second applied the NeutrinoTrackRemoval algorithms:

    python $LARLITE_USERDEVDIR/HitRemoval/NuTrackRemoval/mac/nutrkremoval.py $OUT1

This will create the file nutrkremoval.root in your current directory.

To run with the default settings you will need the following producers:

hits :     gaushit
clusters : pandoraCosmic
vertex :   mcvertex

These can be edited in the python scripts referenced above.
