# generate tutorial data example input
python generate_example_inputs.py

# SynchDP alignment with reference sequence guided
python synchdp.py  -m reference -q ../tutorial/input/querry_set.pkl  -t ../tutorial/input/reference_seq_set.pkl  --date_info ../tutorial/input/date_info.pkl  --omics_info ../tutorial/input/omics_info.csv  -p 0.1 -c 0.65 -w 6 -e 3 -o ../tutorial/out

# SynchDP alignment with De novo reference sequence
python synchdp.py  -m denovo  -q ../tutorial/input/querry_set.pkl  --date_info ../tutorial/input/date_info.pkl  --omics_info ../tutorial/input/omics_info.csv  -r 95 -p 0.01 -c 0.7 -w 6 -e 3.5 -o ../tutorial/outdenovo

