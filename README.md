# Bias-in-scaling-based-glacier-estimates

## Folder structure:
1. Steady Glaciers
..* Code to get the steady state glacier
```python
python tune.py
```
2. Transient Glaciers
..* Code to get the Glacier State with linear ELA change
..1. bed
..* bed files of glaciers
..2. ELA_CHANGE
..* Run the code for linear ELA change
```python
python change_ela.py
```
..3. STEADY_STATE
..* Run the code for obtaining steady state glacier
```python
python steady.py
```
..4. output_steady
..* output file (Volumen Area files) of steady code are stored here
..5. output_change_ela
..* output file (Volumen Area files) of ela change code are stored here

## File headers
1. tau.txt headers : glacier ID, tau, ela, beta, C
2. guess_ela.txt : glacier ID, ela, beta, C

#### note: please check the eigen3 library path on your computer. if needed make changes to MakeFile.


