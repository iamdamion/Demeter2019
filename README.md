## <p align="center">[Functional connectivity fingerprints at rest are similar across youths and adults and vary with genetic similarity](http://www.sciencedirect.com/science/article/pii/S2589004219305474)</p>
[doi.org/10.1016/j.isci.2019.100801](https://doi.org/10.1016/j.isci.2019.100801)  

### Authors:
Damion V. Demeter*, Laura E. Engelhardt, Remington Mallett, Evan M. Gordon, Tehila Nugiel, K. Paige Harden, Elliot M. Tucker-Drob, Jarrod A. Lewis-Peacock, Jessica A. Church
* Correspondence: demeter@utexas.edu

![Graphical Abstract](https://ars.els-cdn.com/content/image/1-s2.0-S2589004219305474-fx1_lrg.jpg)

---
<p align="center"> Scripts and other resources for replication of this work</p>   
### Scripts:
- Main SVM Classifier Script: Twin_Match.py   
- SVM Using Opposite Group Features Script: Opposite_Mask_SVM.py   

- Requirements: (Unfortunately) these scripts are written in python 2.7. Other requirements:
```
numpy
scipy
sklrean
pandas
seaborn
matplotlib
```
- Timecourse files should be in a group folder, with each pair given a family ID, then an individual ID separated with an underscore. Example with two "family" pairs. Family 201 and family 202:
```
/MZ_TWINS/
         /201_1.txt
         /201_2.txt
         /202_1.txt
         /202_2.txt
         /etc
```
- Beyond this, see the -h (help) argument in the scripts for full details of required arguments. 
- See in-code comments for notes about each step, etc. (I try to heavily comment my code....probably too much to be honest).

---
### Required Playlist:
The following playlist was used for data cleaning. The main author of this work can't guarantee your successful use of these scripts without this playlist:   
https://open.spotify.com/playlist/6gDuLpHIlHqsU58J7eARkS?si=C5bKBQSQTK2O2mU14LfaMQ



