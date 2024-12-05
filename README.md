# A coarse-grained model of disordered RNA for simulations of biomolecular condensates

This repository contains codes, raw data and parameter files that used to reproduce the results in "A coarse-grained model of disordered RNA for simulations of biomolecular condensates by Yasuda, I., von Bülow, S., Tesei, G., Yamamoto E., Yasuoka K. & Lindorff-Larsen, K. (2024)".

### Layout

- `calvados/` Codes and installation guide to run CALVADOS 2-RNA simulation as implemented in this work.  
- `analysis/` Scrtipts to run analysis 
- `all-atom/` All-atom structures of UCAAUC from [Bergonzo ,C., Grishaev, A., Bottaro, S., RNA, 28, 937–946 (2022)], and HCG_models from [Pietrek, L. M., Stelzl, L. S., Hummer, G., J. Chem. Theory Comput., 20, 2246–2260 (2024).]
- `single-chain/` Input files into calvados for single chain simulations of polyR30 
- `flory_exponent/` Input files to run single chain simulations of polyR10-100
- `400-chain/` Input files into calvados for 400-chain simulations of polyR30
- `FUSRGG3/`  Input files into calvados for FUSRGG simulations of polyR40
- `MED1/` Input files into calvados for MED1 simulations of polyR40
- `viewer1.ipynb` Jupyter Notebook to produce plots for the paper figures (Figure 2).
- `viewer2.ipynb` Jupyter Notebook to produce plots for the paper figures (Figures 3 and 4).
- `data/` Files to generarte figures

Simulations can be run using `prepare_full.py`
```
python prepare_full.py
cd {name}/{temperature}/
python run.py
```

Continuing simulation can be run in the same directory of previous simulation,
```
python run.py
