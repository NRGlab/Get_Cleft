# GetCleft
The Get_Cleft program detects cavities in PDB-formated macromolecular structure models. It is based on the Surfnet algorithm of Laskowski et al. Get_Cleft is used by IsoMIF as well as FlexAID and the NRGsuite.

To compile GetCleft: 
```
gcc Get_Cleft.c -o Get_Cleft -O3
```

GetCleft is also available to be installed as a python package:
https://pypi.org/project/getcleft-py/
