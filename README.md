# Dissertation
## The following Python packages are required to run the overall pipeline:
*numpy<br>
*pandas<br>
*biobox<br>
*matplotlib<br>
*seaborn<br>
*Modeller(Note: Installation of Modeller requires a license key: https://salilab.org/modeller/registration.html. Note: there are issues with Modeller version 10.3, but any more recent version works correctly.)<br>

There are diferent reauirements for the methods for calculating diferent featres for the vsines. The table below documents the dependencies for each.
| Feature  | Methods  | Requirements |Where to find |
|----------|--------- |--------------|--------------|
| **SASA** | Biobox   | Biobox       |anaconda      |
| **Depth**| Biopython| Biopython    |anaconda      |           
| **pKa**  | PROPKA3  | propka       |anaconda      |

## A few points to note:
1.**`Raw data extraction.ipynb`** is mainly used to extract structural data from `PDB` and `AlphaFold` databases and calculate feature data. When running, the `Uniprot_Entry` column in the **`glycation.csv`** file needs to be used as input. Note: Duplicate proteins should be deleted.<br>

2.**`Alphafold.py`**, **`measure.py`**, **`protein.py`**, **`uniport.py`** and other files are python files needed in the extraction process.<br> 

*The **`Uniprot`** class is used to handle collecting the structural information for the proteins required by mining the UNIPROT database. <br>

*The **`Protein`** class is used to handle extracting PDB files for proteins specified within the `Uniprot` class. Structures are downloaded and patched to ensure good quality structures as used for calculations.<br>

*The **`Alphafold`** class is used to handle extracting AlphaFold data for proteins specified within the `Uniprot class`. Structures are downloaded and patched to ensure good quality structures as used for calculations.<br>

*The **`Measure`** class has been implemented to facilitate the addition of new measurable features.<br>

3.**`Data Processing and Machine Learning.ipynb`** is the main model code of this research, which is used to process raw data, build machine learning models, and visualize them.<br>

4.The **`EXP dataset`** and the **`PAL dataset`** are the datasets used in this research.<br>
