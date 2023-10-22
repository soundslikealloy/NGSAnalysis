# _In-silico_ Next Generation Sequencing analysis (IN-NGSa)

*Contributors: Eloi Martinez-Rabert

_In-silico_ Next Generation Sequencing analysis (IN-NGSa) is a Python-based script that performs the simulation of Next Generation Sequencing (NGS) of a group of genomic sequencies (from whole genomes, partial genomws or only few sequences). With  In addition, it is able to examine the ability of 16S ribosomal RNA gene and ITS regions to differentiate microbes at 'species and/or strain level'. For it, the number of mismatches between the (unique) amplicons are counted as an indicator of the differentiation ability. A summarized workflow of IN-NGSa is shown below.

![IN-NGSa](https://github.com/soundslikealloy/NGSAnalysis/assets/81569132/64c29b56-8e1e-4ce3-862f-36456384923a)
____________________________

## Before having fun...
**:warning: To open the links in a new tab: right click on the link + "Open link in new tab". :warning:**

### :gear: Anaconda Python installation
This code is built up in Python. To execute Python scripts is recommended the installation of **Anaconda**. **Anaconda Python** is a free, open-source platform that allows to write and execute code in the programming language Python ([Python Tutorial](https://docs.python.org/3/tutorial/index.html)). This platform simplifies package installation, managment and development, and alos comes with a large number of libraries/packages that can be you for your projects. To install **Anaconda**, just head to the [Anaconda Documentation website](https://docs.anaconda.com/free/anaconda/install/index.html) and follow the instructions to download teh installer for your operating system.

### Anaconda Navigator
Anaconda Navigator is a desktop graphical user interface that allows you to launch applications and efficiently manage conda packages, environments, and channels without using command-line commands. For more info, click [here](https://docs.anaconda.com/free/navigator/).

### Anaconda Prompt or Terminal
Anaconda Prompt is a command line interface with Anaconda Distribution. Terminal is a command line interface that comes with macOS and Linux. To open it in **Windows**: Click Start, search for _"Anaconda Prompt"_ and click to open. In **macOS**: use Cmd+Space to open Spotlight Search and type _"Navigator"_ to open the program. In **Linux-CentOS**: open Applications > System Tools > Terminal.

### Spyder
Spyder is a Python development environment with many features for working with Python code, such as a text editor, debugger, profiler, and interactive console. You can execute **Spyder** using the **Anaconda Navigator**. You can find Spyder Tutorials [here](https://www.youtube.com/watch?v=E2Dap5SfXkI&list=PLPonohdiDqg9epClEcXoAPUiK0pN5eRoc&ab_channel=SpyderIDE).

### Python packages
A **Python package** is a collection of files containing Python code (i.e., modules). To execute **IN-NGSa**, the following packages must to be installed:
- **Pydna.** Pydna is a Python package providing code for simulation of the creation of recombinant DNA molecules using molecular biology techniques. Pydna provides simulation of Primer design, PCR, Restriction digestion and so on. For more info and tutorials, click [here](https://pydna.readthedocs.io/index.html).
- **Pyteomics.** Pyteomics is a collection of lightweight and handly tools for Python taht help to handle various sorts of proteomics data (such as FASTA files). For more info and tutorials, click [here](https://pyteomics.readthedocs.io/en/latest/).
- **Matlabplotlib.** Matplotlib is a comprehensive library for creating static, animated and interactive visualizations in Python. Fore more info and tutorials, click [here](https://matplotlib.org/).
- **Biopython.** Biopython is a set of freely available tools for Computational Molecular Biology written in Python. For more info and tutorials, click [here](https://biopython.org/).
- **Numpy.** NumPy is the fundamental package for scientific computing in Python. It is a Python library that provides a multidimensional array object, various derived objects (such as masked arrays and matrices), and an assortment of routines for fast operations on arrays, including mathematical, logical, shape manipulation, sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic statistical operations, random simulation and much more. For more info and tutorials, click [here](https://numpy.org/).

### Installation of packages using Anaconda Navigator
You can install any Python package using the **Anaconda Navigator**. For this, execute the navigator and click to **Environments**. In this section you can install new packages and delete the already installed. For more info, click [here](https://docs.anaconda.com/free/navigator/).

### Installation of packages using pip
**pip** is the package installer for Python. In general, pip installs the minimal instalation requirements automatically, but not the optionals requirements. To execute **IN-NGSa**, the minimal requirements are adequate. To install the mentioned packages using pip, you have only to write the following command lines in **Anaconda Prompt or Terminal**:

**Pydna**:
```
pip install pydna
```
**Pyteomics**:
```
pip install pyteomics
```
**Matlabplotlib**:
```
conda install matplotlib
```
**Biopython**:
```
pip install biopython
```
**NumPy**:
```
pip install numpy
```
## :clipboard: Instructions to Download and Run IN-NGSa
1. Download .zip code. Last version: `v1.0.0`. [Download code]([https://github.com](https://github.com/soundslikealloy/NGSAnalysis/archive/refs/tags/v1.0.0.zip)).
2. Extract files to a destination (:bulb: recommendation - Desktop).
3. Open **Anaconda Prompt or Terminal**.
4. Go to the **Code folder<sup>2</sup>** using `cd` command (more info about [Using Terminal](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/use-terminal/?highlight=Using%20Terminal)).
    &#09;<br><sup><sup>2</sup>Code folder: folder with `NGSanalysis.py` file. </sup>
5. Copy your GenBank (Standard) files (_.gb_) to the `in_gb/` folder, and [FASTA files](https://en.wikipedia.org/wiki/FASTA_format) (_.fa_) to the `in_fasta/` folder.
6. Select the feature that you want to test writing the name of its _/product_ in `feature/feature.txt` (be sure that it's the same as in the GenBank file). For example:
   ```
   16S ribosomal RNA
   ```
   If you want to amplify a spacer region (i.e., a region of non-coding DNA between genes), in `feature/feature.txt` write:
   ```
   spacer
   Front product of spacer region
   Back product of spacer region
   ```
   For example, bacterial ITS region between 16S ribosomal RNA and 23S ribosomal RNA:
   ```
   spacer
   16S ribosomal RNA
   23S ribosomal RNA
   ```
7. Write the primers in `primers/main_PrimersList.txt`, using the first line for the forward primer and the second line for the reverse primer. For example, Illumina primers for v3v4 of 16S ribosomal RNA:
   ```
   ACTCCTACGGGAGGCAGCAGT
   GGACTACCAGGGTATCTAATCCTGT
   ```
8. Execute **IN-NGSa** with the command line:
   ```
   python NGSanalysis.py
   ```
   **Optional arguments:**
   <table border="0">
       <tr><td>-h, --help</b></td><td>Show help message and optional arguments.</b></td></tr>
       <tr><td>-unique</td><td>Only unique amplicons are considered in the analysis.</td></tr>
       <tr><td>-nomismatches</td><td>No mismatches analysis is performed.</td></tr>
       <tr><td>-onlymismatches</td><td>Only mismatches figure is saved. FASTA and alignment files are not saved.</td></tr>
   </table>
   
   ```
   python NGSanalysis.py -h
   python NGSanalysis.py -unique
   python NGSanalysis.py -nomismatches
   python NGSanalysis.py -onlymismatches
   python NGSanalysis.py -unique -nomismatches
   python NGSanalysis.py -unique -onlymismatches
   ```
9. Get the results of **IN-NGSa**.
   <br>Amplicons are saved in `out_fasta\` folder as _fa (date).txt_ or _fa_unique (date).txt_.
   <br>Alignment results using MUSCLE[^1] is saved in `out_fasta\` folder as _fa_alignment (date).txt_ or _fa_unique_alignment (date).txt_.
   <br>Figure showing the mismatches between amplicons is saved in `out_misAnalysis\` folder as _alignmentmismatches (date).png_ or _alignmentmismatches_unique (date.png)_.
   <br>Table showing the number of mismatches is also saved in `out_misAnalysis\` folder as _mismatchTable (date).txt_ or _mismatchTable_unique (date).txt_.

## Contact

**Eloi Martinez-Rabert**. :envelope: martinezrabe@wisc.edu or eloi.mrp@gmail.com

### References
[^1]: Edgar, RC (2021). *Nature Reviews Microbiology*. doi: 10.1038/s41467-022-34630-w<br>
