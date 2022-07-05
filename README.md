# Microkinetics-Guided Machine Learning Pathway Search (MMLPS)

## Requirements:

* conda (recommand)
* python 3 (py3.9 tested)
* numpy
* scipy
* mpmath
* pprint
* networkx
* MongoDB and pymongo (optional)

## Usage

1. Add the `$XXX/pymodule` to your `$PYTHONPATH`

2. Run a single branch of MMLPS:

    The following directory structure should be prepared, as shown in the example
    ```
    job1
    ├── sourcedir/
        ├──input
        ├──input_split
        ├──CuZnCHO.pot
    ├── custom_para.py
    ├── pathsample.py
    ├── start.arc
    ```
    
    * `input`:  Parameters for SSW-RS, which sample possible reaction pair (`lasp.in` file)
    * `input_split`: Parameters for DESW, which identify transition state of possible reaction pair (`lasp.in` file)
    * `CuZnCHO.pot`: G-NN potential file
    * `custom_para.py`: Parameter for MMLPS
    * `pathsample.py`: Main program
    * `start.arc`: Starting structure for this branch. It should contain a slab and several molecules.
    
    To run the simulation, run the following command：
    ```bash
    python pathsample.py
    ```
    or you can submit it to the job queue with `jobs.sh`.
    
    An `analyze` directory should be built by `pathsample.py`, which contain all reaction pair sampled.

3. Ananlyze a single branch of MMLPS:
   
   The analyzation is carried out by `pathanalyze.py` in `analyze` directory
   
   To obtain the best reaction pathway, run:
   ```bash
   pathanalyze.py -surface -readallmin 2 -LGibbs 500 -ly ReacNProd.arc
   ```
