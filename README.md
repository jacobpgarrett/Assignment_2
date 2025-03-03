# ME 700 Assignment 2
This assignment serves as an introduction to Matrix Structural Analysis using the direct stiffness method.  It incorporates nonlinear analysis tools to perform elastic critical load analysis and plots the deformed structure of the initial geometry.  I would greatly appreciate any feedback on the functionality/ease of use and readability.

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/jacobpgarrett/Assignment_2/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/jacobpgarrett/Assignment_2)
[![tests](https://github.com/jacobpgarrett/Assignment_2/actions/workflows/tests.yml/badge.svg)](https://github.com/jacobpgarrett/Assignment_2/actions)
---

## Setup
To utilize this package, download the repo using the following command:

```bash
git clone https://github.com/jacobpgarrett/Assignment_2
```

Then navigate to the associated directory.  Then, begin entering the command below to load miniconda:

```bash
module load miniconda
```

Then, set up a mamba environment
```bash
conda create --name MSA python=3.12
```

And then activate said environment:
```bash
conda activate MSA
```

Double-check that the correct version of Python is installed
```bash
python --version
```

Ensure that pip is using the most up-to-date version of setuptools:
```bash
pip install --upgrade pip setuptools wheel
```

Create an editable install of the code:
```bash
pip install -e .
```

Then, you can run pytest to ensure proper code coverage:
```bash
pytest -v --cov=src --cov-report term-missing
```

## Tutorials
There is a pair of tutorials for the Matriux Structural Analysis algorithm in the tutorials folder saved as Jupyter Notebooks.  These tutorials show an example implementation of the nonlinear analysis and matrix structural analysis scripts.  To view geometric stiffness matricies, critical load factors, or other parameters pertaining to nodes or elements, the user can index into the nodevals and elementvals vectors.  To do so, execute the following commands:

```bash
print(elementvals[i].k_g) # Geometric Stiffness matrix
print(elementvals[i].lambda_crit) # Elastic Critical Load Factor
```

# General Information
## Generative AI Use

In this assignment, the forms of Generative AI that were used were vsCode Copilot and ChatGPT, both of which assisted in code writing.

## Contributing
Pull Requests are welcome

## License
[MIT](https://choosealicense.com/licenses/mit/)
