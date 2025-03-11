
# Semidefinite Programming (SDP) Problem Generator and Solver

This project provides a Python-based implementation for formulating and solving a Semidefinite Programming (SDP) problem using the `cvxpy` library. It constructs matrices representing the SDP formulation and generates a textual representation of the constraints.

All theoretical explanations and detailed formulations can be found in `moment_sos_hierarchy.pdf`.

* [Features](#-features)
* [Dependencies](#-dependencies)
* [Usage](#-usage)
* [Function Descriptions](#-function-descriptions)
* [Output](#-output)
* [Example Output](#-example-output)
* [License](#-license)

## 📌 Features

- Constructs symmetric positive semidefinite matrices `Q0`, `Q1`, and `Q2`.
- Maximizes a scalar variable `λ (lambda)` subject to semidefinite and equality constraints.
- Automatically writes the symbolic SDP formulation to a text file.
- Modular code to integrate custom polynomial bases and matrix dimensions.

## 🔧 Dependencies

- Python 3.x
- [cvxpy](https://www.cvxpy.org/)
- numpy

You can install the dependencies using:

```bash
pip install cvxpy numpy
```

## 🚀 Usage

You can run the solver using:

```bash
python sdp_solver.py
```

You may modify the dictionary of monomials and `help1`, `help2`, `help3`, `help4` arrays to reflect your custom basis and constraints.

### Example (inside the script):

```python
print('RELAXATION (SOS)3')

# Monomials used in B0(x)
dico_monomials = {
    '1': 0,
    'x1': 1, 'x2': 2, 'x3': 3,
    'x1**2': 4, 'x2**2': 5, 'x3**2': 6,
    'x1*x2': 7, 'x1*x3': 8, 'x2*x3': 9,
    ...
}

# Lists representing B1(x) multiplied by x1, x2, x3 and B1(x) alone
help1 = [...]
help2 = [...]
help3 = [...]
help4 = [...]

# Matrix sizes
d0 = ...
d1 = ...
d2 = ...

# Run the SDP solver
SDP(dico_monomials, help1, help2, help3, help4, d0, d1, d2, "sdp_problem.txt")
```

## 🧠 Function Descriptions

### `SDP(...)`

Formulates and solves the semidefinite programming problem by:
- Creating matrix variables `Q0`, `Q1`, `Q2`
- Adding polynomial constraints
- Generating a symbolic form of all constraints
- Solving the problem using `cvxpy`
- Printing results and writing constraint formulation to a `.txt` file

### `generation_SDP_texte(C, filename)`

Creates a readable textual representation of the constraints used in the SDP problem.

#### Parameters:
- `C`: List of symbolic constraints
- `filename`: Output text file name

## 📄 Output

After running the script, the output will include:
- The optimal value of lambda
- The values of matrices `Q0`, `Q1`, and `Q2`
- The status of the optimization problem
- A text file (`sdp_problem.txt`) with the constraint formulation

## 📚 Example Output

```
RELAXATION (SOS)3
Optimal value of lambda: 0.84321
Q0 =
[[1.2, 0, ...], [...]]
Q1 =
[[...]]
Q2 =
[[...]]
Problem status: optimal
sdp_problem.txt has been generated.
```

## 📃 License

This project is open-source and free to use under the [MIT](https://github.com/sepanta007/Moment_SOS_Hierarchy/blob/master/LICENSE) License.