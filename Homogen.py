# Homogenous solution

import numpy as np
from scipy.linalg import eig
import sympy as sp

def analyse_system(A):
    t = sp.symbols('t', real=True)                  # Define symbolic time
    eigenvalues, eigenvectors = eig(A, right=True)  # Calculate Eigenvalues and -vectors

    # Adjust eigenvectors to scale the first non-zero element to 1
    for i in range(eigenvectors.shape[1]):
        vec = eigenvectors[:, i]
        non_zero_index = np.nonzero(vec)[0][0]      # Find the first non-zero index
        eigenvectors[:, i] /= vec[non_zero_index]   # Scale by the first non-zero element

    # Handle multiplicity
    unique_eigenvalues, unique_indices, Total = np.unique(eigenvalues, return_index=True, return_counts=True)

    print("Eigenvalues and their eigenvectors, respectively:")
    ci_count = 1  # Counter for constants
    Solution_terms = []
    stability = "stable"
    constant_terms = []  # To store constant descriptions
    for i, eigenvalue in enumerate(unique_eigenvalues):
        multiplicity = Total[i]
        lambda_real = round(np.real(eigenvalue), 3)  # Real part of the eigenvalue
        lambda_imag = round(np.imag(eigenvalue), 3)  # Imaginary part of the eigenvalue
        eigenvector = eigenvectors[:, unique_indices[i]]
        Rounded_eigenvector = tuple(np.round(np.real(v), 3) for v in eigenvector)  # Round each component

        print(f"Eigenvalue {i+1} ({np.round(eigenvalue, 3)}): Multiplicity {multiplicity}, Eigenvector {i+1}: Vector {Rounded_eigenvector}")

        if lambda_real > 0:
            stability = "unstable"
        elif lambda_real == 0:
            if stability != "unstable":
                stability = "stable, but not asymptotically stable"
        else:
            if stability != "unstable" and stability != "stable, but not asymptotically stable":
                stability = "asymptotically stable"

        # Generate terms for the homogeneous solution
        for j in range(multiplicity):
            ci = sp.symbols(f'c{ci_count}', real=True)
            # Handle t power formatting
            if j > 1:
                t_potens = f"t**{j}"
            elif j == 1:
                t_potens = "t"
            else:
                t_potens = ""

            # Form the homogeneous solution term for different cases
            if lambda_imag == 0:  # Real eigenvalue case
                term = f"{ci} * Vektor{Rounded_eigenvector} {t_potens}".strip()
                if lambda_real != 0:  # Include the exponential term if lambda_real is non-zero
                    if lambda_real == 1:  # Format as exp(t)
                        exp_term = "exp(t)"
                    elif lambda_real == -1:  # Format as exp(-t)
                        exp_term = "exp(-t)"
                    else:
                        exp_term = f"exp({lambda_real} * t)"
                    term += f" * {exp_term}"
            else:  # Complex eigenvalue case
                omega = abs(lambda_imag)
                # Format the cosine term properly for 1 or -1 cases
                if omega == 1:
                    cos_term = "cos(t)"
                elif omega == -1:
                    cos_term = "cos(-t)"
                else:
                    cos_term = f"cos({omega} * t)"
                
                # Always include the exponential term for the real part of the eigenvalue
                if lambda_real != 0:
                    if lambda_real == 1:  # Format as exp(t)
                        exp_term = "exp(t)"
                    elif lambda_real == -1:  # Format as exp(-t)
                        exp_term = "exp(-t)"
                    else:
                        exp_term = f"exp({lambda_real} * t)"
                    term = f"{ci} * Vektor{Rounded_eigenvector} {t_potens} * {exp_term} * {cos_term}".strip()
                else:
                    term = f"{ci} * Vektor{Rounded_eigenvector} {t_potens} * {cos_term}".strip()

            # Clean up formatting (remove redundant spaces, handle 1.0 coefficients)
            term = term.replace(" 1.0 ", " ").replace(" -1.0 ", "-").replace(" * ", " ").strip()
            Solution_terms.append(term)
            constant_terms.append(f"c{ci_count} (real number)")  # Add constant description
            ci_count += 1  # Increment counter per iteration

    print()  # Newline
    # Output the homogeneous solution
    solution_string = " + ".join(Solution_terms).replace("  ", " ")
    constants_description = ", ".join(constant_terms)  # Describe the constants
    print(f"Homogeneous solution: y_h(t) = {solution_string}")
    print(f"Constants: {constants_description}")  # Print constants as real numbers
    print()  # Newline
    print(f"Stability: The system is {stability}.")

# Define the matrix A with a complex eigenvalue
A = np.array([[-1 - 1j, 1 + 3j, 1], [-2, 1j, 2], [3, 3, 1j]])
analyse_system(A)