Supplementary files for the article accepted to the Journal of Cheminformatics:

Atomic Ring Invariant and Modified CANON Extended Connectivity Algorithm for Symmetry Perception in Molecular Graphs
and Rigorous Canonicalization of SMILES.

Dr. Dmytro Krotko, Enamine Ltd., 2020. E-mail: d.krotko@gmail.com

PL/SQL source codes for the SMILES parsing and the calculation of the local and nonlocal atomic invariants (smilin.txt),
the Modified CANON algorithm (canon.txt), the deep-first traversal of the molecular graph (cangen.txt) and the Breaking
Ties and the selection of the lexicographically minimal canonical SMILES (canonsmi.txt) are available.

All types, procedures and functions must be compiled in the order in which they are ordered in the package description
file ChemInf.txt.

The script smilin.sql is intended for detailed analysis of any molecular graph and establishment of the local and nonlocal
invariants of its atoms and bonds. Also the script shows the results of canonical code generation by the automorphism
permutation for the complex molecular graph.
