CalEst 1.0 — Structural Calculation Engine
CalEst 1.0 is a structural calculation system built with Borland C++ 5.5.
It reads structural models from .dat input files, performs the analysis, and produces a resolution output file containing the computed results.
A complementary project provides a graphical IDE for structural analysis—allowing you to draw the structure, define loads/forces, and then export/run the analysis with CalEst 1.0.

✨ Key Features

* Engineered in Borland C++ 5.5 for classic, fast native execution.
* File‑driven workflow: ingest .dat model files, compute, and emit a resolution output file.
* Deterministic analysis suitable for versioning and reproducible runs.
* Decoupled UI: modeling can be done with a separate IDE (another GitHub repo) that supports graphical capture of members, nodes, boundary conditions, and forces.

🧩 Project Structure

CalEst/
  src/
    CalEst.h
    CalEst.cpp
    Solvers.h
    Solver.cpp
    EdoError.h
    EdoError.cpp
    AlgLineal.h
    AlgLineal.cpp

📥 Input Format (.dat)
CalEst 1.0 expects a plain‑text .dat file describing:

Nodes (IDs and coordinates)
Elements/Members (connectivity, section properties, material properties)
Boundary Conditions (supports, restraints)
Loads/Forces (nodal forces, distributed loads, load cases)
Analysis Options (if applicable)


Tip: Keep units consistent across geometry, loads, and materials.

Example (illustrative only):
📥 Input Format (.dat)
CalEst 1.0 expects a plain‑text .dat file describing:

Nodes (IDs and coordinates)
Elements/Members (connectivity, section properties, material properties)
Boundary Conditions (supports, restraints)
Loads/Forces (nodal forces, distributed loads, load cases)
Analysis Options (if applicable)


Tip: Keep units consistent across geometry, loads, and materials.

Example (illustrative only):
📥 Input Format (.dat)
CalEst 1.0 expects a plain‑text .dat file describing:

Nodes (IDs and coordinates)
Elements/Members (connectivity, section properties, material properties)
Boundary Conditions (supports, restraints)
Loads/Forces (nodal forces, distributed loads, load cases)
Analysis Options (if applicable)


Tip: Keep units consistent across geometry, loads, and materials.

Example (illustrative only):
📥 Input Format (.dat)
CalEst 1.0 expects a plain‑text .dat file describing:

Nodes (IDs and coordinates)
Elements/Members (connectivity, section properties, material properties)
Boundary Conditions (supports, restraints)
Loads/Forces (nodal forces, distributed loads, load cases)
Analysis Options (if applicable)


Tip: Keep units consistent across geometry, loads, and materials.

Example (illustrative only):

PreMECA 1.0
Viga_continua
CM	KG
3	2	3	1	2	3	1	1	1	2	0	0	0
1	1	2
1	2	3
0	0	
500	0	
1500	0	
1	1	1	0	0	0	0	
2	1	1	0	0	0	0	
3	1	1	0	0	0	0	
1	1	1	1	
Caso_sin_título
0	2	2	0	0	0
1	-10000	250
2	-10000	250
1	-500
2	-500

📤 Output (Resolution File)
CalEst 1.0 writes a resolution output file (text) summarizing:

Global input echo (nodes/elements/materials)
Stiffness assembly and analysis options (optional)
Displacements at nodes
Element forces and reactions
Basic checks and convergence information (if applicable)

Sample excerpt:

ProMECA 1.0
Viga_continua
CM	KG
3	2	3	1	2	3	1	1	1	2	0	0	0
1	1	2	
1	2	3	
0.000000	0.000000	
500.000000	0.000000	
1500.000000	0.000000	
1	1	1	0	0.000000	0.000000	0.000000	
2	1	1	0	0.000000	0.000000	0.000000	
3	1	1	0	0.000000	0.000000	0.000000	
1.000000e+00	1.000000e+00	1.000000e+00	1.000000e+00	
Elemento 1
0.002000	0.000000	0.000000	-0.002000	0.000000	0.000000	
0.000000	0.000000	0.000024	0.000000	0.000000	0.000024	
0.000000	0.000024	0.008000	0.000000	-0.000024	0.004000	
-0.002000	0.000000	0.000000	0.002000	0.000000	0.000000	
0.000000	0.000000	-0.000024	0.000000	0.000000	-0.000024	
0.000000	0.000024	0.004000	0.000000	-0.000024	0.008000	

Elemento 2
0.001000	0.000000	0.000000	-0.001000	0.000000	0.000000	
0.000000	0.000000	0.000006	0.000000	0.000000	0.000006	
0.000000	0.000006	0.004000	0.000000	-0.000006	0.002000	
-0.001000	0.000000	0.000000	0.001000	0.000000	0.000000	
0.000000	0.000000	-0.000006	0.000000	0.000000	-0.000006	
0.000000	0.000006	0.002000	0.000000	-0.000006	0.004000	

CASO_DE_CARGA:1	Caso_sin_título
0	2	2	0	0	0
 1	-10000.000000	250.000000
 2	-10000.000000	250.000000
1	-500.000000
2	-500.000000
0.000000	0.000000	1263020833.333333	
0.000000	0.000000	-5286458333.333333	
0.000000	0.000000	13177083333.333334	
0.000000	33437.500000	0.000000	
0.000000	532343.750000	0.000000	
0.000000	204218.750000	0.000000	
0.000000	33437.500000	0.000000	0.000000	226562.500000	-48281250.000000	
0.000000	305781.250000	48281250.000000	0.000000	204218.750000	0.000000	

