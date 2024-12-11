# Alya_input_setup

This pipeline handles the creation, monitoring, postprocessing, and analysis of Alya electromechanical or electrophysiological simulations. It is highly modular and adaptable to different supercomputing or cluster computing platforms, as well as to different version of Alya, and is designed to be able to deploy and analyse large volumes of simulations. 

Code structure is summarised in the following flowchart:
```mermaid
graph TD;
    A[Input geometric and field data e.g. vtk format]-->B(Mesh Preprocessing)
    B --> C[Standard .CSV format geometric and field data];
    B --> |Optional| D(Generate additional fields e.g. fibre vectors, scars);
    D --> C;
    C -->|Optional| E(Digital Twin inference pipeline);
    E --> F;
    F[Personalised electrophysiology] --> C;
    C --> G(Alya files generation);
    H[Version-specific Alya templates] --> G;
    I[Cluster-specific job submission templates] --> G;
    J[Parameters set json file] -->G;
    K(Sampling methods)-->J;
    G --> L[Alya simulation files];
```


