# Alya_input_setup

This pipeline handles the creation, monitoring, postprocessing, and analysis of Alya electromechanical or electrophysiological simulations. It is highly modular and adaptable to different supercomputing or cluster computing platforms, as well as to different version of Alya, and is designed to be able to deploy and analyse large volumes of simulations. 

Code structure is summarised in the following simplified flowchart:
```mermaid
graph TD;
    A(Input geometric and field data) --> |Mesh Preprocessing|B(Standard .CSV format dataset);
    B --> |Fields generation|C(Additional fields);
    B --> H;
    H(Personalisation) -->|Inference| C;
    B --> |Alya formatting| D(Alya simulation files);
    C --> |Alya formatting| D;
    E(Sampling methods) --> |Generation| F(Parameters set json files);
    F --> |Alya formatting| D;
    D --> |Postprocessing| G(Simulated biomarkers);
```
And a more comprehensive flowchart below:

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
    L --> M(Submit multiple jobs to queue);
    N(Live monitoring scripts) --> M;
    M --> O(Postprocessing);
    O --> P[ECG biomarkers];
    O --> Q[Pressure volume biomarkers];
    O --> R[Displacement biomarkers];
    O --> S[Strain biomarkers];
    P --> T(Sensitivity analyses);
    Q --> T;
    R --> T;
    S --> T;
    P --> V(Evaluation);
    Q --> V;
    R --> V;
    S --> V;
    U[Healthy biomarker ranges] --> V;
    
```


