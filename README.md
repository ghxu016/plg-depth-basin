# Dataset Overview

This dataset contains model data and supporting files for analyzing peak ring materials. Below is the structure and description of the dataset contents:

## Folder Structure
- **C30** and **C40**: Contain model data for different crustal thickness scenarios.

## Files and Data Description

1. **Model Information:**
   - **`modelx`** and **`modely`**: Represent the spatial location of each cell in the model.
   - **`mat`**: Includes data for five time steps, indicating the material type of all cells at each time point.

2. **Tracer Data:**
   - Each cell contains a single tracer.
   - **`xmark`** and **`ymark`**: Indicate the X (radial) and Y (depth) positions of tracers at five different times.
   - **`mass`**: Specifies the mass of each tracer.
   - **`TrP`**: Records the peak shock pressure experienced by each tracer during the impact.

3. **Indexes:**
   - **`RingIndex`**: Index of tracers in peak ring materials subjected to shock pressures < 25 GPa.
   - **`RingIndexAll`**: Index of all tracers in peak ring materials.
   - **`25Index`**: Index of tracers in materials experiencing shock pressures > 25 GPa.

4. **Basin Data:**
   - **`RingDepthBasin`**: Contains basin-specific information relevant to the study.

5. **Python Files:**
   - **`30v45_10v15_PRComp`**: A Python script for comparing four models to visualize the movement of peak ring materials.  
   - **`MatComp_TrP`**: A Python script for comparing two models to analyze the formation of peak ring and the peak shock pressure they experienced.

## Additional Notes
- The **input files** required for running models are located within the dataset folders.
