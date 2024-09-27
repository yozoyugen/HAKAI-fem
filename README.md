# HAKAI-fem
## HAKAI is ...  
Simple finite element method dynamic explicit solver  
There are Matlab(.m) and Julia(.jl) versions  

## Overview  
Input: .inp file. You can use Abaqus Learning Edition for Pre  
Output: .vtk file. You can use ParaView for Post  
  
## Specifications of HAKAI  
Analysis type: Dynamic explicit only  
Element type: Solid hex only  
Contact: All other instances  
Element deletion: Simply judge by plastic strain  
Time step: Fixed  
Mass scaling: Fixed (user input)  
