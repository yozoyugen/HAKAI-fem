# HAKAI-fem

![cae_20240828_001_ss](https://github.com/user-attachments/assets/bb749479-2e07-4cfb-8251-5dac315f9dfa)


## HAKAI is ...  
Simple finite element method dynamic explicit solver  
There are Matlab(.m) and Julia(.jl) versions  

## Overview  
Input: .inp file. You can use Abaqus Learning Edition for Pre  
Output: .vtk file. You can use ParaView for Post  
Windows only  
  
## Specifications of HAKAI  
Analysis type: Dynamic explicit only  
Element type: Solid hex only  
Contact: All other instances  
Element deletion: Simply judge by plastic strain  
Time step: Fixed  
Mass scaling: Fixed (user input)  
  
## Installation  
No need to install  
Just DL "HAKAI-vx.x.x" folders to your PC  
  
## Execution example  
### Matlab  
Matlab command window  
res = HAKAI('input\\Tensile5e.inp')  
  
### Julia  
Command prompt  
julia HAKAI_j.jl input\Tensile5e.inp  
  
## Demo video
Tensile5e.inp, Tensile5e-glmsh.inp  
https://youtu.be/7pOkPujlj1I


