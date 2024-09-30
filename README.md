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
[![''](https://github.com/user-attachments/assets/0954b6b4-574d-4911-a83a-97b43eabf731)](https://youtu.be/7pOkPujlj1I)
  
bullet-impact.inp, bullet-impact-glmsh.inp, 
[![''](https://github.com/user-attachments/assets/46943222-3946-465b-9b9e-a4c64aeb565b)](https://youtu.be/RczUhYBTxVg)  
  
Tensile-test.inp, Tensile-test-glmsh.inp (Dog-bone)  
[![''](https://github.com/user-attachments/assets/430d8c8f-69ab-41f7-9ad8-02c56e010f8b)](https://youtu.be/zrBSUzXMs58)  

  




