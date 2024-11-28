# HAKAI-fem

![HAKAI_icon_20241001_001](https://github.com/user-attachments/assets/6307a870-0d5e-4274-8eec-c5a10ec388ce)

## Web site ...  
https://caelab.work/  
  
## HAKAI is ...  
Simple finite element method dynamic explicit solver  
for elastoplastic and fracture analysis.  
There are Matlab(.m) and Julia(.jl) versions  

## Overview  
Input: .inp file. You can use Abaqus Learning Edition for Pre  
Output: .vtk file. You can use ParaView for Post  
Windows only  
  
## Specifications of HAKAI  
Analysis type: Dynamic explicit only  
Element type: Solid hex only  
Contact: All other instances(No friction)  
Element deletion: Simply judge by plastic strain  
Time step: Fixed  
Mass scaling: Fixed (user input)  
No damping  
  
## Installation  
No need to install  
Just DL "HAKAI-vx.x.x" folders to your PC  
### Demo video  
[![''](https://github.com/user-attachments/assets/35cd1b3b-4d7a-473d-9537-e63199070e8a)](https://youtu.be/XKaCQ7hLbfo)  
  
## Execution example  
### Matlab  
Matlab command window  
res = HAKAI('..\\input\\Tensile5e.inp')  
  
### Julia  
Command prompt  
julia HAKAI_j.jl ..\input\Tensile5e.inp  
  
## Example analysis  
metal-cutting.inp, metal-cutting-glmsh.inp  
[![''](https://github.com/user-attachments/assets/074c40b4-2ee5-4738-9152-57285dd7324f)](https://youtu.be/V5tDbxtWsfY)  
  
Charpy-test.inp, Charpy-test-glmsh.inp  
[![''](https://github.com/user-attachments/assets/34fed318-ca60-479c-9620-71e320c64c21)](https://youtu.be/05grq5s9JtM)  
  
Tensile-test.inp, Tensile-test-glmsh.inp (Dog-bone)  
[![''](https://github.com/user-attachments/assets/430d8c8f-69ab-41f7-9ad8-02c56e010f8b)](https://youtu.be/zrBSUzXMs58)  
  
bullet-impact.inp, bullet-impact-glmsh.inp, 
[![''](https://github.com/user-attachments/assets/46943222-3946-465b-9b9e-a4c64aeb565b)](https://youtu.be/RczUhYBTxVg)  
  
Tensile5e.inp, Tensile5e-glmsh.inp  
[![''](https://github.com/user-attachments/assets/0954b6b4-574d-4911-a83a-97b43eabf731)](https://youtu.be/7pOkPujlj1I)
  
