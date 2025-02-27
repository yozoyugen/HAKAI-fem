## HAKAI is ...  
Simple finite element method dynamic explicit solver for elastoplastic and fracture analysis.  
There are Matlab(.m) and Julia(.jl) versions  

<a href="https://www.youtube.com/watch?v=b7lrxCQCVUU"><img src="https://github.com/user-attachments/assets/b1b8b0d1-12a9-4463-adc6-955fad97c07d" width="400px"></a>
<a href="https://www.youtube.com/watch?v=CabZ8rJxqiU"><img src="https://github.com/user-attachments/assets/6dce3420-cc3c-44d0-bf8b-7c76bc818d1c" width="400px"></a>  
(Click images to show YouTube videos)  
  
## Web site  
[https://caelab.work/](https://caelab.work/)
  
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
Click images to show YouTube videos  
<a href="https://youtu.be/V5tDbxtWsfY"><img src="https://github.com/user-attachments/assets/074c40b4-2ee5-4738-9152-57285dd7324f" width="400px"></a>
<a href="https://youtu.be/05grq5s9JtM"><img src="https://github.com/user-attachments/assets/34fed318-ca60-479c-9620-71e320c64c21" width="400px"></a>  

<a href="https://youtu.be/zrBSUzXMs58"><img src="https://github.com/user-attachments/assets/430d8c8f-69ab-41f7-9ad8-02c56e010f8b" width="400px"></a>
<a href="https://youtu.be/RczUhYBTxVg"><img src="https://github.com/user-attachments/assets/46943222-3946-465b-9b9e-a4c64aeb565b" width="400px"></a>  

<a href="https://youtu.be/7pOkPujlj1I"><img src="https://github.com/user-attachments/assets/0954b6b4-574d-4911-a83a-97b43eabf731" width="400px"></a>

