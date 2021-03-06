<a href="https://imt-nord-europe.fr/en/"><img src="https://github.com/adeboracris/Epanet_Leak_localization/blob/main/IMG/1200px-%C3%89cole_des_Mines_de_Douai.svg.png" width="150">  <a href="https://www.upc.edu/ca"><img src="https://github.com/adeboracris/Epanet_Leak_localization/blob/main/IMG/upc.png" width="300"><a>

# <img src="https://media.giphy.com/media/VgCDAzcKvsR6OM0uWg/giphy.gif" width="50"> Epanet Leak localization  

  ## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [How to use](#How-to-use)
- [Contributors](#Contributors)
  
  
  ## Introduction
    This project refers to the leak localization, and sensor validation method explained in the 'Robust Data-Driven Leak Localization In Water DistributionNetworks Using Pressure Measurements And TopologicalInformation' article.  
  
  This is the research carried out within Universidad Politécnica de Cat
  aluña (UPC) and École nationale supérieure des Mines-Télécom de Lille-Douai (IMT). The PHD candidate is Débora Alves with the theme Leak Supervision in Water Distribution Networks.

 
  
  ## Requirements

* [Matlab](http://www.mathworks.com/)

  ## How to use
  
  The folder <b>Leak_localization</b> contain: 3 data .mat coresponding the Hanoi network data (graph, topologie, measurements) and 2 folder .m that are the base code.
  The code that generates the plot is: main.m. The code will generate the following plot:
  
  <img src="https://github.com/adeboracris/Epanet_Leak_localization/blob/main/IMG/atd.jpg" width="500"> 
  
   The folder <b>Sensor validation</b> contain: a .mat file with some Hanoi measurements information and the main code (Spacial_residual.m) that is used to generate the following plots:
  
  <img src="https://github.com/adeboracris/Epanet_Leak_localization/blob/main/IMG/filtered.jpg" width="400"> 
  
  <img src="https://github.com/adeboracris/Epanet_Leak_localization/blob/main/IMG/Spacial.jpg" width="400"> 

  ## Contributors
  
    
  [Joaquim Blesa](https://www.iri.upc.edu/staff/jblesa)
    
  [Eric Duviella](https://sites.google.com/site/ericduviella/)
    
  [Lala Rajaoarisoa](https://www.researchgate.net/profile/Lala-Rajaoarisoa)
  
 
  
  ==================================
  
  &uparrow; [Back to top](#table-of-contents)
