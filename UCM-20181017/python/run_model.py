##
import os, sys




for uc in  [1,2,3,4][:]:
    
    if uc == 4:
        slab = 3
        #----------
        text = []
        s = open("../namelist.ucm","r")
        for line in s.readlines():
            if "ofiname" in line:
                line = ' ofiname = "./slab'+str(slab)+'_"\n'
            text.append(line)
        s.close() 
           
        s = open("../namelist.ucm","w")
        s.write("".join(text))
        s.close() 
        os.system('sed -i "/ sf_surface_physics*/c\ sf_surface_physics  = 1" ../namelist.ucm')
        os.system('sed -i "/ sltype*/c\ sltype = '+str(slab)+'" ../namelist.ucm')
        #---------
        
        
    else:
        
        #----------
        text = []
        s = open("../namelist.ucm","r")
        for line in s.readlines():
            if "ofiname" in line:
                line = ' ofiname = "./case'+str(uc)+'_"\n'
            text.append(line)
        s.close()    
        s = open("../namelist.ucm","w")
        s.write("".join(text))
        s.close() 
        os.system('sed -i "/utype*/c\ utype = '+str(uc)+'" ../namelist.ucm')
        os.system('sed -i "/ sltype*/c\ sltype     = 3" ../namelist.ucm')
        os.system('sed -i "/ sf_surface_physics*/c\ sf_surface_physics  = 2" ../namelist.ucm')
        #----------
    
    

    
    
    os.system("cd ../; sh compile.sh; cd -")
    










     
 