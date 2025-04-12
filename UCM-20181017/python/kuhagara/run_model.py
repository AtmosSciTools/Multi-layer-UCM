##
import os, sys




text = []
os.system("cp ../../namelist/namelist.ucm.kugahara ../../namelist.ucm")
os.system("cp ../../namelist/URBPARM.TBL.kugahara ../../URBPARM.TBL")


os.system("cd ../../; sh compile.sh; cd -")

sys.exit()
s = open("../../namelist.ucm","r")
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
    
    

    
#os.system("cd ../; sh compile.sh; cd -")
    










     
 