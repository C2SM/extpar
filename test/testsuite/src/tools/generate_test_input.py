#!/usr/bin/env python

"""
COSMO TECHNICAL TESTSUITE

General purpose script to generate new input from INPUT_XXX in current directory
"""

# built-in modules
import os, sys

# information
__author__     = "Xavier Lapillonne, Nicolo Lardelli"
__email__      = "cosmo-wg6@cosmo.org"
__maintainer__ = "xavier.lapillonne@meteoswiss.ch"


def main():

    l_files=['INPUT_ASS','INPUT_DIA','INPUT_DYN','INPUT_INI',\
                 'INPUT_IO','INPUT_ORG','INPUT_PHY','INPUT_SAT']

    #tests specifications one value per test
    #all lists should have same length
    info=['Only dynamic','Radiation (lphys+lrad)',\
          'Trubulence (lphys+lrad+ltur)','Convection (lphys+lrad+ltur+lconv)',\
          'Soil (lphys+lrad+lconv+ltur+lsoil)',\
          'Sso (lphys+lrad+lconv+ltur+lsoil+lsso)', \
          'Precipitation (lphys+lrad+lconv+ltur+lsoil+lsso+lgsp)',\
          'Observations and all physics']
    
          
    lphys=  ['.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.']
    lrad =  ['.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.']
    ltur =  ['.FALSE.','.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.']
    lconv=  ['.FALSE.','.FALSE.','.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.']
    lsoil = ['.FALSE.','.FALSE.','.FALSE.','.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.' ,'.TRUE.']
    lsso =  ['.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.TRUE.' ,'.TRUE.' ,'.TRUE.']
    lgsp =  ['.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.TRUE.' ,'.TRUE.']
    luseobs=['.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.FALSE.','.TRUE.']
    
    

    syso=os.system('rm -rf new_INPUT')
    syso=os.system('mkdir new_INPUT')

    for i in range(len(lphys)):
        print('generating test '+str(i+1))
        for fname in l_files:
            new_fname='new_INPUT/'+fname+'_'+str(i+1)
            syso=os.system('cp '+fname+' '+new_fname)
            if fname=='INPUT_ORG':
                change_line(new_fname,' lphys =','  lphys = '+lphys[i]+',\n')
                change_line(new_fname,' luseobs =','  luseobs = '+luseobs[i]+',\n')
                change_line(new_fname,'$info','! $info = '+info[i]+'\n')
            elif fname=='INPUT_PHY':
                change_line(new_fname,' lgsp =','  lgsp = '+lgsp[i]+',\n')
                change_line(new_fname,' lrad =','  lrad = '+lrad[i]+',\n')
                change_line(new_fname,' lconv =','  lconv = '+lconv[i]+',\n')
                change_line(new_fname,' ltur =','  ltur = '+ltur[i]+',\n')
                change_line(new_fname,' lsoil =','  lsoil = '+lsoil[i]+',\n')
                
                
    print('files written in ./new_INPUT')

#------------------------------------------------------------------
# function change_line(filename,identifier,newline)

def change_line(filename,identifier,newline):
    # change_line(filename,identifier,newline)
    # change line containing identifier in filename with newline

    data=open(filename).readlines()
    new_data=[]
    counter=0
    for line in data:
        if identifier.strip() in line:
            line=newline
            counter+=1
        
        new_data.append(line)

    #check if exaclty one change
    if counter != 1:
        print('Error: %s appears %i times in file %s' %(identifier,counter,filename))
        print('Note: check that white space are not missing in corresponding line %s' %(identifier))
        sys.exit()
        
    fout=open(filename,'w')
    fout.write(''.join(new_data))
    fout.close
    
#------------------------------------------------------------------
# main
if __name__ == "__main__":
      main()


