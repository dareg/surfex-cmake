# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:28:18 2013

@author: seb
"""

import time
import os
import copy
import numpy as np
import myfct2
import re
#import matplotlib.pyplot as plt
import netCDF4
from Scientific.IO.NetCDF import NetCDFFile


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

def get_lastdir(path) :

    return (os.path.split(os.path.dirname(path))[1])
    
#---------------------------------------------------------------------------

def get_nc_var_list(path) :

     if os.path.isfile(path)==True :
   
        f = NetCDFFile(path)
        var_list = f.variables.keys() 
        f.close()

        err = False

        return var_list,err
        
     else :

        print("file "+str(path)+" doesn't exist")

        err = True 

        return '', err

#---------------------------------------------------------------------------

def chg22Dlist(list1,fic) :
    
    for i in range(len(list1)) : 

        list1[i] = ([list1[i],fic])

    return list1

#---------------------------------------------------------------------------
            
def extract_common_files(list1,list2) :
    
    listc  = []
    list1a = []
    list2a = copy.deepcopy(list2)
    
    for f in list1 :

        try :

            ind2 = list2.index(f)

        except : 

            list1a.append(f)

        else :

            listc.append(list2[ind2])   
            list2a.remove(list2[ind2])
    
    return listc, list1a, list2a

#---------------------------------------------------------------------------

def Save_list2txt (path,head,array) :

    k = len(str(len(array)))

    li_data = [ str(i+1).rjust(k)+' '+str(x)+'\n' for i,x in enumerate(array) ]

    open(path,'w').write(head+'\n'+''.join(li_data))

#---------------------------------------------------------------------------   
  
def get_nc_agreg_list(dir) :
    
    agregvar=[]

    files_list = myfct2.get_nc_out_list(dir)

    for fi in files_list :

        varl,err = get_nc_var_list(dir+fi)
        b = []

        for  var in varl :

            if var not in ['xx','yy','time'] :

                b.append([var,fi])

        agregvar.extend(b)
    
    return agregvar    


def get_nc_out_list(dir) :

    listnc = [ s for s in os.listdir(dir) if s.endswith('.OUT.nc') ]

    return listnc   

#---------------------------------------------------------------------------     

def trait_list(step,listM,dirN,dirO,simul,pathrpt,patherr,pathsynn,pathsyno) :
  
    listEc  = []  
    listall = []  
    listnan = []
    
    nbnan = 0
    nberr = 0
    
    if len(listM)<>0 :
  
        for i in range(len(listM)-0) :
            
            """ get namefile correponding to var """  
            fic1 = dirN+listM[i][1]    
            fic2 = dirO+listM[i][1]
            var  = listM[i][0]
	    nanN = False
            nanO = False
            """ return if presence of nan ,list of values, table of differences Emax index ...  """
	    obj,nanN,nanO,tableN,tableO,tableD,Emax,indmax,Emin,indmin,DMoy = myfct2.diff_rpt(fic1,fic2,var)
   
	    if obj[0]=='0' :

		Err = np.nan
		indmax = np.nan
                Cmpnewval = tableN
                Cmpoldval = tableO
		
	    elif obj=='no' :

		Err = 0.0
		Emax = 0.0
		indmax = 0
                Cmpnewval = tableN
                Cmpoldval = tableO

	    else :
	 	
                Cmpnewval = tableN[indmax]
                Cmpoldval = tableO[indmax]

            """ create list of relative distance only for correct array """
            if obj=='XDr' or obj=='XDi'or obj=='0Dr' :

                Err = 0 if (Emax==0.0 and Cmpoldval==0.0) else 100*abs(2*Emax/(Cmpoldval+Cmpnewval))


	    if obj=='XDs' or obj=='XDb' : 
                
                Err = 100*DMoy

 	    if obj=='0Ds' or obj=='0Db' or obj=='0Di' :


		Err = 100*Emax
		DMoy = Err
		indmax = 0


            listEc.append([Err,obj,indmax,Emax,Cmpnewval,Cmpoldval,DMoy,listM[i]]) 


            if (nanN or nanO)==True :

                nbnan = nbnan+1
                
		""" create list files containing nan values """    
                listnan.append([nanN,nanO,Emax-Emin,listM[i]])
            
	    if (Emax<>0. and Emax<>-0.) :
                nberr = nberr+1
            
            """ add field name   """  
     
        """sort list by 1st value  """              
        listEc_s = myfct2.org_sort(listEc)

        #on informe les fichiers de synthese du deroulement de cette epx
        if nbnan>0 :

            open(pathsynn,'a').write(simul+": "+step+" Nan: " + str(nbnan) + "\n")

        if nberr>0 :

	    open(pathsynn,'a').write(simul+": "+step+" Emax<>0: " + str(nberr) \
			    + " " + str(listEc_s[0][0]) + " " + str(listEc_s[0][3])  \
			    + " " + str(listEc_s[1][0]) + " " + str(listEc_s[1][3]) + "\n")

        if nbnan==0 and nberr==0 :

            open(pathsyno,'a').write(simul+": "+step+"\n")
            
            
        
        """ save list of relative error  and nan values"""
        headrpt = '%Err ; NC ; indmax ; DeltaMax ; NewValue ; OldValue ; Mean ; File '
        myfct2.Save_list2txt(pathrpt,headrpt,listEc_s)
        
        headnan='nanpresenceinNew ; nanpresenceinOld; DEltaMax; File'
        myfct2.Save_list2txt(patherr,headnan,listnan)
        
        """ attention too many files open : close -> close()"""        
            
            
#-------utilisees dans myfct2---------------------

def diff_rpt(path1,path2,field) :  

    """ only 2D so far  """
    
    DataN = myfct2.get_ncdf_1data(path1,field)
    DataO = myfct2.get_ncdf_1data(path2,field)

    indmax=np.nan; indmin=np.nan; Moy=np.nan; Emin=np.nan; tableD=np.nan
    nanN=np.nan; nanO=np.nan; Emax=np.nan;

    obj = ''
    if type(DataN)==type(np.array([])) :

	if DataN.dtype.type==np.float_ :
	    if len(DataN)==len(DataO) :
               obj = 'XDr'
	       DataN = np.reshape(DataN,DataO.shape)
               nanN,nanO,tableD,Emax,indmax,Emin,indmin,Moy = diff_XDr(DataN,DataO)
	    else :
               obj = 'no'
               #print(field,'no comparison because 19 vegtypes',len(DataN),len(DataO))

	if DataN.dtype.type==np.float32 :
	    if len(DataN)==len(DataO) :
               obj = 'XDr'
	       DataN = np.reshape(DataN,DataO.shape)
               nanN,nanO,tableD,Emax,indmax,Emin,indmin,Moy = diff_XDr(DataN,DataO)
	    else :
               obj = 'no'
               #print(field,'no comparison because 19 vegtypes',len(DataN),len(DataO))

	if DataN.dtype.type==np.int32 :

	    if len(DataN)==len(DataO) :
               obj = 'XDi'
               nanN,nanO,tableD,Emax,indmax,Emin,indmin,Moy = diff_XDi(DataN,DataO)
	    else :
               obj = 'no'
               #print(field,'no comparison because 19 vegtypes',len(DataN),len(DataO))

        if DataN.dtype.type==np.string_ :

	    if len(DataN)==len(DataO) :
               obj = 'XDs'
               nanN,nanO,Emax,indmax,Moy = myfct2.diff_XDs(DataN,DataO)  
	    else :
               obj = 'no'
               #print(field,'no comparison because 19 vegtypes',len(DataN),len(DataO))
        
        if DataN.dtype.type==np.bool_ :

	    if len(DataN)==len(DataO) :
               obj = 'XDb'
               Emax,indmax,Moy = myfct2.diff_XDb(DataN,DataO)
	       nanN = False; nanO = False
	    else :
               obj = 'no'
               #print(field,'no comparison because 19 vegtypes',len(DataN),len(DataO))
            
    elif type(DataN)==type(DataO) :
        
        """ if scalar """
        if DataN.dtype.type==np.float_ :

            obj = '0Dr'
            nanN,nanO,Emax = myfct2.diff_0Dr(DataN,DataO) 

        if DataN.dtype.type==np.int32 :

            obj = '0Di'
            nanN,nanO,Emax = myfct2.diff_0Dr(DataN,DataO) 

        """ if string """
	if DataN.dtype.type==np.string_  :

            obj = '0Ds'
            nanN,nanO,Emax = myfct2.diff_0Ds(DataN,DataO) 
            
        """ if boolean """
        if DataN.dtype.type==np.bool_ :

            obj='0Db'
            Emax = myfct2.diff_0Db(DataN,DataO)
	    nanN = False; nanO = False
	
    else :
        obj = 'no'
        print(field,'changed type')
           
    if obj=='' : 

        print(field,'not found')
    

    return obj,nanN,nanO,DataN,DataO,tableD,Emax,indmax,Emin,indmin,Moy


def strA2str(strA) :

    str = ''

    for char in strA :

        str=str+char
    
    return str


def org_sort(L,col=0) :

    par = [ item[col] for i,item in enumerate(L)  ]     
    nanl = np.where(np.isnan(par))[0]
    L2 = sublist(L,nanl)
    
    L3 = [ item for i,item in enumerate(L) if i not in nanl ]
    L3 = sorted(L3,key=lambda x:(x[0]),reverse=True)
    L3.extend(L2)
    
    return L3


def get_ncdf_1data(path,var):
   
        f = NetCDFFile(path)
        valnc= f.variables[var] 
        val=valnc.getValue()
        f.close()
        return val


def sublist(L, inds) :

    return [L[i] for i in inds]   
    

def diff_0Db(resn,reso) :

    if (resn==reso) :

        E = 0

    else :

        E = 1

    return E


def diff_XDb(resn,reso) :

    C = 0
    E = 0
    indmax = 0

    for i in range(len(resn)-0) :

	C = C+1

    	if (not(resn[i]==reso[i])) :
	
          E = E+1
	  indmax = i

    mean = E/C

    return E,indmax,mean


def diff_0Ds(strn,stro) :

    E = diff_0Db(strn,stro)
        
    nanN = len(strn)<1
    nanO = len(stro)<1
    
    return nanN,nanO,E


def diff_XDs(strn,stro) :

    E,indmax,mean = diff_XDb(strn,stro)
    
    nanN = ((strn==' ').all()) or ((strn=='').all())
    nanO = ((stro==' ').all()) or ((stro=='').all())
    
    return nanN,nanO,E,indmax,mean


def diff_0Dr(resn,reso) :
    
    nanN = np.isnan(resn)
    nanO = np.isnan(reso)
	
    D = (resn-reso)

    return nanN,nanO,D


def diff_XD(resn,reso) :
        
    """ test si au moins 1 nan """    
    nanN = np.isnan(resn).any()
    nanO = np.isnan(reso).any()
            
    """ Table D without nan values """    
    tableD = abs(np.subtract(resn,reso))
    
    """ mean calculus """
    inan = ~np.isnan(tableD)
    tableDnan = tableD[inan]
    mean = tableDnan.mean()
    
    """ ecart max -> Erel""" 
    indmax = np.unravel_index(np.nanargmax(tableD), tableD.shape)
    Emax = tableD[indmax]

    indmin = np.unravel_index(np.nanargmin(tableD), tableD.shape)
    Emin = tableD[indmin]
        
    return nanN,nanO,tableD,Emax,indmax,Emin,indmin,mean


def diff_XDr(resn,reso) :
            
    """ shift 1e+20 to 0 """ 
    np.putmask(resn, resn >= 1.e+20, 0.)
    np.putmask(reso, reso >= 1.e+20, 0.)

    nanN,nanO,tableD,Emax,indmax,Emin,indmin,mean = diff_XD(resn,reso)

    return nanN,nanO,tableD,Emax,indmax,Emin,indmin,mean


def diff_XDi(resn,reso) :
            
    """ shift 1e+9 to 0 """        
    np.putmask(resn, resn >= 1e+9, 0)
    np.putmask(reso, reso >= 1e+9, 0)
       
    nanN,nanO,tableD,Emax,indmax,Emin,indmin,mean = diff_XD(resn,reso)
        
    return nanN,nanO,tableD,Emax,indmax,Emin,indmin,mean

