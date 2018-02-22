#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 14:58:44 2013
@author: dasprezs
"""

import sys
import myfct2
import os
import numpy as np
import copy
import re
import subprocess
import datetime
import ConfigParser


def Vercomp(dirN='a',dirO='b',simul='c') :


    """ Out files 1) differences list 2) list of nan values  """
  
    config = ConfigParser.RawConfigParser() 

    '''    compdebug.cfg   ///     compR.cfg    '''
    config.read('compR.cfg')      

    verN = config.get('name','verN')
    verO = config.get('name','verO')
    basedir = config.get('directory', 'basedir')
    Pgdfile  = config.get('file','Pgdfile')
    Prepfile = config.get('file','Prepfile')
    Outfile = config.get('file','Outfile')

    dirN = basedir+dirN
    dirO = basedir+dirO
        
    resdir = basedir+'TESTS/PYTHON/results/'

    #fichiers de synthese
    pathnch = resdir+verN+'_and_'+verO+'_conf.nch'

    pathsyno = resdir+verN+'_and_'+verO+'_conf.syno'
    pathsynn = resdir+verN+'_and_'+verO+'_conf.synn'

    pathrpt = resdir+verN+'_and_'+verO+'_conf_'+simul+'.rpt'
    patherr = resdir+verN+'_and_'+verO+'_conf_'+simul+'.nan'
    pathlog = resdir+verN+'_and_'+verO+'_conf_'+simul+'.log'
    
    pathrptPgd = resdir+'PGD_'+verN+'_and_'+verO+'_conf_'+simul+'.rpt'
    patherrPgd = resdir+'PGD_'+verN+'_and_'+verO+'_conf_'+simul+'.nan'
    pathlogPgd = resdir+'PGD_'+verN+'_and_'+verO+'_conf_'+simul+'.log'       
    
    pathrptPrep = resdir+'PREP_'+verN+'_and_'+verO+'_conf_'+simul+'.rpt'
    patherrPrep = resdir+'PREP_'+verN+'_and_'+verO+'_conf_'+simul+'.nan'
    pathlogPrep = resdir+'PREP_'+verN+'_and_'+verO+'_conf_'+simul+'.log'      

    pathrptOut = resdir+'OUT_'+verN+'_and_'+verO+'_conf_'+simul+'.rpt'
    patherrOut = resdir+'OUT_'+verN+'_and_'+verO+'_conf_'+simul+'.nan'
    pathlogOut = resdir+'OUT_'+verN+'_and_'+verO+'_conf_'+simul+'.log'

    pathfiles = resdir+'mainlog.txt'
   
    head_log = '[listN, listO, listM, listnMN, listnMO]'

    """ -------------create list of all outing variables ------""" 
        
    listPgdN,err1 = myfct2.get_nc_var_list(dirN+Pgdfile)
    listPgdO,err2 = myfct2.get_nc_var_list(dirO+Pgdfile)
            
    q1 = False
    q2 = False
    q4 = False

    if err1 or err2 :

        print("Pgd.nc not found ,exiting code")

        if err1 :
            log = str(dirN+Pgdfile)+' not found, conf '+simul,str(datetime.datetime.now())+'\n'

        if err2 :
            log = str(dirO+Pgdfile)+' not found, conf '+simul,str(datetime.datetime.now())+'\n'

        q1 = True
 
    else :

        log = str(Pgdfile)+'  found, conf '+simul+str(datetime.datetime.now())+'\n'

    	open(pathfiles,"a").write(log)

    #if q : sys.exit(0)

    	""" change to 2D list"""       
   	listPgdN = myfct2.chg22Dlist(listPgdN,Pgdfile)
    	listPgdO = myfct2.chg22Dlist(listPgdO,Pgdfile)
    
    	""" compare two lists, create list of matching, orphan 1st and 2nd list """
    	listPgdM,listPgdnMN,listPgdnMO = myfct2.extract_common_files(listPgdN,listPgdO)

    	""" create log file """
    	List_log = [len(listPgdN),len(listPgdO),len(listPgdM),len(listPgdnMN),len(listPgdnMO)]
    	open(pathlogPgd,'w').write(head_log+'\n'+''.join(str(List_log)))

    """------------------------------------------------------------"""   
    
    listPrepN,err1 = myfct2.get_nc_var_list(dirN+Prepfile)
    listPrepO,err2 = myfct2.get_nc_var_list(dirO+Prepfile)

    if err1 or err2 :

        print("PREP.nc not found ,exiting code")

        if err1 :
            log = str(dirN+Prepfile)+' not found, conf '+simul+str(datetime.datetime.now())+'\n'

        if err2 :
            log = str(dirO+Prepfile)+' not found, conf '+simul+str(datetime.datetime.now())+'\n'          
        
        q2 = True

    else :

        log = str(Prepfile)+' found, conf '+simul+str(datetime.datetime.now())+'\n'
       
    	open(pathfiles,"a").write(log)

    	#if q : sys.exit(0)
    
    	""" change to 2D list"""
    	listPrepN = myfct2.chg22Dlist(listPrepN,Prepfile)
    	listPrepO = myfct2.chg22Dlist(listPrepO,Prepfile)  
    
    	""" compare two lists, create list of matching, orphan 1st and 2nd list """
    	listPrepM,listPrepnMN,listPrepnMO = myfct2.extract_common_files(listPrepN,listPrepO)  
    
    	""" create log file """    
    	List_log = [len(listPrepN),len(listPrepO),len(listPrepM),len(listPrepnMN),len(listPrepnMO)]
    	open(pathlogPrep,'w').write(head_log+'\n'+''.join(str(List_log)))

    """-------------------------------------------------------"""  

    """ create list of all outing variables """ 
    #liste des variables dans les sorties "NEW"
    listN  =myfct2.get_nc_agreg_list(dirN)
    #liste des variables dans les sorties "OLD"
    listO = myfct2.get_nc_agreg_list(dirO)
  
    q3 = False

    """ Cas d une liste vide """   
    if len(listN)<1 : 

        log = "Data not found, conf "+simul+str(datetime.datetime.now())+'\n'

	if len(listO)>0 : 
		print("No data found ,exiting code")
		q3 = True

    else :

        log = "Data  found, conf "+simul+str(datetime.datetime.now())+'\n'

    open(pathfiles,"a").write(log)

    #if q3 : sys.exit(0)

    """ compare two lists, create list of matching, orphan 1st and 2nd list """
    #listM: liste des champs communs
    #listnMN: presents dans NEW mais pas dans OLD
    #listnMO: presents dans OLD mais pas dans NEW
    listM,listnMN,listnMO = myfct2.extract_common_files(listN,listO)
 
    List_log = [len(listN),len(listO),len(listM),len(listnMN),len(listnMO)]
    open(pathlog,'w').write(head_log+'\n'+''.join(str(List_log)))


    #recapitulatif des champs presents NEW OLD pour cette exp dans fichier de synthese
    log2 = simul+' - listN:'+str(len(listN))+', listO:'+str(len(listO))+', listM: '+str(len(listM))
    log2 = log2+', listnMN:'+str(len(listnMN))+', listnMO:'+str(len(listnMO))
    open(pathnch,'a').write(log2 + '\n')
   
    """----------------------------------------------------"""
   
    listOutN,err1 = myfct2.get_nc_var_list(dirN+Outfile)
    listOutO,err2 = myfct2.get_nc_var_list(dirO+Outfile)

    if err1 or err2 :

        print("SURFOUT.nc not found ,exiting code")

        if err1 :
            log = str(dirN+Outfile)+' not found, conf '+simul+str(datetime.datetime.now())+'\n'

        if err2 :
            log = str(dirO+Outfile)+' not found, conf '+simul+str(datetime.datetime.now())+'\n'          
        
        q4 = True

    else :

        log = str(Outfile)+' found, conf '+simul+str(datetime.datetime.now())+'\n'
       
    	open(pathfiles,"a").write(log)

    	#if q : sys.exit(0)
    
    	""" change to 2D list"""
    	listOutN = myfct2.chg22Dlist(listOutN,Outfile)
    	listOutO = myfct2.chg22Dlist(listOutO,Outfile)  
    
    	""" compare two lists, create list of matching, orphan 1st and 2nd list """
    	listOutM,listOutnMN,listOutnMO = myfct2.extract_common_files(listOutN,listOutO)  
    
    	""" create log file """    
    	List_log = [len(listOutN),len(listOutO),len(listOutM),len(listOutnMN),len(listOutnMO)]
    	open(pathlogOut,'w').write(head_log+'\n'+''.join(str(List_log)))

    """-------------------------------------------------------"""  

    """ listM   is list of [varname][filename]...  """     
    
    print """reportPgd"""
    if not q1 :

    	myfct2.trait_list('PGD ',listPgdM,dirN,dirO,simul,pathrptPgd,patherrPgd,pathsynn,pathsyno)

    print """reportPrep"""
    if not q2 :

    	myfct2.trait_list('PREP ',listPrepM,dirN,dirO,simul,pathrptPrep,patherrPrep,pathsynn,pathsyno)   
    
    print """reportdata"""
    myfct2.trait_list('OUTPUT ',listM,dirN,dirO,simul,pathrpt,patherr,pathsynn,pathsyno)

    print """reportOut"""
    if not q4 :

    	myfct2.trait_list('SURFOUT',listOutM,dirN,dirO,simul,pathrptOut,patherrOut,pathsynn,pathsyno)  

if __name__ == "__main__":
    Vercomp(sys.argv[1],sys.argv[2],sys.argv[3])

