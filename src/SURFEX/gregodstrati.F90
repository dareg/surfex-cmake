SUBROUTINE GREGODSTRATI(KA,KM,KJ,PH,PAGE)             
! --------------------------------------------------------------------------    ! inspire du programme gregod de crocus lui meme inspire de GMAP

! ==> il faudra revenir au gregod d'origine
! **** *GREGOD * - Conversion Date > Ecart par rapport a une date fixe.
! --------------------------------------------------------------------------             
! Auteur:        
! -------        
! 92-05-27, J.M. Piriou.      
!         
! Modifications:        
! --------------        
!         
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
! --------------------------------------------------------------------------        
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres    
! 2 : a l'heure pres    
! 3 : a la minute pres   
! 4 : a la seconde pres   
! - si kopt=1 : kdat au format AAAAMMQQ 
! - si kopt=2 : kdat au format AAAAMMQQHH 
! - si kopt=3 : kdat au format AAAAMMQQHHMM 
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS    
! En sortie:          
! - si kopt=1 : kgre nombre de jours entre 19000101 et kdat
! - si kopt=2 : kgre nombre d'heures entre 1900010100 et kdat 
! - si kopt=3 : kgre nombre de minutes entre 1.90001E+11 et kdat
! - si kopt=4 : kgre nombre de secondes entre 1.90001E+13 et kdat
! --------------------------------------------------------------------------           
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++           
! ATTENTION A LA PRECISION:        
! 1 Vous compilez les entiers sur 32 bits:    
! Vous devez alors vous limiter a kopt <= 2   
! 2 Vous compilez les entiers sur 64 bits:    
! Vous pouvez utiliser toutes les valeurs de kopt.    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++           


!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE 
INTEGER , INTENT(IN) :: KA,KM,KJ
REAL, INTENT(IN) ::PH
REAL, INTENT(OUT) ::PAGE
INTEGER(KIND=4) :: idebm(12)
INTEGER(KIND=4) :: I0
INTEGER(KIND=4) :: IA100
INTEGER(KIND=4) :: IA4
INTEGER(KIND=4) :: IA400
INTEGER(KIND=4) :: IAAAA
INTEGER(KIND=4) :: IBISSEXT
INTEGER(KIND=4) :: ICONV
INTEGER(KIND=4) :: IFRJOUR
INTEGER(KIND=4) :: II
INTEGER(KIND=4) :: II1
INTEGER(KIND=4) :: IJOURP
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IN
INTEGER(KIND=4) :: IN1
INTEGER(KIND=4) :: IN2           
INTEGER(KIND=4) :: IQQ           
INTEGER(KIND=4) :: KDAT           
INTEGER(KIND=4) :: KGRE           
INTEGER(KIND=4) :: KOPT           
REAL(KIND=JPRB) :: ZHOOK_HANDLE
data idebm/0,31,59,90,120,151,181,212,243,273,304,334/            
!             
! --------------------------------------------------------------------------            
! ** 1 Calcul du nb de jours separant ki2 du 1er janv 1900
!             
! * 1.1 Extraction des quantieme, mois et annee     
 
! on fait calcul dans le cas de l'option x? ! Date de type AAAAMMQQ
IF (LHOOK) CALL DR_HOOK('GREGODSTRATI',0,ZHOOK_HANDLE)
iconv=1    
ifrjour=0    
iqq=KJ
imm=KM     
iaaaa=KA     
! * 1.2 L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est     
! (mult de 4 et non mult de 100) ou (mult de 400)
ia400=400*(iaaaa/400)            
ia100=100*(iaaaa/100)            
ia4=4*(iaaaa/4)            
if((iaaaa == ia400).or.((iaaaa == ia4).and.(iaaaa /= ia100)))then      
ibissext=1            
else            
ibissext=0            
endif            
if ((ibissext == 1).and.(imm > 2)) then      
ijourp=1            
else            
ijourp=0            
endif            
             
in2=idebm(imm)+ijourp+iqq-1               
! * 1.4 Calcul du nb de jours separant les 1er janvier de ii et 2000
i0=2000               
in2=in2+365*(iaaaa-i0)+int((iaaaa-1)/4)-int((i0-1)/4)               &
  -int((iaaaa-1)/100)+int((i0-1)/100)               &
  +int((iaaaa-1)/400)-int((i0-1)/400)                 
! --------------------------------------------------------------------------              
! ** 2 Calcul du nb de jours separant ii1 du 1er janv 2000  
!               
! * 2.1 Extraction des quantieme, mois et annee       
ii1=20000101            
ii=ii1            
iqq=ii-(ii/100)*100            
in=(ii-iqq)/100            
imm=in-(in/100)*100            
iaaaa=(in-imm)/100            
! * 2.2 L'annee est-elle bissextile?       
! Une annee est bissextile ssi elle est     
! (mult de 4 et non mult de 100) ou (mult de 400)
iaaaa=iaaaa            
ia400=400*(iaaaa/400)            
ia100=100*(iaaaa/100)            
ia4=4*(iaaaa/4)            
if((iaaaa == ia400).or.((iaaaa == ia4).and.(iaaaa /= ia100)))then      
ibissext=1            
else            
ibissext=0               
endif               
if ((ibissext == 1).and.(imm > 2)) then         
ijourp=1               
else               
ijourp=0               
endif               
! * 2.3 Nombre de jours ecoules depuis le 1er janv     
in1=idebm(imm)+ijourp+iqq-1               
! * 2.4 Calcul du nb de jours separant les 1er janvier de ii et 2000
i0=2000               
in1=in1+365*(iaaaa-i0)+int((iaaaa-1)/4)-int((i0-1)/4)               &
  -int((iaaaa-1)/100)+int((i0-1)/100)               &
  +int((iaaaa-1)/400)-int((i0-1)/400)                 
! --------------------------------------------------------------------------              
! ** 3 Difference in2-in1 sousforne nbre jours / quantieme de jour           
PAGE=(in2-in1)*iconv+ifrjour+PH/24./3600.              
IF (LHOOK) CALL DR_HOOK('GREGODSTRATI',1,ZHOOK_HANDLE)
END SUBROUTINE GREGODSTRATI
