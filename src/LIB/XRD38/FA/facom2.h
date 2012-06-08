      INTEGER I1234, I5678, I3456, IXNVMS
      LOGICAL LLMOER
C
C       Fonction servant a rendre fatale ou non une erreur detectee,
C       a l'aide du code reponse courant, du niveau de filtrage global,
C       et de l'option d'erreur fatale propre au fichier.
C       s'il n'y a pas de fichier (I5678=0, d'ou dimensionnement de
C          *LERRFA*), le niveau de filtrage joue le role principal.
C
      LLMOER (I1234,I5678)=I1234.EQ.-66.OR.
     S (I1234.NE.0.AND.(FA%NRFAGA.EQ.0.OR.
     S (FA%NRFAGA.EQ.1.AND.FA%LERRFA(I5678))))
C*
C       Fonction "en ligne" donnant le plus haut niveau de messagerie
C       acceptable pour l'unite logique de rang "I3456"
C       (utilisation des niveaux de messagerie global et propre au
C        fichier; s'il n'y a pas de fichier - I3456=0, d'ou le dimensio-
C        nnement de *NIVOMS* a partir de zero, le niveau de filtrage
C        global joue seul)
C
      IXNVMS (I3456)=MIN (2,2*FA%NIMSGA,MAX 
     S                  (2*FA%NIMSGA-2,FA%NIVOMS(I3456)))
C
