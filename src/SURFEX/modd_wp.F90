MODULE modd_wp
  !! necessary for compiling Gelato, but not included in LIB/GELATO for reasons relevant to 
  !! Gelato build process
   INTEGER, PUBLIC, PARAMETER ::   wp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
END MODULE modd_wp
