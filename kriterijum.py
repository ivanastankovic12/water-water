

"""
pravljenje sktripte za odredjivanje vrste vodonicne veze prema Milanovim kriterijumima

python kriterijum.py .cif imena 6 atoma vode
ucitam strukturu
normalizujem veze
odredim ugao vode, mora da bude u opsegu
OO rastojanje mora biti <=4A
odredim kojoj grupi interakcija pripada
stampam 7 parametara koji su vezani za kriterijume i P1/P2
"""


import sys
import os


import math
import numpy as np


##############################################
##############  PROMENI OVO  #################
cif=open(sys.argv[-7], 'r')   
prvavoda1=sys.argv[-6]  # imena atoma sa zarezima izmedju
prvavoda2=sys.argv[-5]  # imena atoma sa zarezima izmedju
prvavoda3=sys.argv[-4]  # imena atoma sa zarezima izmedju
drugavoda1=sys.argv[-3]   # imena atoma sa zarezima izmedju
drugavoda2=sys.argv[-2]   # imena atoma sa zarezima izmedju
drugavoda3=sys.argv[-1]   # imena atoma sa zarezima izmedju

waterdown=96.4
waterup=112.8

## ako je xaxisatom vodonik, onda se njegove koordinate normiraju pre odredjivanja distxint (norm=True), i Hnorm je duzina veze xaxisatom-zeroatom
norm=True
Hnorm= 0.993

####################################
####################################


prvavoda=[prvavoda1,prvavoda2,prvavoda3]
drugavoda=[drugavoda1,drugavoda2,drugavoda3]
print "COKA", prvavoda, drugavoda

def angle3points(x,y,z):
    Vector1 = x - y
    Vector2 = z - y
    a = np.dot(Vector1,Vector2) 
    b = np.linalg.norm(Vector1)
    c = np.linalg.norm(Vector2) 
    d = (a/b)/c
    # ogranicim arc na (-1,1)
    if d > 1:
      d=1
    if d < -1:
      d=-1

    angle = np.arccos(d) 
    angletotal = math.degrees(angle)
    return angletotal 

def dihedral(a,b,c,d):   # arrays
    # ravan (a,b,c)
    (u,D)=planevector(a,b,c)
    # ravan (b,c,d)
    (v,D)=planevector(b,c,d)
    return angle3points(u,[0,0,0],v)

def planevector(x,y,z):  # vektor normale na ravan koju cine 3 tacke
    """This function will return the plane equation parameters a,b,c 
        for each three dots that are in the plane   ax+by+cz+d=0  """
    # defining two vectors that start from the same dot, let's choose dot x
    p_xy = y-x
    p_xz = z-x
    # finding ortogonal vector to the ones we have just created
    norm_p = np.cross(p_xy,p_xz)
    #  ax+by+cz+d = 0  (a,b,c)*(x,y,z)+d=0 for any dot in the plane (x,y,z) let's choose X, norm_p * X = d
    d = -1*np.dot(norm_p,x)
    # direction vector presents a,b and c
    a = norm_p[0]
    b = norm_p[1] 
    c = norm_p[2] 
    n=np.array([a,b,c])
    return n, d # vektor normale, parametar u jednacini ravni






# razdvajanje imena atoma  ## nepotrebno
na=0
for name in prvavoda:
 if name[0]=="O":
  zeroatom=name
 if name[0]=="H":
  na+=1
  if na==1:
   rotatedatom=name
  if na==2:
   xaxisatom=name

na=0
for name in drugavoda:
 if name[0]=="O":
  distatom=name
 if name[0]=="H":
  na+=1
  if na==1:
   H3=name
  if na==2:
   H4=name



### ulazim u cif i odredjujem koordinate na osnovu imena atoma
  
##################
############
### citanje cifa
cell=[]
zeroatomc=[]
distatomc=[]
xaxisatomc=[]
rotatedatomc=[]
H3c=[]
H4c=[]
selinteracting=[]
selalign=[]
parse=False
parsecoord=False
Ms=[]   ## lista matrica transformacija iz jednog cifa

for linecif in cif:

 ### citanje vektora celije
 if linecif.startswith("_cell_length_") or linecif.startswith("_cell_angle_") :

  parse=False
  a=linecif.split()[1]
  if linecif.split()[1][-1]==")":
   n= linecif.split()[1].index('(')
   a= linecif.split()[1][0:n]
  cell.append(float(a))
 if linecif.startswith("_atom_site_fract_z") :
  parsecoord=True 



 ### citanje koordinata atoma, stavljanje u listu  -- nepotrebno
 if parsecoord==True:
  if linecif.startswith(zeroatom+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    zeroatomc.append(float(c))
  if linecif.startswith(distatom+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    distatomc.append(float(c))
  if linecif.startswith(xaxisatom+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    xaxisatomc.append(float(c))
  if linecif.startswith(rotatedatom+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    rotatedatomc.append(float(c))
  if linecif.startswith(H3+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    H3c.append(float(c))
  if linecif.startswith(H4+" ") and "1_555" not in linecif:
   for c in linecif.split()[2:5]:
    if c[-1]==")":
     c= c[0:c.index('(')]
    H4c.append(float(c))


  ## ovde imam sve sto mi treba - 2 vode za imenima atoma i koordinatama 
  ### selekcije atoma su liste recnika sa tipom atoma i koordinatama
  for atomsel in drugavoda:
   if linecif.startswith(atomsel+" ") and "1_555" not in linecif:
    dict = {'type': atomsel, 'coord':[]}
    for c in linecif.split()[2:5]:
     if c[-1]==")":
      c= c[0:c.index('(')]
     dict['coord'].append(float(c))
    selinteracting.append(dict)

  for atomsel in prvavoda:
   if linecif.startswith(atomsel+" ") and "1_555" not in linecif:
    dict = {'type': atomsel, 'coord':[]}
    for c in linecif.split()[2:5]:
     if c[-1]==")":
      c= c[0:c.index('(')]
     dict['coord'].append(float(c))
    selalign.append(dict)






 #### citanje transformacija i transformisanje u petlji
 if linecif.startswith("_symmetry_equiv_pos_as_xyz"):
  parse=True
 if parse == True:
  ### trazim linije koje pocinju brojem, posle linije _symmetry_equiv_pos_as_xyz. one sadrze transformacije za dobijanje celije. kada vise ne pocinju brojem, stanem

   if  linecif.split()[0].isdigit():

    algebraicmatrix= linecif.split()[1].split(',')
    R=[]
    T=[]
    Rx=np.array([0,0,0])
    Ry=np.array([0,0,0])
    Rz=np.array([0,0,0])
    for el in algebraicmatrix:
           C= 1
           ### razdvajam reci na slova. ako ima brojeva to je translacija


           if "x" in el:
            if any(i.isdigit() for i in el):
              # za translaciju -  broj je do jednog mesta posle znaka deljenja
              Tx=el[0:el.index('/')+2]
              Tx=Tx.split('/')
              Tx=float(Tx[0])/float(Tx[1])
              # ostatak je za rotaciju 
              el=el[el.index('/')+2:]
            else:
             Tx=0
           

            ## rotacija
            for ell in el:
             if ell=="-" or ell=="+":
            
              #### vracanje u celiju. .cif fajl ima takve koordinate da su neki atomi van celije vec u asimetricnoj jedinici sam da bi molekul bio ceo. loordinatni pocetak je u cosku celije.  vracam u celiju - koordinatni pocetak je u roglju celije, pa ako je neka koordinata negativna treba da se doda odgovarajuci vektor celije
              if ell=="-":
               C=-1
               Tx=Tx+1    ### vracanje u celiju
              if ell=="+":
               C=+1
             else:
              if ell=="x":
               Rx=Rx+np.array([C,0,0])
              if ell=="y":
               Rx=Rx+np.array([0,C,0])
              if ell=="z":
               Rx=Rx+np.array([0,0,C])

           if "y" in el:
            if any(i.isdigit() for i in el):
              # za translaciju -  broj je do jednog mesta posle znaka deljenja
              Ty=el[0:el.index('/')+2]
              Ty=Ty.split('/')
              Ty=float(Ty[0])/float(Ty[1])
              # ostatak je za rotaciju 
              el=el[el.index('/')+2:]
            else:
             Ty=0
           

            ## rotacija
            for ell in el:
             if ell=="-" or ell=="+":
              #### vracanje u celiju. .cif fajl ima takve koordinate da su neki atomi van celije vec u asimetricnoj jedinici sam da bi molekul bio ceo. loordinatni pocetak je u cosku celije.  vracam u celiju - koordinatni pocetak je u roglju celije, pa ako je neka koordinata negativna treba da se doda odgovarajuci vektor celije
              if ell=="-":
               C=-1
               Ty=Ty+1 ### vracanje u celiju
              if ell=="+":
               C=+1
             else:
              if ell=="x":
               Ry=Ry+np.array([C,0,0])
              if ell=="y":
               Ry=Ry+np.array([0,C,0])
              if ell=="z":
               Ry=Ry+np.array([0,0,C])


           if "z" in el:
            if any(i.isdigit() for i in el):
              # za translaciju -  broj je do jednog mesta posle znaka deljenja
              Tz=el[0:el.index('/')+2]
              Tz=Tz.split('/')
              Tz=float(Tz[0])/float(Tz[1])
              # ostatak je za rotaciju 
              el=el[el.index('/')+2:]
            else:
             Tz=0
           

            ## rotacija
            for ell in el:
             if ell=="-" or ell=="+":
              #### vracanje u celiju. .cif fajl ima takve koordinate da su neki atomi van celije vec u asimetricnoj jedinici sam da bi molekul bio ceo. loordinatni pocetak je u cosku celije.  vracam u celiju - koordinatni pocetak je u roglju celije, pa ako je neka koordinata negativna treba da se doda odgovarajuci vektor celije
              if ell=="-":
               C=-1
               Tz=Tz+1  ### vracanje u celiju
              if ell=="+":
               C=+1
             else:
              if ell=="x":
               Rz=Rz+np.array([C,0,0])
              if ell=="y":
               Rz=Rz+np.array([0,C,0])
              if ell=="z":
               Rz=Rz+np.array([0,0,C])

    R.append(Rx)
    R.append(Ry)
    R.append(Rz)
    T.append(Tx)
    T.append(Ty)
    T.append(Tz)


    Ms.append([np.array(R),np.array(T)])

cif.close()
###################
#################
###################
 
 
### racunanje vektora celije
ax=cell[0]
bx=cell[1] * np.cos(math.radians(cell[5]))
by=cell[1] * np.sin(math.radians(cell[5]))
cx= cell[2] * np.cos(math.radians(cell[4]))
cy= cell[2] * (np.cos(math.radians(cell[3]))  -np.cos(math.radians(cell[4])) * np.cos(math.radians(cell[5])) )/ np.sin(math.radians(cell[5]))  
V1=cell[0] * cell[1] *cell[2]
V21=  1 - (np.cos(math.radians(cell[3])))**2  - (np.cos(math.radians(cell[4])))**2- (np.cos(math.radians(cell[5])))**2 
V22=  2 * np.cos(math.radians(cell[3]))*np.cos(math.radians(cell[4]))*np.cos(math.radians(cell[5]))
V= V1 * np.sqrt(V21 + V22)
cz=V / (cell[0] * cell[1]  * np.sin(math.radians(cell[5])) )
cellvectors=[np.array([ax,0,0]),np.array([bx,by,0]),np.array([cx,cy,cz])]

  


## da li su uglovi u vodi dobri  - prvo prebacim u kartezijanske koordinate, trebaju mi vektori celije
water=angle3points(np.dot(rotatedatomc,cellvectors),np.dot(zeroatomc,cellvectors),np.dot(xaxisatomc,cellvectors))
watersec= angle3points(np.dot(H3c,cellvectors), np.dot(distatomc,cellvectors), np.dot(H4c,cellvectors))

if (water>=float(waterdown) and water<=float(waterup)) and (watersec>=float(waterdown) and watersec<=float(waterup)) :  

   ##############################################
   #############################################
   #### PRESLIKAVANJE  u frakcionim koordinatama
   ### za svako preslikavanje, racunam rastojanje OO, preko kartezijanskih koordinata
   ## trazim onu sliku koja je najbliza originalnoj, bez normalizacije vodonika, samo preko OO rastojanja

   stop=False
   brojac=0
   mindistance=100

   zeroatomCART=np.dot(zeroatomc,cellvectors) # skalarni proizvod sa starom bazom
   for M in Ms:
    distatomtrans=M[0].dot(distatomc)+M[1] # prvo rotacija, pa translacija
    for p in   range(-5,5):
     for q in  range(-5,5):
      for r in range(-5,5):
       vec=np.array([p,q,r])
       distatomtrans2=distatomtrans+vec  # transliram duz vektora celije. .matrica iz .cif ne daje sve asimetricne jedinice u istoj celiji, neke interakcije mogu da nedostaju. zbog toga se translira do 3 puta u svim pravcima
       distatomtransCART=np.dot(distatomtrans2,cellvectors)    # tek onda prevedem u kartezijanske - skalarni proizvod sa starom bazom
       ## rastojanje  OO
       distcalc=np.linalg.norm(distatomtransCART-zeroatomCART)
       if distcalc<=4 and distcalc>0 :


        #if distcalc<mindistance:
        # mindistance=distcalc

        brojac+=1
     
        stop=True
          
        ###  preslikam celu selinteracting za nadjene matrice M i vec, u frakcionim koordinatama, pa prevedem u kartezijanske 
        selinteracting2=[]
        for atom in  selinteracting  :
         dict = {'type': atom['type'], 'coord':[]}
         dict['coord']=M[0].dot(atom['coord'])+M[1] +vec
         dict['coord']=np.dot(dict['coord'],cellvectors)
         selinteracting2.append(dict)


        ###  prevedem selalign u kartezijanske, pazi da se ne prepise u petlji
        selalign2=[] ## pravim novu listu za selekciju atoma
        for atom in selalign:
         dict = {'type': atom['type'], 'coord':np.dot(atom['coord'],cellvectors)} ## pravim novi dictionary
         selalign2.append(dict)
    

         

        """
        ###############################
        ##### normalizacija svih vodonika, u kartezijanskim koordinatama
        if norm==True:
         zeroatomNORM=zeroatomCART+ (xaxisatomCART-zeroatomCART)/np.linalg.norm(xaxisatomCART-zeroatomCART)*Hnorm

         for atom in  selalign2  :
         if atom['type'][0]=="H":
          atom['coord']=atom['coord']  zeroatomCART+ (xaxisatomCART-zeroatomCART)/np.linalg.norm(xaxisatomCART-zeroatomCART)*Hnorm
        else:
         xaxisatomNORM=xaxisatomCART  ## bez normiranja vodonika
        """


 
        #### pisanje .xyz fajla   ## sa normalizovanim vodonicima
        #outputfile = open(outputfolder+"/"+pdb+".OO-"+str(dist)+".pdb", 'w')

        #outputfile = open(str(brojac)+"test.xyz", 'w')
        #outputfile.write("%s%s" %(len(selalign2)+len(selinteracting2), "\n")) 
        #outputfile.write("%s" %("\n")) 

        #### stampanje .xyz u standardnom izlazu
        print   "****************************************************" 

        print   ".xyz file of the interacting water molecules:"
        print   "****************************************************" 
 
        print("%s" %(len(selalign2)+len(selinteracting2))) 
        print("" ) 

        index=1

        for atom in selalign2 :  ## kasnije stavim sel1 i sel2
         if atom['type'][0]=="O":
             kiseonik=atom['coord']
        for atom in selinteracting2 :  ## kasnije stavim sel1 i sel2
         if atom['type'][0]=="O":
             kiseonik2=atom['coord']


        for atom in selalign2 :  ## kasnije stavim sel1 i sel2
         #outputfile.write("%s%7i%5s%4s%2s%4s%12.3f%8.3f%8.3f %s" %("ATOM", index, atom['type'] ,"UNK","X",  0, atom['coord'][0],atom['coord'][1],atom['coord'][2], "  1.00  0.00"+"\n"))   ## stampanje .pdb fajla
         ### normalizacija vodonika
         if atom['type'][0]=="H":
              atom['coord']=kiseonik+(atom['coord']-kiseonik)/np.linalg.norm(atom['coord']-kiseonik)*Hnorm
         #outputfile.write("%3s%17.6f%16.6f%16.6f%s" %( atom['type'] , atom['coord'][0],atom['coord'][1],atom['coord'][2],"\n"))   ## stampanje .xyz fajla
         print("%3s%17.6f%16.6f%16.6f" %( atom['type'] , atom['coord'][0],atom['coord'][1],atom['coord'][2]))   ## stampanje .xyz fajla

         index +=1


        for atom in selinteracting2 :  ## kasnije stavim sel1 i sel2
         ### normalizacija vodonika
         if atom['type'][0]=="H":
             atom['coord']=kiseonik2+(atom['coord']-kiseonik2)/np.linalg.norm(atom['coord']-kiseonik2)*Hnorm
         ################
         print("%3s%17.6f%16.6f%16.6f" %( atom['type'] , atom['coord'][0],atom['coord'][1],atom['coord'][2]))   ## stampanje .xyz fajla
         #outputfile.write("%3s%17.6f%16.6f%16.6f%s" %( atom['type'] , atom['coord'][0],atom['coord'][1],atom['coord'][2],"\n"))   ## stampanje .xyz fajla
         #outputfile.write("%s%7i%5s%4s%2s%4s%12.3f%8.3f%8.3f %s" %("ATOM", index, atom['type'] ,"UNK","X",  0, atom['coord'][0],atom['coord'][1],atom['coord'][2], "  1.00  0.00")+"\n")   ## stampanje .pdb fajla
         index +=1
        #outputfile.close()
        print   "****************************************************" 
        print   "" 
        print   "" 

        ###############################
        #############################




  ### normalizovano kroz pisanje .xyz fajla


        ##################### obelezavanje atoma prema Milanovim kriterijumima
        ### preslikavanje   dict['coord']=M[0].dot(atom['coord'])+M[1] +vec
        ####  prebacivanje u kartezijanske    dict['coord']=np.dot(dict['coord'],cellvectors)
        ### za sada imam   distatomtransCART i zeroatomCART



        
        ### obelezavanje Ha1 i Ob - oni sa najkracim OH rastojanjem
        OH=100
        for atomalign in selalign2:
         atomalign.update( {'label' : []} )  #  ### dodavanje obelezavanja u selekcije
         for atominteracting in selinteracting2:
          atominteracting.update( {'label' : []} )  #  ### dodavanje obelezavanja u selekcije
          if (atomalign['type'][0]=="H" and  atominteracting['type'][0]=="O") or (atomalign['type'][0]=="O" and  atominteracting['type'][0]=="H"):
           if np.linalg.norm(atomalign['coord']-atominteracting['coord'])<OH:
            OH=np.linalg.norm(atomalign['coord']-atominteracting['coord'])
            if atomalign['type'][0]=="H":
             Ha1=atomalign
             MOLa=selalign2  #### obelezim molekul a i molekul b
            else:
             Ob=atomalign
             MOLb=selalign2
            if atominteracting['type'][0]=="H":
             Ha1=atominteracting
             MOLa=selinteracting2
            else:
             Ob=atominteracting
             MOLb=selinteracting2



        Ha1.update( {'label' : "Ha1"} )  #  ### dodavanje obelezavanja u selekcije
        Ob.update( {'label' : "Ob"} )  #  ### dodavanje obelezavanja u selekcije


       
        ##### odredjujem ostale atome u molekulu a
        for atom in MOLa:
         if atom['type'][0]=="O":
          Oa=atom
          Oa.update( {'label' : "Oa"} )  #  ### dodavanje obelezavanja u selekcije
         if atom['type'][0]=="H" and  atom['label']==[]:
          Ha2=atom
          Ha2.update( {'label' : "H2a"} )  #  ### dodavanje obelezavanja u selekcije



        #### atom Hb1 je oni koji je blizi Ha1
        HH=100
        for atom in MOLb:
         if atom['type'][0]=="H":
          if np.linalg.norm(atom['coord']-Ha1['coord'])< HH:
           Hb1=atom
           HH= np.linalg.norm(atom['coord']-Ha1['coord'])
        Hb1.update( {'label' : "Hb1"} )  #  ### dodavanje obelezavanja u selekcije
        ##### odredjujem ostale atome u molekulu b
        
        for atom in MOLb:
         if atom['type'][0]=="H" and  atom['label']==[] :
          Hb2=atom
          Hb2.update( {'label' : "Hb2"} )  #  ### dodavanje obelezavanja u selekcije


        #print "COKA", Ha1
        #print "COKA", Ob
        #print "COKA", Oa
        #print "COKA", Ha2
        #print "COKA", Hb2
        #print "COKA", Hb1







        #### odredjivanje ostalih parametara:
        ### za sada imam OH i HH

        
        P1P2=  angle3points(planevector(Oa['coord'], Ha1['coord'],  Ha2['coord'])[0],[0,0,0],planevector(Ob['coord'], Hb1['coord'],  Hb2['coord'])[0]) 
        alpha=angle3points(Oa['coord'], Ha1['coord'],Ob['coord'])
        v1=angle3points(Oa['coord']- Ha1['coord'],[0,0,0], Ob['coord']-Hb1['coord'])
        v2=angle3points(Oa['coord']- Ha1['coord'],[0,0,0], Ob['coord']-Hb2['coord'])
        T=dihedral(  Ha2['coord'], Oa['coord'], Ha1['coord'], Ob['coord']   )
        OO= np.linalg.norm(Oa['coord']-Ob['coord'])

        print   ""

        print   "Geometrical parameters:"
        print "HOH:", water
        print "HOH:", watersec
        print "P1P2:",P1P2
        print "alpha:",alpha
        print "OH:",OH
        print "HH:",HH
        print "v1:",v1
        print "v2:",v2
        print "T:",T
        print "OO:",OO


        ############### odredjivanje vrste interakcije:
        print   ""

        print   "The interaction group is:"

        if OH-OO>=0 or HH<=OH-0.6 or HH<1.5 or OO<2.5 :##  za odbojne
         print "repulsive"
        else:
         if v1>160 and v1<=180 and alpha>80 and alpha <140 and T>40:    ##### paralelna I
          print "parallel I"
         if v2>160 and v2<=180 and alpha>80 and alpha <120 and T>40:    ##### paralelna I
          print "parallel II"
         if OH<HH :    ##### paralelna I
          print "classical"
         else:
          print "other"   #### ????
          




       if stop==True:   # samo jedan pogodak na 4A
         break
     
   if stop==False:
     print "OO>4"

else:
     print " water angles not in range 96.4-112.8"


