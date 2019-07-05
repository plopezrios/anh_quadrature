ANH_QUADRATURE
==============
Tool to obtain the ground state of a particle in a one-dimensional
potential expressed as a polynomial, and to construct Gauss-Hermite
quadrature grids that exactly evaluate expectation values of local
operators that can be expressed as a finite-order polynomial.  Plots
of the potential and eigenfunctions, of the quadrature grids, of the
local operators, and of their expectation values as a function of grid
size are produced.  Note that the local operators used in the
convergence tests are hard-coded.

This tool implementes the grid construction method presented in Phys.
Rev. B 97, 054104 (2018) (https://doi.org/10.1103/PhysRevB.97.054104;
see https://arxiv.org/abs/1711.06265 for preprint) and was used to
make the plots in the paper.  The one-dimensional Schroedinger
equation solver in ANH_QUADRATURE is also used in VIB_LEVELS.

ANH_QUADRATURE is a stand-alone Fortran utility which uses a LAPACK
library.  The code can be compiled with:

```
gfortran -o anh_quadrature anh_quadrature.f90 -llapack
```

The compile.sh script is provided for reference only.

The code asks for the polynomial coefficients defining the potential
(in a.u., or whatever rescaled units one wishes to use) on standard
input.  Optionally, one can provide an effective mass in Dalton
[=a.m.u.; 1 a.u. by default], and a centre and frequency for the
eigensolver [autodetected by default].

Usage example
=============
The following is output from a run which corresponds to one of the
examples in Phys. Rev. B 97, 054104 (2018).  Besides the standard
output, plot files 1D_schroedinger.agr, grid.agr, functions.agr, and
quadrature_int.agr are produced which can be opened with xmgrace.

```
./anh_quadrature 
 Enter coefficients c_k of expansion of V(u) in natural powers,

   V(u) = sum_k=0^n c_k u^k       (a.u., in one line):

0 0 0.5 0.25 0.125

 Effective mass in Da [e.g., 0.5 for H2; empty for 1 a.u.]:


 Enter ucentre and omega (a.u.) [single line; empty to autodetect]:

 Centre          =   0.20491302987311194     
 Omega (a.u.)    =    1.2100000000000004     
 Expansion order = 20
 E0 (a.u.)       =   0.53736996858686670     
 Virial ratio    =    1.0000000004735652     

 Analytical grids
 ================
 Symmetric 2-point grid:
   u_1, P_1 =  -0.86331842188252750       0.50000000000000000     
   u_2, P_2 =   0.45349236213630367       0.50000000000000000     
 Symmetric 3-point grid:
   u_1, P_1 =   -1.3174038515493862       0.17513134732760280     
   u_2, P_2 =  -0.20491302987311194       0.64973730534479435     
   u_3, P_3 =   0.90757779180316223       0.17513134732760280     
 Symmetric 4-point grid:
   u_1, P_1 =   -1.6728743419302781        5.2627927856569667E-002
   u_2, P_2 =  -0.68553160457016371       0.44737207214343033     
   u_3, P_3 =   0.27570554482393983       0.44737207214343033     
   u_4, P_4 =    1.2630482821840543        5.2627927856569667E-002
 Non-symmetric 2-point grid:
   u_1, P_1 =  -0.89678344922915088       0.47523136824145723     
   u_2, P_2 =   0.42164600212617820       0.52476863175854283     

 Uniform-grid integration of test functions
 ==========================================
 Interval: [  -6.5358400984100440      :   6.0293842890189335      ]
 Integral values and variances and symmetrized variances:
   <f_1> =   0.47548701007905492       0.44460741820315736       0.44460741820315736     
   <f_2> =    9.2751789892397712E-002   2.1631804058894679        2.4258580123392108     
   <f_3> =   0.62038961856706554        25.066888410562925        25.066888410562925     
   <f_4> =   0.64989448256010141       0.42372322048232913       0.42372322048232913     

 Numerical grids
 ===============
 Grid size 2:
   u_1, P_1 =  -0.89678344786635200       0.47523136993041748     
   u_2, P_2 =   0.42164600314251671       0.52476863006958241     
   Integration tests:
     <f_1> =   0.47548701042526337        9.7864520627241797E-002   9.7864520627241797E-002
     <f_2> =  -0.25399110822338516       0.11211236369097816        1.0529024876551803E-002
     <f_3> =   0.29404595495925162        5.0367828155015897E-004   5.0367828155015897E-004
     <f_4> =   0.74969977382785646       0.14622686719475744       0.14622686719475744     

 Grid size 3:
   u_1, P_1 =   -1.3816507786449372       0.15627209827372848     
   u_2, P_2 =  -0.24020298812204380       0.64621042257882511     
   u_3, P_3 =   0.84155646105679505       0.19751747914744633     
   Integration tests:
     <f_1> =   0.47548701007557637       0.44460741820111160       0.44460741820111160     
     <f_2> =    9.2751789890529873E-002  0.19056178466032378       0.45323939111930012     
     <f_3> =   0.23041793027656890        3.1770303545864922E-002   3.1770303545864922E-002
     <f_4> =   0.62636374008929363       0.48691882931728636       0.48691882931728636     

 Grid size 4:
   u_1, P_1 =   -1.7559330430052289        4.5020732177790714E-002
   u_2, P_2 =  -0.73779061253581235       0.41501378730789001     
   u_3, P_3 =   0.22122439556554252       0.47540065094158596     
   u_4, P_4 =    1.1641541760442802        6.4564829572733268E-002
   Integration tests:
     <f_1> =   0.47548701007908034       0.44460741820378386       0.44460741820378386     
     <f_2> =    9.2751789893014761E-002   1.6088253902719434        1.8715029967232437     
     <f_3> =   0.62038961859275554        3.2266123837644574        3.2266123837644574     
     <f_4> =   0.65775944312325185       0.40482826317151605       0.40482826317151605     

 Grid size 5:
   u_1, P_1 =   -2.0612482285260212        1.2178207337137913E-002
   u_2, P_2 =   -1.1384320241462429       0.20071754178059287     
   u_3, P_3 =  -0.26295488375462228       0.50361168372590581     
   u_4, P_4 =   0.58054268631029060       0.26400299584032677     
   u_5, P_5 =    1.4292104522638269        1.9489571316036760E-002
   Integration tests:
     <f_1> =   0.47548701009105548       0.44460741816686117       0.44460741816686117     
     <f_2> =    9.2751789846775720E-002   2.1631804065267719        2.4258580129066796     
     <f_3> =   0.62038961849533758        14.303020488807311        14.303020488807311     
     <f_4> =   0.65518890227318272       0.41381450263589731       0.41381450263589731     

 Grid size 6:
   u_1, P_1 =   -2.3198058543611739        3.1761491160420237E-003
   u_2, P_2 =   -1.4722697978112791        8.3372641135014650E-002
   u_3, P_3 =  -0.66637613570715082       0.36297362293701857     
   u_4, P_4 =   0.11615045215478108       0.42052657630362489     
   u_5, P_5 =   0.87633178894559438       0.12436561365029650     
   u_6, P_6 =    1.6559447120041784        5.5853968580034030E-003
   Integration tests:
     <f_1> =   0.47548701007907035       0.44460741820378513       0.44460741820378513     
     <f_2> =    9.2751789893024114E-002   2.1631804070041327        2.4258580134554220     
     <f_3> =   0.62038961859272612        23.223967446066094        23.223967446066094     
     <f_4> =   0.64316894443670869       0.43768735907759654       0.43768735907759654     

 Grid size 7:
   u_1, P_1 =   -2.5447363816014610        8.0849877992583002E-004
   u_2, P_2 =   -1.7575813427693188        3.1522728733606384E-002
   u_3, P_3 =   -1.0111913852566679       0.21054893534808258     
   u_4, P_4 =  -0.27954296520610183       0.42212485876513600     
   u_5, P_5 =   0.43186506649915241       0.28095360240284128     
   u_6, P_6 =    1.1287802185822493        5.2498736273247144E-002
   u_7, P_7 =    1.8551767491794859        1.5426396971608785E-003
   Integration tests:
     <f_1> =   0.47548701007906635       0.44460741820380045       0.44460741820380045     
     <f_2> =    9.2751789893034980E-002   2.1631804070039697        2.4258580134552372     
     <f_3> =   0.62038961859263952        25.066890334742432        25.066890334742432     
     <f_4> =   0.65836836111483765       0.40650086522525780       0.40650086522525780     

 Grid size 8:
   u_1, P_1 =   -2.7444405265415099        2.0209975349326704E-004
   u_2, P_2 =   -2.0065315970014757        1.1169447114958795E-002
   u_3, P_3 =   -1.3108339544655760       0.10648017256516737     
   u_4, P_4 =  -0.62503361192880258       0.32366729137144784     
   u_5, P_5 =    4.7864515172000112E-002  0.37610622876656141     
   u_6, P_6 =   0.70256050310064266       0.16148576812693333     
   u_7, P_7 =    1.3497623357106303        2.0474589820811830E-002
   u_8, P_8 =    2.0336177331041267        4.1440248062605124E-004
   Integration tests:
     <f_1> =   0.47548701007921379       0.44460741820359428       0.44460741820359428     
     <f_2> =    9.2751789892678918E-002   2.1631804070047114        2.4258580134565579     
     <f_3> =   0.62038961859348096        25.066890334735106        25.066890334735106     
     <f_4> =   0.64107701252288285       0.44097201137386133       0.44097201137386133     

 Grid size 9:
   u_1, P_1 =   -2.9251364365263774        4.9556575001350425E-005
   u_2, P_2 =   -2.2279908954038361        3.7606441677641195E-003
   u_3, P_3 =   -1.5753101000619349        4.8910414281020831E-002
   u_4, P_4 =  -0.93087813598324964       0.20863712601526019     
   u_5, P_5 =  -0.29293995966757619       0.36799289758215775     
   u_6, P_6 =   0.33095604703537596       0.27979063134981719     
   u_7, P_7 =   0.93961658418365346        8.3226399236427911E-002
   u_8, P_8 =    1.5467368861095707        7.5233381436294721E-003
   u_9, P_9 =    2.1956798626605658        1.0899264892113745E-004
   Integration tests:
     <f_1> =   0.47548701007906868       0.44460741820378119       0.44460741820378119     
     <f_2> =    9.2751789893016370E-002   2.1631804070040910        2.4258580134553722     
     <f_3> =   0.62038961859272712        25.066890334741796        25.066890334741796     
     <f_4> =   0.65707455562795303       0.40951134631943120       0.40951134631943120     

 Grid size 10:
   u_1, P_1 =   -3.0963896621537117        1.1373588896857764E-005
   u_2, P_2 =   -2.4326485547795746        1.1821408768893813E-003
   u_3, P_3 =   -1.8159455235592759        2.0586035020221292E-002
   u_4, P_4 =   -1.2078291274381250       0.11834503702439671     
   u_5, P_5 =  -0.60209152430014290       0.29272431627970996     
   u_6, P_6 =   -4.6751628680178616E-003  0.34166543301687197     
   u_7, P_7 =   0.57844464221497516       0.18316088614120410     
   u_8, P_8 =    1.1497205862064073        3.9647700653712155E-002
   u_9, P_9 =    1.7242985286334800        2.6488572232541050E-003
   u_10, P_10 =    2.3444634611175181        2.8220174843388717E-005
   Integration tests:
     <f_1> =   0.47548701007907002       0.44460741820378502       0.44460741820378502     
     <f_2> =    9.2751789893024406E-002   2.1631804070041354        2.4258580134554233     
     <f_3> =   0.62038961859272601        25.066890334741373        25.066890334741373     
     <f_4> =   0.64262480010956213       0.43736180662969937       0.43736180662969937     
```
