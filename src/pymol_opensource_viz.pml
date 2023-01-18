# pymol visualization script

bg_color white
set ray_trace_fog = 0
set depth_cue=0

######################################
# All Mutated Residues
######################################

reinitialize

fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 439
show sticks, chain B and resi 440
show sticks, chain B and resi 458
show sticks, chain B and resi 492
show sticks, chain B and resi 501

color brightorange, chain B and resi 439
color brightorange, chain B and resi 440
color brightorange, chain B and resi 458
color brightorange, chain B and resi 492
color brightorange, chain B and resi 501

set cartoon_side_chain_helper, on

set_view (\
     0.612492561,    0.004922007,   -0.790461779,\
    -0.785307229,    0.117967658,   -0.607760489,\
     0.090257555,    0.993001461,    0.076118134,\
     0.000702232,   -0.001529753, -903.013366699,\
   188.205261230,  215.084365845,  225.523193359,\
   244.930374146, 1561.282714844,  -20.000000000 )

png ~/pymol_sars_open/overall_spike_ACE2.png, width=1395, height=2040, dpi=300

set_view (\
     0.454311877,   -0.140507132,   -0.879695654,\
    -0.858511627,    0.194572613,   -0.474442273,\
     0.237826779,    0.970763862,   -0.032230958,\
     0.001055326,   -0.005033240, -178.611297607,\
   195.007263184,  212.024856567,  275.484619141,\
  -353.806640625,  712.439086914,  -20.000000000 )

set_view (\
     0.612492561,    0.004922007,   -0.790461779,\
    -0.785307229,    0.117967658,   -0.607760489,\
     0.090257555,    0.993001461,    0.076118134,\
     0.000695243,   -0.001880370, -172.155731201,\
   191.704284668,  215.823516846,  269.894714355,\
  -485.764892578,  830.587585449,  -20.000000000 )

png ~/pymol_sars_open/zoomin_spike_ACE2.png, width=1395, height=1000, dpi=300


######################################
# N501Y
######################################

reinitialize
fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 501
color brightorange, chain B and resi 501
set cartoon_side_chain_helper, on

set_view (\
    -0.188978806,    0.054364681,   -0.980465770,\
    -0.782804728,   -0.611160755,    0.116996348,\
    -0.592872083,    0.789635003,    0.158055186,\
     0.000000000,    0.000000000, -110.238113403,\
   161.081115723,  209.512115479,  275.955505371,\
    93.333839417,  127.142417908,  -20.000000000 )

png ~/pymol_sars_open/N501Y_pre.png, dpi=300

wizard mutagenesis

set_view (\
    -0.188978806,    0.054364681,   -0.980465770,\
    -0.782804728,   -0.611160755,    0.116996348,\
    -0.592872083,    0.789635003,    0.158055186,\
     0.000000000,    0.000000000, -110.238113403,\
   161.081115723,  209.512115479,  275.955505371,\
    93.333839417,  127.142417908,  -20.000000000 )

png ~/pymol_sars_open/N501Y_post.png, dpi=300

######################################
# K458T
######################################

reinitialize
fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 458
color brightorange, chain B and resi 458
set cartoon_side_chain_helper, on

set_view (\
    -0.113049634,    0.130972430,    0.984920323,\
     0.750862598,    0.660454810,   -0.001642140,\
    -0.650708914,    0.739353299,   -0.173005939,\
     0.000000000,    0.000000000, -138.429046631,\
   181.927658081,  186.251770020,  276.408813477,\
   107.911582947,  168.946533203,  -20.000000000 )

png ~/pymol_sars_open/K458T_pre.png, dpi=300

wizard mutagenesis

set_view (\
    -0.113049634,    0.130972430,    0.984920323,\
     0.750862598,    0.660454810,   -0.001642140,\
    -0.650708914,    0.739353299,   -0.173005939,\
     0.000000000,    0.000000000, -138.429046631,\
   181.927658081,  186.251770020,  276.408813477,\
   107.911582947,  168.946533203,  -20.000000000 )

png ~/pymol_sars_open/K458T_post.png, dpi=300


######################################
# L492I
######################################

reinitialize
fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 492
color brightorange, chain B and resi 492
set cartoon_side_chain_helper, on

set_view (\
     0.858329177,   -0.005835568,   -0.513066828,\
    -0.511126161,    0.077858686,   -0.855969310,\
     0.044941403,    0.996945083,    0.063846298,\
     0.000258982,   -0.000288218, -110.723846436,\
   179.743041992,  209.037048340,  275.132873535,\
    80.686882019,  140.867202759,  -20.000000000 )

png ~/pymol_sars_open/L492I_pre.png, dpi=300

wizard mutagenesis

set_view (\
     0.858329177,   -0.005835568,   -0.513066828,\
    -0.511126161,    0.077858686,   -0.855969310,\
     0.044941403,    0.996945083,    0.063846298,\
     0.000258982,   -0.000288218, -110.723846436,\
   179.743041992,  209.037048340,  275.132873535,\
    80.686882019,  140.867202759,  -20.000000000 )

png ~/pymol_sars_open/L492I_post.png, dpi=300

######################################
# N439A
######################################

reinitialize
fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 439
color brightorange, chain B and resi 439
set cartoon_side_chain_helper, on

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000756673,   -0.000629012,  -90.180236816,\
   160.862762451,  214.211868286,  269.007904053,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars_open/L439A_pre.png, dpi=300

wizard mutagenesis

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000756673,   -0.000629012,  -90.180236816,\
   160.862762451,  214.211868286,  269.007904053,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars_open/L439A_post.png, dpi=300

######################################
# N440G
######################################

reinitialize
fetch 7df4
show cartoon
hide lines

color slate, chain A
color firebrick, chain B or chain C or chain D

show sticks, chain B and resi 440
color brightorange, chain B and resi 440
set cartoon_side_chain_helper, on

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000738518,   -0.000597056,  -90.316703796,\
   159.189788818,  215.771575928,  267.230072021,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars_open/N440G_pre.png, dpi=300

wizard mutagenesis

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000738518,   -0.000597056,  -90.316703796,\
   159.189788818,  215.771575928,  267.230072021,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars_open/N440G_post.png, dpi=300


