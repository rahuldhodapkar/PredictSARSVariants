# pymol visualization script

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
    -0.085034601,   -0.012658661,   -0.996291578,\
    -0.882019579,   -0.464156508,    0.081181020,\
    -0.463468999,    0.885661840,    0.028303608,\
     0.001862928,    0.001709359,  -98.103118896,\
   172.093582153,  209.744323730,  276.604003906,\
    77.050148010,  120.048622131,  -20.000000000 )

png ~/pymol_sars_open/N501_pre.png, dpi=300

wizard mutagenesis

set_view (\
    -0.085034601,   -0.012658661,   -0.996291578,\
    -0.882019579,   -0.464156508,    0.081181020,\
    -0.463468999,    0.885661840,    0.028303608,\
     0.001862928,    0.001709359,  -98.103118896,\
   172.093582153,  209.744323730,  276.604003906,\
    77.050148010,  120.048622131,  -20.000000000 )

png ~/pymol_sars_open/N501_post.png, dpi=300


sele chain B and resi 501 expand 5
# manually expand selection by 5 angstroms

# run "clean" simulation to test affinity changes

# Energy = 173.57

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

# manually zoom in on structure and use the mutagenesis wizard

sele chain B and resi 458 expand 5
# manually expand selection by 5 angstroms

# run "clean" simulation to test affinity changes

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
     0.809001446,   -0.028291758,   -0.587126017,\
    -0.587715030,   -0.021358497,   -0.808784068,\
     0.010341465,    0.999370039,   -0.033906490,\
     0.000258982,   -0.000288218, -110.723846436,\
   179.743041992,  209.037048340,  275.132873535,\
   -28.253112793,  249.807235718,  -20.000000000 )

# manually zoom in on structure and use the mutagenesis wizard

sele chain B and resi 492 expand 5
clean sele

# Energy = 28.66

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

png ~/pymol_sars/n439a_pre.png, 0, 0, -1, ray=0

sele chain B and resi 439 expand 5
clean sele

# Energy = 19.73

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000756673,   -0.000629012,  -90.180236816,\
   160.862762451,  214.211868286,  269.007904053,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars/n439a_post.png, 0, 0, -1, ray=0

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

png ~/pymol_sars/n440g_pre.png, 0, 0, -1, ray=0

sele chain B and resi 440 expand 5
clean sele

# Energy = -81.79

set_view (\
     0.674618959,    0.223594114,   -0.703485966,\
    -0.712679744,   -0.050963290,   -0.699636936,\
    -0.192287996,    0.973347127,    0.124968916,\
    -0.000738518,   -0.000597056,  -90.316703796,\
   159.189788818,  215.771575928,  267.230072021,\
    57.988162994,  121.866149902,  -20.000000000 )

png ~/pymol_sars/n440g_post.png, 0, 0, -1, ray=0


