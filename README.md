# PhDThesis
Versions 1.0 and 2.0 of this program (with all their functions) is implemented by Cristina Gómez Santamaría, Universidad Pontificia Bolivariana, Colombia.
Macro and Microcell scenarios are supported, with manually configurable arbitrary array antenna configuration. For several base stations please carefully configure a cuadratic number of base stations, i.e., 4 (obtaining a 2x2 matrix), 9 (3x3 matrix), 16 (4x4 matrix)... just for easy topology configuration and generation of random scatterers.
This implementation is based on the following publications:
- “A generic model for MIMO wireless propagation channels in Macro- and Microcells”, A.Molisch, IEEE Transactions on Signal Processing, Vol.52, No.1, Jan 2004.
- “The COST259 directional channel model - Part I: overview and methodology”, A.Molisch et al, IEEE Transactions on Wireless Communications, Vol.5, No.12, Dec 2006.
- “The COST259 directional channel model - Part II: Macrocells”, H.Asplund et al, IEEE Transactions on Wireless Communications, Vol.5, No.12, Dec 2006.
- “Wireless flexible personalized communications. COST259: European Co-operation in Mobile Radio Research”, L.Correia, ISBN 978-0-471-49836-0.

Two scenarios are implemented so far
SemiUrban_300MHz: both NLOS and LOS single link MIMO simulations 
are supported. For outdoor LOS multiple link MIMO simulation is supported as well. 
IndoorHall_5GHz: OLOS single link MIMO simulation is supported.

If you use the COST 2100 channel model for publications, please refer to: 
L. Liu, J. Poutanen, F. Quitin, K. Haneda, F. Tufvesson, P. De Doncker,
P. Vainikainen and C. Oestges, ìThe COST 2100 MIMO channel model,î
IEEE Wireless Commun., vol 19, issue 6, pp 92-99, Dec. 2012.

Further details about the COST 2100 channel model can be found in:
Roberto Verdone (Editor), Alberto Zanella (Editor)
Pervasive Mobile and Ambient Wireless Communications Pervasive Mobile 
and Ambient Wireless Communications, ISBN 978-1-4471-2315-6,
Springer, 2012. 

If you use the SemiUrban_300MHz scenario, 
further information can be found in:
1. Meifang Zhu, Gunnar Eriksson, and Fredrik Tufvesson, 
"The COST 2100 Channel Model: Parameterization and Validation 
Based on Outdoor MIMO Measurements at 300 MHz", 
IEEE Transactions on Wireless Commun..
2. Meifang Zhu and Fredrik Tufvesson, "Virtual Multi-link Propagation 
Investigation of an Outdoor Scenario At 300 MHz," Proceedings of the 
7th European Conference on Antennas and Propagation (EUCAP), Gothenburg,
Sweden, April 2013.

If you use the IndoorHall_5GHz scenario, 
further information can be found in:
1. V. M. Kolmonen, P. Almers, J. Salmi, J. Koivunen, A. Richter,
F. Tufvesson, A. Molisch, P. Vainikainen, "A dynamic dual-link 
wideband MIMO channel sounder for 5.3 GHz," IEEE Transactions on 
Instrumentation and Measurement, Vol. 59, No. 4, pp. 873-883, 2010.
2. J. Poutanen, K. Haneda, L. Lin, C. Oestges, F. Tufvesson , 
P. Vainikainen, "Parameterization of the COST 2100 MIMO channel 
modeling in indoor scenarios," Proceedings of the 5th European 
Conference on Antennas and Propagation (EUCAP), Rome, Italy, 
pp. 3606-3610, April 2011.



The script "demo_model.m" provides an example for testing the COST 2100 model. It selects scenario, link type, and initiates the simulation.
The output of demo_model.m is channel functions, with size dependent on the parameter choice and setup: 

1) SISO_omni: Transfer function for SISO omni-directional antenna
create_IR_omni: users have to set up the frequency separation, delta_f
 
2) MIMO_omni: Transfer function for MIMO omini-directional antenna
 create_IR_omni_MIMO: users have to set up the frequency separation, delta_f.
 Only 2 by 2 MIMO system is implemented.
 
3) MIMO_dipole: Transfer function for a theoretical antenna response for 
 any size of lambda/2-spaced linear dipole antenna arrays. An Ntx-by-Nrx theoretical 
 antenna array response is generated and the correponding 
 channel transfer function is simulated.
 
4) MIMO_measured: Transfer function for any measured MIMO antenna response
 get_H: users have to provide the full antenna response at the BS and 
 MS sides, and also the rotation of the antenna arrays. The antenna 
 response mat file have to be the same format as 'antSample.mat' file.


The modifications are listed in document changelog.

Meifang Zhu, 2013.02.06


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 The code is developed under GPL. The COST2100 code was originally written by Lingfeng Liu, UniverstiÈ catholique de Louvain (UCL). 

email: lingfeng.liu@uclouvain.be
tel: +32 (0) 10 47 81 05
address: B‚timent Maxwell, Place du Levent 3, 1348 Louvain-la-Neuve, Belgium


1,Function list
calc_dist
calc_pathloss
cost2100
demo_model
draw_circ
draw_ellpsoid
get_channel
get_channel_los
get_cluster_local
get_common
get_cluster
get_dmc
get_H
get_IR
get_mpc
get_para
get_VR
get_VRLOS
get_VRtable
rotate_matrix
setFontsize
update_chan
visual_pddp
visual_channel

2,Other files list
unit.txt: list of parameter unit defined in the model
antSample.mat: antenna array radiation pattern sample file

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
3, How to run the code

The main routine of the code is cost2100

The visualization routines are named as 'visual_*'

To understand each function,  type help function_name in matlab

The test script is demo_model.m
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
4, How to customize the channel implementation

Please read get_para for more instructions

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
5, Copyright
Copyright (C)2008 LIU Ling-Feng, UniversitÈ catholique de Louvain, Belgium
This program, cost2100, is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at your 
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

6, Acknowledgements
Thanks for Helmut Hofstetter for his previous code on the COST 273 
channel model which inspires me a lot.

Thanks for Claude Oestges and Nicolai Czink for their feedbacks and theory 
supports.

Thanks for Katsuyuki Haneda, Juho Puotanen, Fredrik Tufvesson, and Meifang 
Zhu for their active participations and code testing.

Also thanks for my wife Qin, always be patient during my coding time.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
7, Some further references

[1] L.M. Correia, Mobile Broadband Multimedia Networks. Academic Press, 2006.
[2] N. Czink and C. Oestges, ìThe COST 273 MIMO channel model: Three kinds of clusters,î IEEE 10th Int. Sym., ISSSTAí08, pp. 282ñ286, 2008.
[3] L. Liu, N. Czink, and C. Oestges, Implementing the COST 273 MIMO channel model, NEWCOM 2009

